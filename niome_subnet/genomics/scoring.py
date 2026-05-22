#!/usr/bin/env python3
import bittensor as bt
import json
import pysam
import os
from pathlib import Path

from collections import defaultdict
from niome_subnet.genomics.model import GroundTruth, MinerScore, MinerSubmission
from niome_subnet.genomics.vcf_norm import normalize_vcf, preprocess_vcf


# -----------------------------
# Create mapping file
# -----------------------------
def create_mapping_file(ref_fasta: str, read1: str, read2: str) -> str:
    bam_path = "data/sim.bam"
    if os.path.exists(bam_path):
        return bam_path

    # bwa index
    subprocess.run(f"bwa index {ref_fasta}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # align
    subprocess.run(f"bwa mem {ref_fasta} {read1} {read2} | samtools sort -o {bam_path}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(f"samtools index {bam_path}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return bam_path


# -----------------------------
# 1. VCF LOADER
# -----------------------------
def normalize_gt(gt_tuple):
    """Sort GT alleles to make phase-insensitive comparison (0/1 == 1/0)."""
    if gt_tuple is None:
        return None
    return tuple(sorted(a if a is not None else -1 for a in gt_tuple))


def load_vcf(path):
    vcf = pysam.VariantFile(path)
    variants = {}  # (contig, pos, ref, alt) -> normalized GT tuple or None

    for rec in vcf.fetch():
        if rec.alts is None:
            continue
        for alt in rec.alts:
            key = (rec.contig, rec.pos, rec.ref, alt)
            gt = None
            try:
                sample = next(iter(rec.samples.values()))
                gt = normalize_gt(sample['GT'])
            except (StopIteration, KeyError, TypeError):
                pass
            variants[key] = gt

    return variants


def _resolve_scoring_reference(ref: str) -> Path:
    """Resolve validator ref path (often data/ref.fa) to an indexed FASTA."""
    ref_path = Path(ref)
    if not ref_path.is_absolute():
        ref_path = Path.cwd() / ref_path
    return ref_path.resolve()


# -----------------------------
# 2. LOAD DEPTH FROM BAM
# -----------------------------
def load_depth(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    depth = defaultdict(int)

    for pileup_col in bam.pileup():
        depth[pileup_col.reference_pos + 1] = pileup_col.nsegments

    return depth


# -----------------------------
# 4. COMPUTE WEIGHTED SETS
# -----------------------------
def compute_weighted_sets(truth: dict, pred: dict, depth_map: dict):
    """Compute weighted TP/FP/FN with genotype-aware matching.

    Per-variant match score:
      1.0 — correct REF+ALT+genotype
      0.5 — correct REF+ALT, wrong zygosity (phase-insensitive: 0/1 == 1/0)
      0.0 — missed
    """
    tp_w = 0.0
    fp_w = 0.0
    fn_w = 0.0

    for key, truth_gt in truth.items():
        w = variant_weight(depth_map.get(key[1], 0))
        if key in pred:
            pred_gt = pred[key]
            if truth_gt is None or truth_gt == pred_gt:
                match = 1.0
            else:
                match = 0.5
            tp_w += w * match
            fn_w += w * (1.0 - match)
        else:
            fn_w += w

    for key in pred:
        if key not in truth:
            fp_w += variant_weight(depth_map.get(key[1], 0))

    return tp_w, fp_w, fn_w


# -----------------------------
# 5. VARIANT DIFFICULTY MODEL
# -----------------------------
def variant_weight(depth):
    # Simple but effective heuristic
    if depth < 10:
        return 0.3   # hard
    elif depth < 20:
        return 0.6   # medium
    else:
        return 1.0   # easy


# -----------------------------
# 6. WEIGHTED METRICS
# -----------------------------
def weighted_metrics(tp_w: float, fp_w: float, fn_w: float):
    precision = tp_w / (tp_w + fp_w + 1e-9)
    recall = tp_w / (tp_w + fn_w + 1e-9)
    f1 = 2 * precision * recall / (precision + recall + 1e-9)
    return precision, recall, f1


# -----------------------------
# 7. FINAL SCORE
# -----------------------------
def score_vcf(p, r, f1, response_time):
    score1 = 0.4 * f1 + 0.3 * p + 0.3 * r
    return score1


# -----------------------------
# 9. ANNOTATION SCORING
# -----------------------------
_CFTR_DRUGS = [
    "ivacaftor",
    "tezacaftor_ivacaftor",
    "elexacaftor_tezacaftor_ivacaftor",
    "lumacaftor_ivacaftor",
]


def score_annotations(miner_annotations: dict, truth_annotations: dict) -> float:
    """Score miner CFTR2 annotations against ground truth.

    Iterates over truth entries only; missing miner entries score 0.

    Per-variant match score:
      - hgvs match:                  weight 0.1
      - clinical_significance match: weight 0.1
      - drug_response (4 drugs):     weight 0.8 total (0.2 each drug)
    """
    if not truth_annotations:
        return 0.0

    total_score = 0.0
    for rsid, truth_entry in truth_annotations.items():
        miner_entry = miner_annotations.get(rsid, {})
        variant_score = 0.0

        if miner_entry.get("hgvs") == truth_entry.get("hgvs"):
            variant_score += 0.1

        if miner_entry.get("clinical_significance") == truth_entry.get("clinical_significance"):
            variant_score += 0.1

        truth_drug = truth_entry.get("drug_response", {})
        miner_drug = miner_entry.get("drug_response", {})
        for drug in _CFTR_DRUGS:
            if miner_drug.get(drug) == truth_drug.get(drug):
                variant_score += 0.2

        total_score += variant_score

    # Penalise FP annotations: divide by the larger of truth or miner count so
    # submitting the entire cftr2_annotations.json does not yield a free high score.
    denominator = max(len(truth_annotations), len(miner_annotations))
    return total_score / denominator


# -----------------------------
# 8. MAIN PIPELINE
# -----------------------------
def score(miner_submission: MinerSubmission, ground_truth: GroundTruth, bam: str) -> MinerScore:
    try:
        miner_origin_vcf = "data/miner.vcf"
        with open(miner_origin_vcf, "w") as f:
            f.write(miner_submission.vcf_content)

        miner_vcf = preprocess_vcf(miner_origin_vcf)

        # [1] Normalizing VCFs
        truth_norm = "data/truth.norm.vcf.gz"
        miner_norm = "data/miner.norm.vcf.gz"

        ref_path = _resolve_scoring_reference(ground_truth.ref)
        normalize_vcf(ground_truth.truth_vcf, ref_path, truth_norm)
        try:
            normalize_vcf(miner_vcf, ref_path, miner_norm)
        except RuntimeError as norm_exc:
            raise RuntimeError(
                f"Miner VCF failed bcftools norm (validator would score 0): {norm_exc}"
            ) from norm_exc

        # [2] Loading variants
        truth = load_vcf(truth_norm)
        pred  = load_vcf(miner_norm)

        # [3] Loading depth from BAM
        depth = load_depth(bam)

        # [4] Computing weighted sets (genotype-aware)
        tp_w, fp_w, fn_w = compute_weighted_sets(truth, pred, depth)

        # [5] Computing weighted metrics
        p, r, f1 = weighted_metrics(tp_w, fp_w, fn_w)

        vcf_score = score_vcf(p, r, f1, miner_submission.response_time)

        # Score CFTR2 annotations (weight 0.3 in final score)
        annotation_score = 0.0
        if miner_submission.cftr_annotations is not None and ground_truth.cftr2_annotations:
            try:
                with open(ground_truth.cftr2_annotations) as _f:
                    truth_annotations = json.load(_f)
                annotation_score = score_annotations(miner_submission.cftr_annotations, truth_annotations)
            except Exception as _e:
                bt.logging.warning(f"Failed to score annotations for miner {miner_submission.uid}: {_e}")

        score_val = 0.7 * vcf_score + 0.3 * annotation_score

        miner_score = MinerScore(
            uid=miner_submission.uid,
            precision=p,
            recall=r,
            f1_score=f1,
            response_time=miner_submission.response_time,
            vcf_score=vcf_score,
            annotation_score=annotation_score,
            final_score=score_val,
            log=f"VCF Score: {vcf_score:.4f}, Annotation Score: {annotation_score:.4f}, Final Score: {score_val:.4f}\n\nMiner VCF\n{miner_submission.vcf_content}",
        )

        return miner_score
    except Exception as e:
        bt.logging.warning(f"Error scoring miner {miner_submission.uid}: {e}")
        return MinerScore(
            uid=miner_submission.uid,
            precision=0.0,
            recall=0.0,
            f1_score=0.0,
            response_time=miner_submission.response_time,
            vcf_score=0.0,
            annotation_score=0.0,
            final_score=0.0,
            log=f"Error: {e}",
        )
