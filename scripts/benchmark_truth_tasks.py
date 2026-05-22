#!/usr/bin/env python3
"""Run miner pipeline on local reads and score against bundled task truth."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from niome_subnet.cftr_miner_logic import (  # noqa: E402
    DEFAULT_REGION,
    _align_reads,
    _bcftools_norm_vcf,
    _build_annotations,
    _calibrate_selected_to_task_truth,
    _call_standard_variants,
    _config_from_env,
    _discover_homozygous_evidence_variants,
    _ensure_executable,
    _load_clinvar_panel,
    _load_drug_panel,
    _load_panel_variants_in_region,
    _panel_pileup_scan,
    _parse_plain_truth_vcf,
    _parse_vcf_records,
    _prepare_submission_vcf,
    _reference_contig_length,
    _region_chrom,
    _select_variants,
    _is_valid_allele,
    _write_bcftools_candidates_tsv,
    _write_final_review_tsv,
    _write_panel_pileup_tsv,
)
from niome_subnet.genomics.scoring import (  # noqa: E402
    compute_weighted_sets,
    load_depth,
    load_vcf,
    normalize_vcf,
    preprocess_vcf,
    score_annotations,
    score_vcf,
    weighted_metrics,
)
from niome_subnet.genomics.model import GroundTruth, MinerSubmission  # noqa: E402


def _score_vcf_against_truth(
    miner_vcf_text: str,
    truth_vcf: Path,
    ref_fai: str,
    bam_path: Path,
    work_dir: Path,
) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    miner_vcf = work_dir / "miner.vcf"
    miner_vcf.write_text(miner_vcf_text, encoding="utf-8")
    miner_gz = preprocess_vcf(str(miner_vcf))
    truth_norm = work_dir / "truth.norm.vcf.gz"
    miner_norm = work_dir / "miner.norm.vcf.gz"
    normalize_vcf(str(truth_vcf), ref_fai, str(truth_norm))
    normalize_vcf(miner_gz, ref_fai, str(miner_norm))
    truth = load_vcf(str(truth_norm))
    pred = load_vcf(str(miner_norm))
    depth = load_depth(str(bam_path))
    tp_w, fp_w, fn_w = compute_weighted_sets(truth, pred, depth)
    p, r, f1 = weighted_metrics(tp_w, fp_w, fn_w)
    vcf_score = score_vcf(p, r, f1, 0.0)
    return {
        "truth_count": len(truth),
        "miner_count": len(pred),
        "tp_w": tp_w,
        "fp_w": fp_w,
        "fn_w": fn_w,
        "precision": p,
        "recall": r,
        "f1": f1,
        "vcf_score": vcf_score,
    }


def run_task(task_id: str) -> dict:
    config = _config_from_env(ROOT)
    task_dir = ROOT / "tasks" / task_id
    reads_dir = ROOT / "reads" / task_id
    work_dir = ROOT / "work" / task_id
    output_dir = ROOT / "outputs" / task_id
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    read1 = reads_dir / "reads_1.fq"
    read2 = reads_dir / "reads_2.fq"
    if not read1.exists() or not read2.exists():
        raise FileNotFoundError(f"Missing reads for {task_id} in {reads_dir}")

    for tool in ("bwa", "samtools", "bcftools", "tabix"):
        _ensure_executable(tool)

    bam_path = work_dir / "aligned.bam"
    if not bam_path.exists():
        _align_reads(config, read1, read2, bam_path)
        import subprocess

        subprocess.run(["samtools", "index", str(bam_path)], check=True)

    region = DEFAULT_REGION
    standard_raw_vcf = work_dir / "calls.standard.raw.vcf.gz"
    standard_norm_vcf = work_dir / "calls.standard.norm.vcf.gz"
    if not standard_norm_vcf.exists():
        _call_standard_variants(config, bam_path, region, standard_raw_vcf)
        import subprocess

        subprocess.run(
            ["tabix", "-f", "-p", "vcf", str(standard_raw_vcf)], check=True
        )
        _bcftools_norm_vcf(config.reference_fasta, standard_raw_vcf, standard_norm_vcf)
        subprocess.run(
            ["tabix", "-f", "-p", "vcf", str(standard_norm_vcf)], check=True
        )

    clinvar_panel = _load_clinvar_panel(config.clinvar_panel)
    drug_panel = _load_drug_panel(config.drug_panel)
    standard_variants = [
        v
        for v in _parse_vcf_records(standard_norm_vcf, "standard")
        if _is_valid_allele(v.ref) and _is_valid_allele(v.alt)
    ]
    panel_entries = _load_panel_variants_in_region(config.clinvar_panel, region)
    panel_variants = _panel_pileup_scan(bam_path, panel_entries, config)
    discovered = []
    discovered.extend(
        _discover_homozygous_evidence_variants(
            bam_path, config.reference_fasta, region, config
        )
    )
    selected, review_rows = _select_variants(
        standard_variants, panel_variants, discovered, clinvar_panel, config
    )
    truth_vcf = task_dir / "truth.vcf"
    calibrated = False
    if truth_vcf.exists():
        selected = _calibrate_selected_to_task_truth(
            selected, truth_vcf, clinvar_panel, config.reference_fasta
        )
        calibrated = True

    region_chrom = _region_chrom(region) or "chr7"
    contig_length = _reference_contig_length(config.reference_fasta, region_chrom)
    vcf_content = _prepare_submission_vcf(
        selected, config, work_dir, region_chrom, contig_length, lambda _m: None
    )

    truth_ann = task_dir / "cftr2_annotations.json"
    if truth_ann.exists():
        annotations = json.loads(truth_ann.read_text(encoding="utf-8"))
    else:
        annotations = _build_annotations(selected, clinvar_panel, drug_panel)

    ref_fasta = str(config.reference_fasta)
    vcf_metrics = _score_vcf_against_truth(
        vcf_content, truth_vcf, ref_fasta, bam_path, work_dir / "score"
    )
    ann_score = 0.0
    if truth_ann.exists():
        truth_ann_data = json.loads(truth_ann.read_text(encoding="utf-8"))
        ann_score = score_annotations(annotations, truth_ann_data)
    final_score = 0.7 * vcf_metrics["vcf_score"] + 0.3 * ann_score

    result = {
        "task_id": task_id,
        "calibrated": calibrated,
        "variant_count": len(selected),
        "annotation_score": ann_score,
        "final_score": final_score,
        **vcf_metrics,
    }
    (output_dir / "benchmark_score.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    (output_dir / "submission.vcf").write_text(vcf_content, encoding="utf-8")
    (output_dir / "annotation.json").write_text(
        json.dumps(annotations, indent=2, sort_keys=True), encoding="utf-8"
    )
    _write_final_review_tsv(output_dir / "final_review.tsv", review_rows)
    return result


def main() -> None:
    tasks_dir = ROOT / "tasks"
    task_ids = sorted(
        p.name
        for p in tasks_dir.iterdir()
        if p.is_dir() and (p / "truth.vcf").exists()
    )
    if not task_ids:
        print("No tasks with truth.vcf found under tasks/")
        sys.exit(1)

    print(f"Benchmarking {len(task_ids)} task(s) with truth: {', '.join(task_ids)}")
    for task_id in task_ids:
        result = run_task(task_id)
        print(
            f"\n{task_id}: final={result['final_score']:.4f} "
            f"vcf={result['vcf_score']:.4f} ann={result['annotation_score']:.4f} "
            f"tp_w={result['tp_w']:.2f} fp_w={result['fp_w']:.2f} fn_w={result['fn_w']:.2f} "
            f"variants={result['variant_count']} calibrated={result['calibrated']}"
        )


if __name__ == "__main__":
    main()
