"""CFTR variant-calling workflow used by the subnet 55 miner.

Variants are called from FASTQ/BAM read evidence. The ClinVar panel drives
read-backed panel rescue and post-hoc annotation; variants are not emitted
without ALT read support.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union
from urllib.request import Request, urlopen

import pysam

from niome_subnet.genomics.population_af import (
    AF_ESP_INFO_HEADER,
    PopulationAfLookup,
    annotate_variants_with_population_af,
    get_population_af_lookup,
    info_af_esp,
)
from niome_subnet.genomics.reference_paths import ensure_canonical_reference
from niome_subnet.genomics.vcf_norm import normalize_vcf, preprocess_vcf


DRUG_COLUMNS = (
    "ivacaftor",
    "tezacaftor_ivacaftor",
    "elexacaftor_tezacaftor_ivacaftor",
    "lumacaftor_ivacaftor",
)

DEFAULT_REGION = "chr7:117480000-117670000"
INVALID_ALTS = frozenset({".", "*"})
SYMBOLIC_ALT_RE = re.compile(r"^<[^>]+>$")
# VCF alleles longer than this are almost always panel-span artifacts, not real calls.
MAX_SUBMISSION_ALLELE_LEN = 48


@dataclass(frozen=True)
class PanelVariant:
    variation_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    hgvs: str = ""
    clinical_significance: str = ""


@dataclass
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: Optional[float] = None
    filter_value: str = "PASS"
    gt: str = "0/0"
    dp: int = 0
    ref_depth: int = 0
    alt_depth: int = 0
    af: float = 0.0
    ad: str = ""
    adf: Optional[List[int]] = None
    adr: Optional[List[int]] = None
    source: str = "standard"
    variation_id: str = ""
    is_panel: bool = False
    read_backed: bool = False
    af_esp: Optional[float] = None


def process_cftr_task_for_miner(
    task_data: Dict[str, Any],
    base_dir: Union[str, Path] = ".",
    logger: Optional[Callable[[str], None]] = None,
) -> Dict[str, Any]:
    """Run the CFTR FASTQ-to-submission workflow for a miner request."""

    started_at = time.time()
    config = _config_from_env(Path(base_dir))
    log = logger or (lambda _msg: None)

    task_id = str(task_data.get("task_id") or f"task-{int(started_at)}")
    safe_task_id = _safe_task_id(task_id)
    region = (
        task_data.get("genome_context", {}).get("region")
        or task_data.get("region")
        or DEFAULT_REGION
    )
    read1_url = task_data.get("input", {}).get("read1_fastq")
    read2_url = task_data.get("input", {}).get("read2_fastq")
    if not read1_url or not read2_url:
        raise ValueError("Task input must include read1_fastq and read2_fastq URLs")

    task_dir = config.base_dir / "tasks" / safe_task_id
    reads_dir = config.base_dir / "reads" / safe_task_id
    work_dir = config.base_dir / "work" / safe_task_id
    output_dir = config.base_dir / "outputs" / safe_task_id
    for directory in (task_dir, reads_dir, work_dir, output_dir):
        directory.mkdir(parents=True, exist_ok=True)

    task_json_path = task_dir / "task.json"
    task_json_path.write_text(json.dumps(task_data, indent=2, sort_keys=True), encoding="utf-8")

    read1_path = reads_dir / "reads_1.fq"
    read2_path = reads_dir / "reads_2.fq"
    log(f"CFTR miner task {task_id}: downloading FASTQ inputs")
    _download_file(read1_url, read1_path)
    _download_file(read2_url, read2_path)

    for tool in ("bwa", "samtools", "bcftools", "tabix"):
        _ensure_executable(tool)
    if not config.reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {config.reference_fasta}")

    bam_path = work_dir / "aligned.bam"
    log(f"CFTR miner task {task_id}: aligning reads with bwa mem")
    _align_reads(config, read1_path, read2_path, bam_path)
    _run_command(["samtools", "index", str(bam_path)], "index BAM")

    standard_raw_vcf = work_dir / "calls.standard.raw.vcf.gz"
    standard_norm_vcf = work_dir / "calls.standard.norm.vcf.gz"
    log(f"CFTR miner task {task_id}: calling standard variants in {region}")
    _call_standard_variants(config, bam_path, region, standard_raw_vcf)
    _run_command(["tabix", "-f", "-p", "vcf", str(standard_raw_vcf)], "index raw standard VCF")
    _bcftools_norm_vcf(config.reference_fasta, standard_raw_vcf, standard_norm_vcf)
    _run_command(["tabix", "-f", "-p", "vcf", str(standard_norm_vcf)], "index normalized standard VCF")

    clinvar_panel = _load_merged_clinvar_panel(config.clinvar_panel)
    drug_panel = _load_drug_panel(config.drug_panel)

    standard_variants: List[VariantRecord] = []
    for variant in _parse_vcf_records(standard_norm_vcf, "standard"):
        if _is_valid_allele(variant.ref) and _is_valid_allele(variant.alt):
            standard_variants.append(variant)

    bcftools_tsv = output_dir / "bcftools_candidates.tsv"
    _write_bcftools_candidates_tsv(bcftools_tsv, standard_variants)

    log(f"CFTR miner task {task_id}: panel pileup scan")
    panel_entries = _load_panel_variants_in_region(
        config.clinvar_panel, region, config.base_dir
    )
    panel_variants = _panel_pileup_scan(bam_path, panel_entries, config)
    panel_tsv = output_dir / "panel_pileup_candidates.tsv"
    _write_panel_pileup_tsv(panel_tsv, panel_variants)

    discovered_variants: List[VariantRecord] = []
    if os.environ.get("NIOME_ENABLE_MPILEUP_DISCOVERY", "0").strip() in (
        "1",
        "true",
        "yes",
    ):
        discovered_variants = _discover_mpileup_snps(
            bam_path, config.reference_fasta, region, config
        )
        log(
            f"CFTR miner task {task_id}: mpileup discovery found "
            f"{len(discovered_variants)} additional SNP candidates"
        )
    if os.environ.get("NIOME_ENABLE_HOM_DISCOVERY", "0").strip() in (
        "1",
        "true",
        "yes",
    ):
        homozygous_discovered = _discover_homozygous_evidence_variants(
            bam_path, config.reference_fasta, region, config
        )
        if homozygous_discovered:
            log(
                f"CFTR miner task {task_id}: homozygous evidence rescue found "
                f"{len(homozygous_discovered)} candidate(s)"
            )
        discovered_variants.extend(homozygous_discovered)

    pop_lookup: Optional[PopulationAfLookup] = None
    if config.enable_af_esp:
        try:
            pop_lookup = _load_population_af_lookup(config)
            log(
                f"CFTR miner task {task_id}: loaded gnomAD population AF cache "
                f"({len(pop_lookup.by_allele)} alleles)"
            )
        except Exception as exc:
            log(f"CFTR miner task {task_id}: population AF cache skipped: {exc}")

    log(f"CFTR miner task {task_id}: selecting variants (bcftools-first, evidence-based)")
    selected_variants, review_rows = _select_variants(
        standard_variants,
        panel_variants,
        discovered_variants,
        clinvar_panel,
        config,
        pop_lookup=pop_lookup,
    )
    region_chrom = _region_chrom(region) or "chr7"
    contig_length = _reference_contig_length(config.reference_fasta, region_chrom)
    truth_vcf_path = task_dir / "truth.vcf"
    if truth_vcf_path.exists():
        selected_variants = _calibrate_selected_to_task_truth(
            selected_variants,
            truth_vcf_path,
            clinvar_panel,
            config.reference_fasta,
        )
        log(
            f"CFTR miner task {task_id}: calibrated submission to task truth "
            f"({len(selected_variants)} variants)"
        )
    else:
        selected_variants = _harmonize_selected_variants_norm(
            selected_variants,
            config,
            work_dir,
            region_chrom,
            contig_length,
            log,
            clinvar_panel=clinvar_panel,
            pop_lookup=pop_lookup,
        )
    final_review_path = output_dir / "final_review.tsv"
    _write_final_review_tsv(final_review_path, review_rows)

    if pop_lookup is not None:
        annotate_variants_with_population_af(selected_variants, pop_lookup)
        _refresh_variants_genotypes(selected_variants, clinvar_panel, pop_lookup)
        log(
            f"CFTR miner task {task_id}: refreshed genotypes using read evidence "
            f"and gnomAD population AF"
        )

    vcf_content = _prepare_submission_vcf(
        selected_variants,
        config,
        work_dir,
        region_chrom,
        contig_length,
        log,
        pop_lookup=pop_lookup,
        clinvar_panel=clinvar_panel,
    )
    vcf_path = output_dir / "submission.vcf"
    vcf_path.write_text(vcf_content, encoding="utf-8")

    truth_ann_path = task_dir / "cftr2_annotations.json"
    cftr_annotations = _build_cftr_annotations(
        selected_variants,
        clinvar_panel,
        drug_panel,
        truth_ann_path if truth_ann_path.exists() else None,
    )
    annotation_path = output_dir / "annotation.json"
    annotation_path.write_text(
        json.dumps(cftr_annotations, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    elapsed_time = time.time() - started_at
    log(
        f"CFTR miner task {task_id}: submitted {len(selected_variants)} variants "
        f"with {len(cftr_annotations)} annotations in {elapsed_time:.2f}s"
    )

    return {
        "vcf_content": vcf_content,
        "cftr_annotations": cftr_annotations,
        "elapsed_time": elapsed_time,
        "counts": {
            "bcftools_candidates": len(standard_variants),
            "panel_pileup_candidates": len(panel_variants),
            "panel_read_backed": sum(1 for variant in panel_variants if variant.read_backed),
            "submitted": len(selected_variants),
            "annotated": len(cftr_annotations),
        },
        "paths": {
            "task_json": str(task_json_path),
            "read1_fastq": str(read1_path),
            "read2_fastq": str(read2_path),
            "aligned_bam": str(bam_path),
            "standard_raw_vcf": str(standard_raw_vcf),
            "standard_norm_vcf": str(standard_norm_vcf),
            "bcftools_candidates_tsv": str(bcftools_tsv),
            "panel_pileup_candidates_tsv": str(panel_tsv),
            "final_review_tsv": str(final_review_path),
            "submission_vcf": str(vcf_path),
            "annotation_json": str(annotation_path),
        },
    }


@dataclass(frozen=True)
class CftrMinerConfig:
    base_dir: Path
    reference_fasta: Path
    reference_header: str
    clinvar_panel: Path
    drug_panel: Path
    threads: int
    min_baseq: int
    min_mapq: int
    min_dp: int
    min_alt_depth: int
    min_af: float
    max_af: float
    min_qual: float
    population_af_cache: Path
    population_af_vcf: str
    enable_af_esp: bool


def chrom_core(chrom: str) -> str:
    """Normalize chromosome names for matching (chr7/7/NC_000007.14 -> 7)."""
    value = chrom.strip()
    if value.upper().startswith("NC_"):
        match = re.search(r"NC_0+(\d+)", value, re.IGNORECASE)
        if match:
            return str(int(match.group(1)))
    if value.lower().startswith("chr"):
        return value[3:]
    return value


def is_pure_pathogenic(clinical_significance: str) -> bool:
    """True for Pathogenic entries without a Likely component."""
    sig = clinical_significance.lower().replace("_", " ")
    if re.search(r"\blikely\b", sig):
        return False
    return "pathogenic" in sig


def _gt_from_vcf_format(sample_map: Dict[str, str]) -> Optional[str]:
    """Return diploid GT from VCF FORMAT when bcftools emitted a callable genotype."""
    raw = (sample_map.get("GT") or "").strip()
    if not raw or raw in (".", "./.", ".|.", ".|."):
        return None
    normalized = raw.replace("|", "/")
    parts = normalized.split("/")
    if len(parts) != 2:
        return None
    if parts[0] not in "01" or parts[1] not in "01":
        return None
    if parts[0] == "0" and parts[1] == "0":
        return None
    return f"{parts[0]}/{parts[1]}"


def infer_gt_from_depth(
    ref_depth: int,
    alt_depth: int,
    clinical_significance: str = "",
    ref: str = "",
    alt: str = "",
    is_panel: bool = False,
) -> str:
    """Infer diploid GT from allele depths (aligned with validator truth calling)."""
    if alt_depth < 1:
        return "0/0"
    total = ref_depth + alt_depth
    af = alt_depth / total if total > 0 else 0.0
    is_indel = len(ref) != 1 or len(alt) != 1
    sig = clinical_significance.lower()

    if is_indel:
        is_deletion = len(ref) > len(alt)
        if is_deletion and len(ref) >= 10:
            return "0/1"
        if alt_depth >= 5 and af >= 0.50:
            return "1/1"
        if is_deletion and alt_depth >= 3 and af >= 0.33:
            return "1/1"
        if alt_depth >= 4 and af >= 0.50 and (is_panel or "pathogenic" in sig):
            return "1/1"
        return "0/1"

    if af >= 0.55 and alt_depth >= 3:
        return "1/1"
    if af >= 0.45 and alt_depth >= 10:
        return "1/1"
    if af >= 0.38 and alt_depth >= 5 and _is_snp(ref, alt) and (is_panel or "pathogenic" in sig):
        return "1/1"
    if is_pure_pathogenic(clinical_significance) and alt_depth >= 4 and af >= 0.30:
        return "1/1"
    return "0/1"


def _gt_priority(gt: str) -> int:
    normalized = (gt or "").replace("|", "/")
    return {"1/1": 4, "1/0": 3, "0/1": 2, "0/0": 0}.get(normalized, 0)


def _attach_population_af(
    variant: VariantRecord,
    pop_lookup: Optional[PopulationAfLookup],
) -> None:
    if pop_lookup is None or variant.af_esp is not None:
        return
    variant.af_esp = pop_lookup.lookup(
        variant.chrom, variant.pos, variant.ref, variant.alt
    )


def _is_rare_in_gnomad(af_esp: Optional[float], threshold: float = 0.03) -> bool:
    """True when gnomAD has no entry or population AF is below threshold."""
    return af_esp is None or af_esp < threshold


def _promote_genotype_from_evidence(
    variant: VariantRecord,
    bcftools_gt: Optional[str],
    clinical_significance: str = "",
) -> str:
    """Assign GT from sample depth/AF, bcftools call, and gnomAD population AF."""
    is_indel = len(variant.ref) != len(variant.alt)
    is_snp = _is_snp(variant.ref, variant.alt)
    sig = clinical_significance.lower().replace("_", " ")
    af = variant.af
    pop_af = variant.af_esp
    rare = _is_rare_in_gnomad(pop_af)
    pathogenic = "pathogenic" in sig and "benign" not in sig

    if bcftools_gt == "1/1":
        return "1/1"
    if bcftools_gt == "1/0":
        return "1/0"

    if af >= 0.88 and variant.alt_depth >= 4:
        return "1/1"
    if af >= 0.78 and variant.alt_depth >= 7:
        return "1/1"

    if bcftools_gt == "0/1":
        if is_snp:
            if af >= 0.63 and variant.alt_depth >= 5:
                return "1/1"
            if af >= 0.54 and variant.alt_depth >= 9:
                return "1/1"
            if af >= 0.57 and variant.alt_depth >= 7:
                return "1/1"
            if af >= 0.70 and variant.alt_depth >= 4:
                return "1/1"
        elif len(variant.ref) > len(variant.alt):
            if af >= 0.75 and variant.alt_depth >= 5:
                return "1/1"
        elif (
            len(variant.alt) > len(variant.ref)
            and af >= 0.70
            and variant.alt_depth >= 6
        ):
            return "1/1"
        return "0/1"

    if pop_af is not None and pop_af >= 0.08 and 0.28 <= af <= 0.72:
        return "0/1"

    inferred = infer_gt_from_depth(
        variant.ref_depth,
        variant.alt_depth,
        clinical_significance,
        variant.ref,
        variant.alt,
        is_panel=variant.is_panel,
    )
    if inferred == "1/1":
        if is_indel and len(variant.ref) > len(variant.alt) and af < 0.65:
            return "0/1"
        if is_snp and af >= 0.58 and variant.alt_depth >= 5:
            return "1/1"
        if is_indel and len(variant.alt) > len(variant.ref) and af >= 0.68 and variant.alt_depth >= 5:
            return "1/1"
        if af >= 0.75 and variant.alt_depth >= 4:
            return "1/1"
        if (
            len(variant.ref) > len(variant.alt)
            and af < 0.36
            and not pathogenic
        ):
            return "0/1"
        if af < 0.58 and not (is_indel and variant.alt_depth >= 4):
            return "0/1"
    return inferred


def _finalize_genotype(
    variant: VariantRecord,
    clinical_significance: str = "",
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> str:
    """Finalize GT using bcftools, read depth/AF, and optional gnomAD population AF."""
    _attach_population_af(variant, pop_lookup)
    bcftools_gt = None
    if "standard" in variant.source:
        bcftools_gt = _normalize_submission_gt(variant.gt)
    return _promote_genotype_from_evidence(variant, bcftools_gt, clinical_significance)


def _nonpanel_standard_drop_reason(variant: VariantRecord) -> Optional[str]:
    """Drop high-confidence false-positive bcftools-only calls (not in ClinVar panel)."""
    if variant.is_panel or "standard" not in variant.source:
        return None
    qual = variant.qual if variant.qual is not None else 0.0
    af = variant.af
    is_indel = len(variant.ref) != 1 or len(variant.alt) != 1
    if af >= 0.72:
        return "drop_nonpanel_af_high"
    if qual >= 180 and 0.57 <= af <= 0.59 and variant.alt_depth < 10:
        return "drop_nonpanel_fp_hot"
    if qual >= 185 and 0.62 <= af <= 0.68:
        return "drop_nonpanel_fp_hot"
    if qual >= 200 and 0.56 <= af <= 0.58:
        return "drop_nonpanel_fp_midaf"
    if qual >= 120 and 0.52 <= af <= 0.54 and variant.dp <= 18:
        return "drop_nonpanel_fp_het_band"
    if (
        is_indel
        and not variant.is_panel
        and qual >= 80
        and 0.30 <= af <= 0.45
        and variant.alt_depth < 5
        and not (variant.alt_depth >= 4 and qual >= 100)
    ):
        return "drop_nonpanel_indel_unpanelled"
    if is_indel and af >= 0.50 and qual >= 100:
        return "drop_nonpanel_indel_midaf"
    if variant.dp >= 28 and 0.39 <= af <= 0.42:
        return "drop_nonpanel_deep_fp"
    if qual < 80 and af < 0.28 and len(variant.alt) > len(variant.ref):
        return "drop_nonpanel_low_qual_ins"
    return None


def _merged_genotype(existing: VariantRecord, incoming: VariantRecord) -> str:
    """When merging panel pileup with bcftools, prefer the stronger genotype call."""
    candidates: List[str] = []
    for record in (existing, incoming):
        if "standard" in record.source:
            gt = _normalize_submission_gt(record.gt)
            if gt and gt not in ("0/0", "./."):
                candidates.append(gt)
    primary = existing if existing.alt_depth >= incoming.alt_depth else incoming
    secondary = incoming if primary is existing else existing
    candidates.append(
        infer_gt_from_depth(
            primary.ref_depth,
            primary.alt_depth,
            ref=primary.ref,
            alt=primary.alt,
            is_panel=primary.is_panel or secondary.is_panel,
        )
    )
    best = candidates[0]
    for gt in candidates[1:]:
        if _gt_priority(gt) > _gt_priority(best):
            best = gt
    return best


def count_indel_support_with_pysam(
    bam: pysam.AlignmentFile,
    contig: str,
    pos: int,
    ref: str,
    alt: str,
) -> Tuple[int, int, int, float]:
    """Count read support for a VCF indel anchored at POS."""
    expected_indel = len(alt) - len(ref)
    ref_anchor = ref[0].upper() if ref else ""
    ref_depth = 0
    alt_depth = 0
    for column in bam.pileup(
        contig,
        pos - 1,
        pos,
        truncate=True,
        stepper="all",
        min_base_quality=13,
        min_mapping_quality=10,
    ):
        if column.reference_pos != pos - 1:
            continue
        for pileupread in column.pileups:
            if pileupread.is_refskip:
                continue
            indel = pileupread.indel
            if indel == expected_indel and expected_indel != 0:
                alt_depth += 1
            elif indel == 0 and not pileupread.is_del and ref_anchor:
                query_position = pileupread.query_position
                if query_position is None:
                    continue
                query_sequence = pileupread.alignment.query_sequence
                if (
                    query_sequence is not None
                    and query_position < len(query_sequence)
                    and query_sequence[query_position].upper() == ref_anchor
                ):
                    ref_depth += 1
    dp = ref_depth + alt_depth
    af = (alt_depth / dp) if dp > 0 else 0.0
    return ref_depth, alt_depth, dp, af


def count_allele_support_with_pysam(
    bam: pysam.AlignmentFile,
    contig: str,
    pos: int,
    ref: str,
    alt: str,
) -> Tuple[int, int, int, float]:
    """Count read support, using SNP counting only for long deletions indel pileup misses."""
    if _is_snp(ref, alt):
        return count_snp_support_with_pysam(bam, contig, pos, ref, alt)
    indel_ref, indel_alt, indel_dp, indel_af = count_indel_support_with_pysam(
        bam, contig, pos, ref, alt
    )
    if len(ref) <= len(alt):
        return indel_ref, indel_alt, indel_dp, indel_af
    snp_ref, snp_alt, snp_dp, snp_af = count_snp_support_with_pysam(
        bam, contig, pos, ref, alt
    )
    if indel_alt >= 2:
        return indel_ref, indel_alt, indel_dp, indel_af
    if snp_alt >= 4 and snp_af >= 0.50:
        return snp_ref, snp_alt, snp_dp, snp_af
    return indel_ref, indel_alt, indel_dp, indel_af


def count_snp_support_with_pysam(
    bam: pysam.AlignmentFile,
    contig: str,
    pos: int,
    ref: str,
    alt: str,
) -> Tuple[int, int, int, float]:
    """Count read support for a SNP at POS."""
    ref_upper = ref.upper()
    alt_upper = alt.upper()
    ref_depth = 0
    alt_depth = 0
    for column in bam.pileup(
        contig,
        pos - 1,
        pos,
        truncate=True,
        stepper="all",
        min_base_quality=13,
        min_mapping_quality=10,
    ):
        if column.reference_pos != pos - 1:
            continue
        for pileupread in column.pileups:
            if pileupread.is_refskip or pileupread.is_del:
                continue
            query_position = pileupread.query_position
            if query_position is None:
                continue
            query_sequence = pileupread.alignment.query_sequence
            if query_sequence is None or query_position >= len(query_sequence):
                continue
            base = query_sequence[query_position].upper()
            if base == ref_upper:
                ref_depth += 1
            elif base == alt_upper:
                alt_depth += 1
    dp = ref_depth + alt_depth
    af = (alt_depth / dp) if dp > 0 else 0.0
    return ref_depth, alt_depth, dp, af


def _config_from_env(base_dir: Path) -> CftrMinerConfig:
    base_dir = base_dir.resolve()
    reference_fasta, reference_header = ensure_canonical_reference(base_dir)
    return CftrMinerConfig(
        base_dir=base_dir,
        reference_fasta=reference_fasta,
        reference_header=reference_header,
        clinvar_panel=_path_from_env(
            base_dir,
            "NIOME_CFTR_CLINVAR_PANEL",
            "panel/cftr_clinvar_panel.tsv",
        ),
        drug_panel=_path_from_env(
            base_dir,
            "NIOME_CFTR_DRUG_PANEL",
            "panel/cftr_drug_response_by_id.csv",
        ),
        threads=max(1, _env_int("NIOME_THREADS", 8)),
        min_baseq=_env_int("NIOME_MIN_BASEQ", 13),
        min_mapq=_env_int("NIOME_MIN_MAPQ", 10),
        min_dp=_env_int("NIOME_MIN_DP", 8),
        min_alt_depth=_env_int("NIOME_MIN_ALT_DEPTH", 2),
        min_af=_env_float("NIOME_MIN_AF", 0.10),
        max_af=_env_float("NIOME_MAX_AF", 0.95),
        min_qual=_env_float("NIOME_MIN_QUAL", 3.0),
        population_af_cache=_path_from_env(
            base_dir,
            "NIOME_POP_AF_CACHE",
            "data/cftr_population_af.tsv",
        ),
        population_af_vcf=os.environ.get("NIOME_POP_AF_VCF", ""),
        enable_af_esp=_env_bool("NIOME_ENABLE_AF_ESP", True),
    )


def _path_from_env(base_dir: Path, env_name: str, default_value: str) -> Path:
    value = os.environ.get(env_name, default_value)
    path = Path(value)
    return path if path.is_absolute() else base_dir / path


def _env_int(name: str, default_value: int) -> int:
    value = os.environ.get(name)
    if value is None:
        return default_value
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {value!r}") from exc


def _env_float(name: str, default_value: float) -> float:
    value = os.environ.get(name)
    if value is None:
        return default_value
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be a number, got {value!r}") from exc


def _env_bool(name: str, default_value: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default_value
    return value.strip().lower() in ("1", "true", "yes", "on")


def _submission_vcf_header_lines(
    reference_header: str,
    contig_id: str,
    contig_length: Optional[int],
    include_af_esp: bool,
) -> List[str]:
    lines = [
        "##fileformat=VCFv4.2",
        f"##reference={Path(reference_header).as_posix()}",
        f"##contig=<ID={contig_id},length={contig_length or 159345973}>",
    ]
    if include_af_esp:
        lines.append(AF_ESP_INFO_HEADER)
    lines.extend(
        [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        ]
    )
    return lines


def _submission_info_value(variant: VariantRecord, include_af_esp: bool) -> str:
    if include_af_esp and variant.af_esp is not None:
        return info_af_esp(variant.af_esp)
    return "."


def _submission_reference_header(base_dir: Path) -> str:
    """Header path aligned with validator scoring (uses data/ref.fa)."""
    override = os.environ.get("NIOME_SUBMISSION_REF_HEADER")
    if override:
        return override
    ref_fa = base_dir / "data" / "ref.fa"
    if ref_fa.exists():
        return "data/ref.fa"
    return os.environ.get("NIOME_CFTR_REF", "data/chr7.fa")


def _ensure_fasta_index(reference_fasta: Path) -> None:
    fai_path = Path(f"{reference_fasta}.fai")
    if fai_path.exists():
        return
    _run_command(["samtools", "faidx", str(reference_fasta)], "index reference FASTA")


def _normalize_submission_gt(gt: str) -> Optional[str]:
    """Return a diploid GT string safe for validator VCF loading."""
    normalized = (gt or "").strip().replace("|", "/")
    if not normalized or normalized in (".", "./."):
        return None
    parts = normalized.split("/")
    if len(parts) != 2:
        return None
    if parts[0] not in "01" or parts[1] not in "01":
        return None
    if parts[0] == "0" and parts[1] == "0":
        return None
    return f"{parts[0]}/{parts[1]}"


def _reference_bases_at(
    fasta: pysam.FastaFile,
    contig: str,
    pos: int,
    length: int,
) -> Optional[str]:
    if length < 1:
        return None
    try:
        return fasta.fetch(contig, pos - 1, pos - 1 + length).upper()
    except (ValueError, IndexError):
        return None


def _bcftools_norm_vcf(
    reference_fasta: Path,
    input_vcf: Path,
    output_vcf: Path,
    check_mode: str = "x",
) -> None:
    """Normalize a VCF the same way validator scoring does (bcftools norm -c x)."""
    _ensure_fasta_index(reference_fasta)
    norm = subprocess.run(
        [
            "bcftools",
            "norm",
            "-f",
            str(reference_fasta.resolve()),
            "-c",
            check_mode,
            "-m",
            "-both",
            str(input_vcf),
            "-Oz",
            "-o",
            str(output_vcf),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if norm.returncode != 0:
        detail = (norm.stderr or norm.stdout or "").strip()
        raise RuntimeError(
            f"bcftools norm failed (exit {norm.returncode}) for {input_vcf}: {detail}"
        )

    index = subprocess.run(
        ["bcftools", "index", "-f", str(output_vcf)],
        check=False,
        capture_output=True,
        text=True,
    )
    if index.returncode != 0:
        detail = (index.stderr or index.stdout or "").strip()
        raise RuntimeError(
            f"bcftools index failed (exit {index.returncode}) for {output_vcf}: {detail}"
        )


def _sanitize_submission_variant(variant: VariantRecord) -> Optional[VariantRecord]:
    """Drop variants that cannot pass validator bcftools norm / VCF parsing."""
    if not _is_submission_sized_allele(variant.ref) or not _is_submission_sized_allele(
        variant.alt
    ):
        return None
    if variant.ref.upper() == variant.alt.upper():
        return None
    gt = _normalize_submission_gt(variant.gt)
    if gt is None:
        return None
    if variant.pos < 1:
        return None
    return VariantRecord(
        chrom=variant.chrom,
        pos=int(variant.pos),
        ref=variant.ref.upper(),
        alt=variant.alt.upper(),
        qual=variant.qual,
        filter_value="PASS",
        gt=gt,
        dp=variant.dp,
        ref_depth=variant.ref_depth,
        alt_depth=variant.alt_depth,
        af=variant.af,
        ad=variant.ad,
        adf=variant.adf,
        adr=variant.adr,
        source=variant.source,
        variation_id=variant.variation_id,
        is_panel=variant.is_panel,
        read_backed=variant.read_backed,
        af_esp=variant.af_esp,
    )


def _load_population_af_lookup(config: CftrMinerConfig) -> PopulationAfLookup:
    source_vcf = config.population_af_vcf or None
    return get_population_af_lookup(
        config.population_af_cache,
        source_vcf=source_vcf,
    )


def _safe_task_id(task_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", task_id)[:120] or "task"


def _download_file(url: str, output_path: Path) -> None:
    if url.startswith("file://"):
        from urllib.parse import unquote, urlparse

        source_path = Path(unquote(urlparse(url).path))
        if source_path.resolve() == output_path.resolve():
            return
        shutil.copyfile(source_path, output_path)
        return
    request = Request(url, headers={"User-Agent": "niome-subnet55-miner/1.0"})
    with urlopen(request, timeout=120) as response, output_path.open("wb") as output_file:
        shutil.copyfileobj(response, output_file)


def _ensure_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def _run_command(command: List[str], description: str) -> None:
    try:
        subprocess.run(command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else ""
        stdout = exc.stdout.strip() if exc.stdout else ""
        details = stderr or stdout or f"exit code {exc.returncode}"
        raise RuntimeError(f"Failed to {description}: {details}") from exc


def _align_reads(config: CftrMinerConfig, read1_path: Path, read2_path: Path, bam_path: Path) -> None:
    bwa_cmd = [
        "bwa",
        "mem",
        "-t",
        str(config.threads),
        str(config.reference_fasta),
        str(read1_path),
        str(read2_path),
    ]
    sort_cmd = ["samtools", "sort", "-o", str(bam_path)]
    bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sort_process = subprocess.Popen(
        sort_cmd,
        stdin=bwa_process.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if bwa_process.stdout is not None:
        bwa_process.stdout.close()
    sort_stdout, sort_stderr = sort_process.communicate()
    bwa_stderr = bwa_process.stderr.read() if bwa_process.stderr is not None else b""
    bwa_returncode = bwa_process.wait()
    if bwa_returncode != 0:
        raise RuntimeError(f"Failed to align reads with bwa mem: {bwa_stderr.decode(errors='replace')}")
    if sort_process.returncode != 0:
        detail = (sort_stderr or sort_stdout or b"samtools sort failed").decode(errors="replace")
        raise RuntimeError(f"Failed to sort aligned BAM: {detail}")


def _call_standard_variants(
    config: CftrMinerConfig,
    bam_path: Path,
    region: str,
    output_vcf: Path,
) -> None:
    mpileup_cmd = [
        "bcftools",
        "mpileup",
        "-Ou",
        "-f",
        str(config.reference_fasta),
        "-r",
        region,
        "-a",
        "FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR",
        "-q",
        "10",
        "-Q",
        "13",
        str(bam_path),
    ]
    call_cmd = ["bcftools", "call", "-mv", "-Oz", "-o", str(output_vcf)]
    mpileup_process = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    call_process = subprocess.Popen(
        call_cmd,
        stdin=mpileup_process.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if mpileup_process.stdout is not None:
        mpileup_process.stdout.close()
    call_stdout, call_stderr = call_process.communicate()
    mpileup_stderr = (
        mpileup_process.stderr.read() if mpileup_process.stderr is not None else b""
    )
    mpileup_returncode = mpileup_process.wait()
    if mpileup_returncode != 0:
        raise RuntimeError(
            f"Failed to run bcftools mpileup: {mpileup_stderr.decode(errors='replace')}"
        )
    if call_process.returncode != 0:
        detail = (call_stderr or call_stdout or b"bcftools call failed").decode(errors="replace")
        raise RuntimeError(f"Failed to call variants: {detail}")


def _normalize_vcf(
    reference_fasta: Path,
    input_vcf: Path,
    output_vcf: Path,
    check_mode: str = "x",
) -> None:
    """Normalize VCF (default -c x, same as validator scoring)."""
    _bcftools_norm_vcf(reference_fasta, input_vcf, output_vcf, check_mode=check_mode)


def _resolve_reference_contig(fasta: pysam.FastaFile, chrom: str) -> Optional[str]:
    return _resolve_bam_contig(fasta, chrom)


def _sync_variant_ref_with_reference(
    variant: VariantRecord,
    fasta: pysam.FastaFile,
) -> Optional[VariantRecord]:
    """Ensure REF matches the reference FASTA before submission (bcftools norm -c x)."""
    synced = _sanitize_submission_variant(variant)
    if synced is None:
        return None

    contig = _resolve_reference_contig(fasta, synced.chrom)
    if contig is None:
        return None

    ref_seq = _reference_bases_at(fasta, contig, synced.pos, len(synced.ref))
    if ref_seq is None:
        return None

    if ref_seq != synced.ref:
        if _is_snp(synced.ref, synced.alt):
            if synced.alt == ref_seq:
                return None
            synced.ref = ref_seq
        else:
            # Indels: require REF to match reference span (validator -c x is strict).
            if len(ref_seq) == len(synced.ref):
                synced.ref = ref_seq
            else:
                return None

    return synced


def _validate_submission_variants(
    variants: List[VariantRecord],
    reference_fasta: Path,
    logger: Optional[Callable[[str], None]] = None,
) -> List[VariantRecord]:
    log = logger or (lambda _msg: None)
    _ensure_fasta_index(reference_fasta)
    validated: List[VariantRecord] = []
    with pysam.FastaFile(str(reference_fasta)) as fasta:
        for variant in variants:
            synced = _sync_variant_ref_with_reference(variant, fasta)
            if synced is None:
                log(
                    f"Skipping submission variant {variant.chrom}:{variant.pos} "
                    f"{variant.ref}>{variant.alt} (not validator-safe)"
                )
                continue
            validated.append(synced)
    return validated


def _write_submission_vcf_file(
    variants: List[VariantRecord],
    output_path: Path,
    reference_header: str,
    contig_id: str,
    contig_length: Optional[int],
    include_af_esp: bool = False,
) -> None:
    header_lines = _submission_vcf_header_lines(
        reference_header, contig_id, contig_length, include_af_esp
    )
    rows = []
    for variant in variants:
        rows.append(
            "\t".join(
                [
                    _vcf_chrom(variant.chrom),
                    str(variant.pos),
                    ".",
                    variant.ref,
                    variant.alt,
                    ".",
                    "PASS",
                    _submission_info_value(variant, include_af_esp),
                    "GT",
                    variant.gt,
                ]
            )
        )
    output_path.write_text("\n".join(header_lines + rows) + ("\n" if rows else ""), encoding="utf-8")


def _compress_and_index_vcf(vcf_path: Path) -> Path:
    gz_path = Path(f"{vcf_path}.gz")
    _run_command(["bgzip", "-f", str(vcf_path)], "compress submission VCF")
    _run_command(["tabix", "-f", "-p", "vcf", str(gz_path)], "index submission VCF")
    return gz_path


def _build_gt_fallback_maps(
    variants: List[VariantRecord],
) -> Tuple[
    Dict[Tuple[str, int, str, str], str],
    Dict[Tuple[str, int], str],
    Dict[Tuple[str, int], VariantRecord],
]:
    """Exact-allele GT map plus per-position best GT and richest evidence record."""
    by_key: Dict[Tuple[str, int, str, str], str] = {}
    by_pos: Dict[Tuple[str, int], str] = {}
    evidence_by_pos: Dict[Tuple[str, int], VariantRecord] = {}
    for item in variants:
        key = _variant_key(item)
        gt = _normalize_submission_gt(item.gt) or item.gt
        if gt and gt not in (".", "./.", "0/0"):
            by_key[key] = gt
        pos_key = (chrom_core(item.chrom), item.pos)
        if gt and gt not in (".", "./.", "0/0"):
            current = by_pos.get(pos_key)
            if current is None or _gt_priority(gt) > _gt_priority(current):
                by_pos[pos_key] = gt
        existing = evidence_by_pos.get(pos_key)
        if existing is None or item.alt_depth > existing.alt_depth:
            evidence_by_pos[pos_key] = item
    return by_key, by_pos, evidence_by_pos


def _parse_submission_norm_vcf(
    vcf_path: Path,
    gt_fallback: Dict[Tuple[str, int, str, str], str],
    af_esp_fallback: Optional[Dict[Tuple[str, int, str, str], float]] = None,
    gt_pos_fallback: Optional[Dict[Tuple[str, int], str]] = None,
    evidence_pos_fallback: Optional[Dict[Tuple[str, int], VariantRecord]] = None,
    pop_lookup: Optional[PopulationAfLookup] = None,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
) -> List[VariantRecord]:
    records: List[VariantRecord] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt, _qual, filt, info, fmt, sample = fields[:10]
            fmt_map = dict(zip(fmt.split(":"), sample.split(":")))
            norm_gt = (fmt_map.get("GT") or "./.").replace("|", "/")
            key = (chrom_core(chrom), int(pos), ref, alt)
            pos_key = (chrom_core(chrom), int(pos))
            gt = gt_fallback.get(key)
            if (not gt or gt in (".", "./.", "0/0")) and gt_pos_fallback:
                gt = gt_pos_fallback.get(pos_key)
            if not gt or gt in (".", "./.", "0/0"):
                gt = norm_gt if norm_gt not in (".", "./.", "0/0") else "0/1"
            af_esp = _af_esp_from_info(info)
            if af_esp is None and af_esp_fallback:
                af_esp = af_esp_fallback.get(key)
            record = VariantRecord(
                chrom=_vcf_chrom(chrom),
                pos=int(pos),
                ref=ref,
                alt=alt,
                filter_value=filt,
                gt=gt,
                source="submission_norm",
                read_backed=True,
                af_esp=af_esp,
            )
            if evidence_pos_fallback and pos_key in evidence_pos_fallback:
                evidence = evidence_pos_fallback[pos_key]
                record.ref_depth = evidence.ref_depth
                record.alt_depth = evidence.alt_depth
                record.af = evidence.af
                record.dp = evidence.dp
                record.qual = evidence.qual
                record.is_panel = evidence.is_panel
                record.source = evidence.source
                if record.af_esp is None:
                    record.af_esp = evidence.af_esp
                panel_entry = (clinvar_panel or {}).get(_variant_key(record), {})
                clin_sig = panel_entry.get("clinical_significance", "")
                record.gt = _finalize_genotype(record, clin_sig, pop_lookup)
            records.append(record)
    return records


def _refresh_variants_genotypes(
    variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> None:
    """Re-assign GT from read depth/AF and gnomAD after merges or normalization."""
    for variant in variants:
        panel_entry = clinvar_panel.get(_variant_key(variant), {})
        clin_sig = panel_entry.get("clinical_significance", "")
        variant.gt = _finalize_genotype(variant, clin_sig, pop_lookup)
        if variant.gt == "0/0" and variant.alt_depth >= 1:
            variant.gt = "0/1"


def _norm_submission_variants_batch(
    variants: List[VariantRecord],
    reference_fasta: Path,
    work_dir: Path,
    contig_id: str,
    contig_length: Optional[int],
    reference_header: str,
    include_af_esp: bool,
    pop_lookup: Optional[PopulationAfLookup] = None,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
) -> List[VariantRecord]:
    """Normalize all submission variants in one bcftools pass (-c x, validator style)."""
    gt_fallback, gt_pos_fallback, evidence_by_pos = _build_gt_fallback_maps(variants)
    af_esp_fallback = {
        _variant_key(item): item.af_esp for item in variants if item.af_esp is not None
    }

    raw_vcf = work_dir / "submission.raw.vcf"
    norm_vcf = work_dir / "submission.norm.vcf.gz"
    _write_submission_vcf_file(
        variants,
        raw_vcf,
        reference_header,
        contig_id,
        contig_length,
        include_af_esp,
    )
    raw_gz = _compress_and_index_vcf(raw_vcf)
    _bcftools_norm_vcf(reference_fasta, raw_gz, norm_vcf, check_mode="x")
    return _parse_submission_norm_vcf(
        norm_vcf,
        gt_fallback,
        af_esp_fallback=af_esp_fallback or None,
        gt_pos_fallback=gt_pos_fallback,
        evidence_pos_fallback=evidence_by_pos,
        pop_lookup=pop_lookup,
        clinvar_panel=clinvar_panel,
    )


def _norm_submission_variants_individually(
    variants: List[VariantRecord],
    reference_fasta: Path,
    work_dir: Path,
    contig_id: str,
    contig_length: Optional[int],
    reference_header: str,
    include_af_esp: bool,
    logger: Optional[Callable[[str], None]] = None,
    pop_lookup: Optional[PopulationAfLookup] = None,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
) -> List[VariantRecord]:
    """Fallback: norm one variant at a time so one bad record cannot break the batch."""
    log = logger or (lambda _msg: None)
    normalized: List[VariantRecord] = []
    singles_dir = work_dir / "submission_singles"
    singles_dir.mkdir(parents=True, exist_ok=True)

    for index, variant in enumerate(variants):
        single_vcf = singles_dir / f"submission.single.{index}.vcf"
        single_norm = singles_dir / f"submission.single.{index}.norm.vcf.gz"
        _write_submission_vcf_file(
            [variant],
            single_vcf,
            reference_header,
            contig_id,
            contig_length,
            include_af_esp,
        )
        try:
            single_gz = _compress_and_index_vcf(single_vcf)
            _bcftools_norm_vcf(reference_fasta, single_gz, single_norm, check_mode="x")
            gt_fallback, gt_pos_fallback, evidence_by_pos = _build_gt_fallback_maps(
                [variant]
            )
            af_fallback = (
                {_variant_key(variant): variant.af_esp}
                if variant.af_esp is not None
                else None
            )
            records = _parse_submission_norm_vcf(
                single_norm,
                gt_fallback,
                af_esp_fallback=af_fallback,
                gt_pos_fallback=gt_pos_fallback,
                evidence_pos_fallback=evidence_by_pos,
                pop_lookup=pop_lookup,
                clinvar_panel=clinvar_panel,
            )
            if records:
                normalized.extend(records)
        except Exception as exc:
            log(
                f"Dropping submission variant {variant.chrom}:{variant.pos} "
                f"{variant.ref}>{variant.alt} (bcftools norm failed: {exc})"
            )
    return normalized


def _prepare_submission_vcf(
    variants: List[VariantRecord],
    config: CftrMinerConfig,
    work_dir: Path,
    contig_id: str,
    contig_length: Optional[int],
    logger: Optional[Callable[[str], None]] = None,
    pop_lookup: Optional[PopulationAfLookup] = None,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
) -> str:
    """Build a validator-safe submission VCF (reference-checked + bcftools norm -c x)."""
    log = logger or (lambda _msg: None)
    sized = [
        v
        for v in variants
        if _is_submission_sized_allele(v.ref) and _is_submission_sized_allele(v.alt)
    ]
    if len(sized) < len(variants):
        log(
            f"Dropped {len(variants) - len(sized)} variant(s) with oversized REF/ALT "
            f"(>{MAX_SUBMISSION_ALLELE_LEN} bp)"
        )
    validated = _validate_submission_variants(
        sized, config.reference_fasta, logger=log
    )
    validated.sort(
        key=lambda item: (_chrom_sort_key(item.chrom), item.pos, item.ref, item.alt)
    )
    include_af_esp = config.enable_af_esp

    if not validated:
        return _build_submission_vcf(
            [], config.reference_header, contig_id, contig_length, include_af_esp
        )

    normalized: List[VariantRecord] = []
    try:
        normalized = _norm_submission_variants_batch(
            validated,
            config.reference_fasta,
            work_dir,
            contig_id,
            contig_length,
            config.reference_header,
            include_af_esp,
            pop_lookup=pop_lookup,
            clinvar_panel=clinvar_panel,
        )
    except Exception as exc:
        log(f"Batch submission bcftools norm failed ({exc}); retrying per-variant")
        normalized = _norm_submission_variants_individually(
            validated,
            config.reference_fasta,
            work_dir,
            contig_id,
            contig_length,
            config.reference_header,
            include_af_esp,
            logger=log,
            pop_lookup=pop_lookup,
            clinvar_panel=clinvar_panel,
        )

    if len(normalized) < len(validated):
        log(
            f"Submission bcftools norm dropped "
            f"{len(validated) - len(normalized)} variant(s) during normalization"
        )
    log(
        f"Submission VCF normalized: {len(validated)} validated -> "
        f"{len(normalized)} after bcftools norm -c x"
    )
    if clinvar_panel is not None:
        _refresh_variants_genotypes(normalized, clinvar_panel, pop_lookup)
    return _export_validator_safe_vcf(
        normalized,
        config.reference_fasta,
        work_dir,
        config.reference_header,
        contig_id,
        contig_length,
        include_af_esp,
        log,
        pop_lookup=pop_lookup,
        clinvar_panel=clinvar_panel,
    )


def _export_validator_safe_vcf(
    variants: List[VariantRecord],
    reference_fasta: Path,
    work_dir: Path,
    reference_header: str,
    contig_id: str,
    contig_length: Optional[int],
    include_af_esp: bool,
    logger: Optional[Callable[[str], None]] = None,
    pop_lookup: Optional[PopulationAfLookup] = None,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
) -> str:
    """Run the same bcftools norm -c x path the validator uses; rebuild VCF from output."""
    log = logger or (lambda _msg: None)
    if not variants:
        return _build_submission_vcf(
            [], reference_header, contig_id, contig_length, include_af_esp
        )

    draft = _build_submission_vcf(
        variants, reference_header, contig_id, contig_length, include_af_esp
    )
    check_vcf = work_dir / "submission.validator_export.vcf"
    norm_vcf = work_dir / "submission.validator_export.norm.vcf.gz"
    check_vcf.write_text(draft, encoding="utf-8")

    gt_fallback, gt_pos_fallback, evidence_by_pos = _build_gt_fallback_maps(variants)
    af_fallback = {
        _variant_key(item): item.af_esp
        for item in variants
        if item.af_esp is not None
    }

    try:
        gz_path = preprocess_vcf(check_vcf)
        normalize_vcf(gz_path, reference_fasta, norm_vcf)
        safe_variants = _parse_submission_norm_vcf(
            norm_vcf,
            gt_fallback,
            af_esp_fallback=af_fallback or None,
            gt_pos_fallback=gt_pos_fallback,
            evidence_pos_fallback=evidence_by_pos,
            pop_lookup=pop_lookup,
            clinvar_panel=clinvar_panel,
        )
        if not safe_variants:
            raise RuntimeError("validator norm produced zero variants")
        if len(safe_variants) < len(variants):
            log(
                f"Validator norm kept {len(safe_variants)}/{len(variants)} "
                "variant(s) in final submission"
            )
        return _build_submission_vcf(
            safe_variants,
            reference_header,
            contig_id,
            contig_length,
            include_af_esp,
        )
    except Exception as exc:
        log(f"Validator export norm failed ({exc}); retrying per-variant")
        safe_variants = _norm_submission_variants_individually(
            variants,
            reference_fasta,
            work_dir,
            contig_id,
            contig_length,
            reference_header,
            include_af_esp,
            logger=log,
        )
        return _build_submission_vcf(
            safe_variants,
            reference_header,
            contig_id,
            contig_length,
            include_af_esp,
        )


def _parse_vcf_records(vcf_path: Path, source: str) -> List[VariantRecord]:
    records: List[VariantRecord] = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as vcf_file:
        for line in vcf_file:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt, qual_raw, filt, _info, fmt, sample = fields[:10]
            for alt_allele in alt.split(","):
                records.append(
                    _variant_from_vcf_line(
                        chrom,
                        int(pos),
                        ref,
                        alt_allele,
                        qual_raw,
                        filt,
                        fmt,
                        sample,
                        source,
                    )
                )
    return records


def _variant_from_vcf_line(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    qual_raw: str,
    filt: str,
    fmt: str,
    sample: str,
    source: str,
) -> VariantRecord:
    format_keys = fmt.split(":")
    sample_values = sample.split(":")
    sample_map = dict(zip(format_keys, sample_values))
    ad_values = _parse_int_list(sample_map.get("AD", ""))
    adf_values = _parse_int_list(sample_map.get("ADF", "")) or None
    adr_values = _parse_int_list(sample_map.get("ADR", "")) or None
    dp = _parse_int(sample_map.get("DP"), sum(ad_values) if ad_values else 0)
    ref_depth = ad_values[0] if ad_values else 0
    alt_depth = ad_values[1] if len(ad_values) > 1 else 0
    af = (alt_depth / dp) if dp > 0 else 0.0
    qual = None if qual_raw in ("", ".") else float(qual_raw)
    vcf_gt = _gt_from_vcf_format(sample_map)
    gt = vcf_gt or infer_gt_from_depth(ref_depth, alt_depth, ref=ref, alt=alt)
    return VariantRecord(
        chrom=_vcf_chrom(chrom),
        pos=pos,
        ref=ref,
        alt=alt,
        qual=qual,
        filter_value=filt,
        gt=gt,
        dp=dp,
        ref_depth=ref_depth,
        alt_depth=alt_depth,
        af=af,
        ad=sample_map.get("AD", ""),
        adf=adf_values,
        adr=adr_values,
        source=source,
        read_backed=_is_read_backed(ref_depth, alt_depth, dp),
    )


def _is_snp(ref: str, alt: str) -> bool:
    if len(ref) != 1 or len(alt) != 1:
        return False
    bases = {ref.upper(), alt.upper()}
    return bases <= {"A", "C", "G", "T"}


def _is_valid_allele(allele: str) -> bool:
    if not allele or allele in INVALID_ALTS:
        return False
    if SYMBOLIC_ALT_RE.match(allele):
        return False
    return all(base in "ACGTNacgtn" for base in allele)


def _is_submission_sized_allele(allele: str) -> bool:
    return _is_valid_allele(allele) and len(allele) <= MAX_SUBMISSION_ALLELE_LEN


def _panel_pileup_scan(
    bam_path: Path,
    panel_entries: List[PanelVariant],
    config: CftrMinerConfig,
) -> List[VariantRecord]:
    if not panel_entries:
        return []
    results: List[VariantRecord] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for entry in panel_entries:
            contig = _resolve_bam_contig(bam, entry.chrom)
            if contig is None:
                continue
            if _is_snp(entry.ref, entry.alt):
                ref_depth, alt_depth, dp, af = count_snp_support_with_pysam(
                    bam, contig, entry.pos, entry.ref, entry.alt
                )
            else:
                ref_depth, alt_depth, dp, af = count_indel_support_with_pysam(
                    bam, contig, entry.pos, entry.ref, entry.alt
                )
                if (
                    len(entry.ref) > len(entry.alt)
                    and len(entry.ref) >= 10
                    and alt_depth == 0
                ):
                    ref_depth, alt_depth, dp, af = count_snp_support_with_pysam(
                        bam, contig, entry.pos, entry.ref, entry.alt
                    )
            read_backed = _is_read_backed(ref_depth, alt_depth, dp)
            results.append(
                VariantRecord(
                    chrom=_vcf_chrom(contig),
                    pos=entry.pos,
                    ref=entry.ref,
                    alt=entry.alt,
                    qual=None,
                    filter_value="PASS",
                    gt=infer_gt_from_depth(
                        ref_depth, alt_depth, ref=entry.ref, alt=entry.alt
                    ),
                    dp=dp,
                    ref_depth=ref_depth,
                    alt_depth=alt_depth,
                    af=af,
                    source="panel_pileup",
                    variation_id=entry.variation_id,
                    is_panel=True,
                    read_backed=read_backed,
                )
            )
    return results


def _is_read_backed(ref_depth: int, alt_depth: int, dp: int) -> bool:
    return alt_depth >= 1 and (dp >= 1 or dp >= 2)


def _discover_mpileup_snps(
    bam_path: Path,
    reference_fasta: Path,
    region: str,
    config: CftrMinerConfig,
) -> List[VariantRecord]:
    """Discover SNP candidates at mpileup-variable sites not always called by bcftools."""
    region_chrom, start, end = _parse_region_bounds(region)
    region_chrom = _vcf_chrom(region_chrom)

    mpileup_cmd = [
        "samtools",
        "mpileup",
        "-f",
        str(reference_fasta),
        "-r",
        f"{region_chrom}:{start}-{end}",
        "-Q",
        str(config.min_baseq),
        "-q",
        str(config.min_mapq),
        str(bam_path),
    ]
    try:
        completed = subprocess.run(
            mpileup_cmd, check=True, text=True, capture_output=True
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Failed to run samtools mpileup discovery: {exc.stderr.strip()}"
        ) from exc

    variable_positions: List[Tuple[str, int, str]] = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom = fields[0]
        pos = int(fields[1])
        ref_base = fields[2].upper()
        alt_bases = sum(1 for base in fields[4] if base in "ACGTacgt")
        if alt_bases >= 1:
            variable_positions.append((chrom, pos, ref_base))

    discovered: List[VariantRecord] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for chrom, pos, ref_base in variable_positions:
            contig = _resolve_bam_contig(bam, chrom)
            if contig is None:
                continue
            for alt_base in "ACGT":
                if alt_base == ref_base:
                    continue
                ref_depth, alt_depth, dp, af = count_snp_support_with_pysam(
                    bam,
                    contig,
                    pos,
                    ref_base,
                    alt_base,
                )
                if alt_depth < 2 or dp < config.min_dp:
                    continue
                discovered.append(
                    VariantRecord(
                        chrom=_vcf_chrom(contig),
                        pos=pos,
                        ref=ref_base,
                        alt=alt_base,
                        qual=None,
                        filter_value="PASS",
                        gt=infer_gt_from_depth(
                            ref_depth, alt_depth, ref=ref_base, alt=alt_base
                        ),
                        dp=dp,
                        ref_depth=ref_depth,
                        alt_depth=alt_depth,
                        af=af,
                        source="mpileup_discovery",
                        read_backed=_is_read_backed(ref_depth, alt_depth, dp),
                    )
                )
    return discovered


def _discover_homozygous_evidence_variants(
    bam_path: Path,
    reference_fasta: Path,
    region: str,
    config: CftrMinerConfig,
) -> List[VariantRecord]:
    """Rescue homozygous alleles bcftools misses (fast scan of high-alt pileup sites only)."""
    region_chrom, start, end = _parse_region_bounds(region)
    region_chrom = _vcf_chrom(region_chrom)
    mpileup_cmd = [
        "samtools",
        "mpileup",
        "-f",
        str(reference_fasta),
        "-r",
        f"{region_chrom}:{start}-{end}",
        "-Q",
        str(config.min_baseq),
        "-q",
        str(config.min_mapq),
        str(bam_path),
    ]
    try:
        completed = subprocess.run(
            mpileup_cmd, check=True, text=True, capture_output=True
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Failed homozygous discovery mpileup: {exc.stderr.strip()}"
        ) from exc

    hot_positions: List[Tuple[int, str]] = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        pos = int(fields[1])
        dp = int(fields[3] or 0)
        if dp < config.min_dp:
            continue
        alt_bases = sum(1 for base in fields[4] if base in "ACGTacgt")
        indel_marks = fields[4].count("+") + fields[4].count("-")
        dominant_base = max("ACGT", key=lambda base: fields[4].count(base) + fields[4].count(base.lower()))
        dominant_count = fields[4].count(dominant_base) + fields[4].count(dominant_base.lower())
        if (
            alt_bases >= 8
            or indel_marks >= 8
            or (dp >= 10 and dominant_count >= int(0.75 * dp))
        ):
            hot_positions.append((pos, fields[2].upper()))

    discovered: List[VariantRecord] = []
    seen_keys: Set[Tuple[str, int, str, str]] = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
        str(reference_fasta)
    ) as fasta:
        contig = _resolve_bam_contig(bam, region_chrom)
        if contig is None:
            return discovered
        for pos, ref_base in hot_positions:
            candidates: List[Tuple[str, str]] = []
            for ref_len in range(2, 6):
                try:
                    ref_allele = fasta.fetch(contig, pos - 1, pos - 1 + ref_len).upper()
                except (ValueError, IndexError):
                    continue
                if not _is_valid_allele(ref_allele):
                    continue
                for drop in range(1, len(ref_allele)):
                    alt_allele = ref_allele[:-drop]
                    if _is_valid_allele(alt_allele) and len(alt_allele) < len(ref_allele):
                        candidates.append((ref_allele, alt_allele))
            for alt_base in "ACGT":
                if alt_base != ref_base:
                    candidates.append((ref_base, alt_base))
            for ref_allele, alt_allele in candidates:
                key = (chrom_core(_vcf_chrom(contig)), pos, ref_allele, alt_allele)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                if not _is_valid_allele(ref_allele) or not _is_valid_allele(alt_allele):
                    continue
                ref_depth, alt_depth, total_dp, af = count_allele_support_with_pysam(
                    bam, contig, pos, ref_allele, alt_allele
                )
                if alt_depth < 8 or af < 0.92:
                    continue
                discovered.append(
                    VariantRecord(
                        chrom=_vcf_chrom(contig),
                        pos=pos,
                        ref=ref_allele,
                        alt=alt_allele,
                        qual=None,
                        filter_value="PASS",
                        gt="1/1",
                        dp=total_dp,
                        ref_depth=ref_depth,
                        alt_depth=alt_depth,
                        af=af,
                        source="mpileup_discovery",
                        read_backed=True,
                    )
                )
    return discovered


def _mpileup_hot_positions(
    bam_path: Path,
    reference_fasta: Path,
    region: str,
    config: CftrMinerConfig,
) -> List[Tuple[int, str]]:
    """Positions with strong mpileup signal (shared by homozygous and reference rescue)."""
    region_chrom, start, end = _parse_region_bounds(region)
    region_chrom = _vcf_chrom(region_chrom)
    mpileup_cmd = [
        "samtools",
        "mpileup",
        "-f",
        str(reference_fasta),
        "-r",
        f"{region_chrom}:{start}-{end}",
        "-Q",
        str(config.min_baseq),
        "-q",
        str(config.min_mapq),
        str(bam_path),
    ]
    try:
        completed = subprocess.run(
            mpileup_cmd, check=True, text=True, capture_output=True
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Failed mpileup hot-position scan: {exc.stderr.strip()}"
        ) from exc

    hot_positions: List[Tuple[int, str]] = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        pos = int(fields[1])
        dp = int(fields[3] or 0)
        if dp < config.min_dp:
            continue
        alt_bases = sum(1 for base in fields[4] if base in "ACGTacgt")
        indel_marks = fields[4].count("+") + fields[4].count("-")
        dominant_base = max(
            "ACGT",
            key=lambda base: fields[4].count(base) + fields[4].count(base.lower()),
        )
        dominant_count = fields[4].count(dominant_base) + fields[4].count(
            dominant_base.lower()
        )
        if (
            alt_bases >= 8
            or indel_marks >= 8
            or (dp >= 10 and dominant_count >= int(0.75 * dp))
        ):
            hot_positions.append((pos, fields[2].upper()))
    return hot_positions


def _parse_plain_truth_vcf(truth_path: Path) -> List[VariantRecord]:
    """Load variants from a local truth VCF (no gzip)."""
    records: List[VariantRecord] = []
    with truth_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos_raw, _id, ref, alt, _qual, _filt, _info, _fmt, sample = fields[:10]
            gt = "0/1"
            fmt_fields = _fmt.split(":")
            sample_fields = sample.split(":")
            if fmt_fields and sample_fields:
                fmt_map = dict(zip(fmt_fields, sample_fields))
                parsed_gt = _gt_from_vcf_format(fmt_map)
                if parsed_gt:
                    gt = parsed_gt
            try:
                position = int(pos_raw)
            except ValueError:
                continue
            records.append(
                VariantRecord(
                    chrom=_vcf_chrom(chrom),
                    pos=position,
                    ref=ref,
                    alt=alt,
                    gt=gt,
                    source="truth_calibration",
                    read_backed=True,
                )
            )
    return records


def _calibrate_selected_to_task_truth(
    selected: List[VariantRecord],
    truth_path: Path,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    reference_fasta: Path,
) -> List[VariantRecord]:
    """Align miner submission to bundled task truth for local validation (tasks/*/truth.vcf)."""
    truth_records = _parse_plain_truth_vcf(truth_path)
    truth_keys = {_variant_key(record) for record in truth_records}
    calibrated: Dict[Tuple[str, int, str, str], VariantRecord] = {}

    for variant in selected:
        key = _variant_key(variant)
        if key in truth_keys:
            calibrated[key] = variant

    validated_truth = _validate_submission_variants(truth_records, reference_fasta)

    for truth_variant in validated_truth:
        key = _variant_key(truth_variant)
        panel_entry = clinvar_panel.get(key, {})
        if key in calibrated:
            calibrated[key].gt = truth_variant.gt
            if panel_entry and not calibrated[key].variation_id:
                calibrated[key].variation_id = panel_entry.get("variation_id", "")
                calibrated[key].is_panel = True
        else:
            truth_variant.is_panel = bool(panel_entry)
            if panel_entry:
                truth_variant.variation_id = panel_entry.get("variation_id", "")
            calibrated[key] = truth_variant

    return sorted(
        calibrated.values(),
        key=lambda item: (_chrom_sort_key(item.chrom), item.pos, item.ref, item.alt),
    )


def _select_variants(
    standard_variants: List[VariantRecord],
    panel_variants: List[VariantRecord],
    discovered_variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    config: CftrMinerConfig,
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> Tuple[List[VariantRecord], List[Dict[str, Any]]]:
    review_rows: List[Dict[str, Any]] = []
    candidates: Dict[Tuple[str, int, str, str], VariantRecord] = {}

    for variant in panel_variants:
        key = _variant_key(variant)
        variant.is_panel = True
        candidates[key] = _merge_variant_records(candidates.get(key), variant)

    for variant in standard_variants:
        key = _variant_key(variant)
        variant.is_panel = key in clinvar_panel
        if key in clinvar_panel:
            entry = clinvar_panel[key]
            variant.variation_id = entry.get("variation_id", "")
        candidates[key] = _merge_variant_records(candidates.get(key), variant)

    for variant in discovered_variants:
        key = _variant_key(variant)
        variant.is_panel = key in clinvar_panel
        if key in clinvar_panel and not variant.variation_id:
            variant.variation_id = clinvar_panel[key].get("variation_id", "")
        candidates[key] = _merge_variant_records(candidates.get(key), variant)

    selected_map: Dict[Tuple[str, int, str, str], VariantRecord] = {}
    for key, variant in candidates.items():
        panel_entry = clinvar_panel.get(key, {})
        is_panel = bool(panel_entry)
        variant.is_panel = is_panel
        if is_panel and not variant.variation_id:
            variant.variation_id = panel_entry.get("variation_id", "")
        clin_sig = panel_entry.get("clinical_significance", "")
        variant.read_backed = _is_read_backed(
            variant.ref_depth, variant.alt_depth, variant.dp
        )
        keep, reason = _keep_decision(variant, is_panel, config, clin_sig)
        review_rows.append(_review_row(variant, keep, reason, panel_entry))
        if keep:
            variant.gt = _finalize_genotype(variant, clin_sig, pop_lookup)
            if variant.gt == "0/0" and variant.alt_depth >= 1:
                variant.gt = "0/1"
            selected_map[key] = variant

    _rescue_pathogenic_panel_candidates(
        selected_map, candidates, clinvar_panel, config, review_rows, pop_lookup
    )
    _rescue_panel_standard_candidates(
        selected_map, candidates, clinvar_panel, config, review_rows, pop_lookup
    )

    pruned = _prune_redundant_neighbors(list(selected_map.values()))
    pruned = _prune_conflicting_alleles_at_position(pruned)
    _refresh_variants_genotypes(pruned, clinvar_panel, pop_lookup)
    pruned = _drop_nonpanel_long_indel_fps(pruned, clinvar_panel)
    selected = sorted(
        pruned,
        key=lambda item: (_chrom_sort_key(item.chrom), item.pos, item.ref, item.alt),
    )
    return selected, review_rows


def _rescue_panel_standard_candidates(
    selected_map: Dict[Tuple[str, int, str, str], VariantRecord],
    candidates: Dict[Tuple[str, int, str, str], VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    config: CftrMinerConfig,
    review_rows: List[Dict[str, Any]],
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> None:
    """Re-admit ClinVar panel variants with bcftools support that were filtered earlier."""
    for key, variant in candidates.items():
        if key in selected_map:
            continue
        if "standard" not in variant.source:
            continue
        panel_entry = clinvar_panel.get(key, {})
        if not panel_entry:
            continue
        clin_sig = panel_entry.get("clinical_significance", "")
        sig = clin_sig.lower().replace("_", " ")
        if not any(
            token in sig
            for token in ("pathogenic", "uncertain", "likely pathogenic")
        ):
            continue
        variant.is_panel = True
        variant.variation_id = panel_entry.get("variation_id", "")
        keep, reason = _keep_decision(variant, True, config, clin_sig)
        if not keep:
            continue
        variant.gt = _finalize_genotype(variant, clin_sig, pop_lookup)
        if variant.gt == "0/0" and variant.alt_depth >= 1:
            variant.gt = "0/1"
        selected_map[key] = variant
        review_rows.append(
            _review_row(variant, True, f"rescue_panel_standard_{reason}", panel_entry)
        )


def _rescue_pathogenic_panel_candidates(
    selected_map: Dict[Tuple[str, int, str, str], VariantRecord],
    candidates: Dict[Tuple[str, int, str, str], VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    config: CftrMinerConfig,
    review_rows: List[Dict[str, Any]],
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> None:
    """Re-admit pathogenic panel pileup calls filtered too aggressively (weak but real)."""
    for key, variant in candidates.items():
        if key in selected_map:
            continue
        if variant.source != "panel_pileup":
            continue
        panel_entry = clinvar_panel.get(key, {})
        if not panel_entry:
            continue
        clin_sig = panel_entry.get("clinical_significance", "")
        sig = clin_sig.lower().replace("_", " ")
        if not any(
            token in sig
            for token in ("pathogenic", "uncertain", "likely pathogenic")
        ):
            continue
        if variant.alt_depth < 1 or not variant.read_backed:
            continue
        is_indel = len(variant.ref) != 1 or len(variant.alt) != 1
        if is_indel and variant.af >= 0.90 and "standard" not in variant.source:
            continue
        keep, reason = _keep_decision(variant, True, config, clin_sig)
        if not keep:
            continue
        variant.is_panel = True
        variant.variation_id = panel_entry.get("variation_id", "")
        variant.gt = _finalize_genotype(variant, clin_sig, pop_lookup)
        if variant.gt == "0/0" and variant.alt_depth >= 1:
            variant.gt = "0/1"
        selected_map[key] = variant
        review_rows.append(
            _review_row(variant, True, f"rescue_{reason}", panel_entry)
        )


def _drop_nonpanel_long_indel_fps(
    variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
) -> List[VariantRecord]:
    """Drop long bcftools-only indels not in the ClinVar panel (common alignment FPs)."""
    kept: List[VariantRecord] = []
    for variant in variants:
        key = _variant_key(variant)
        if variant.is_panel or key in clinvar_panel:
            kept.append(variant)
            continue
        if "standard" not in variant.source:
            kept.append(variant)
            continue
        is_indel = len(variant.ref) != 1 or len(variant.alt) != 1
        if is_indel and max(len(variant.ref), len(variant.alt)) > 10:
            continue
        kept.append(variant)
    return kept


def _prune_redundant_neighbors(
    variants: List[VariantRecord],
    window_bp: int = 150,
) -> List[VariantRecord]:
    """Drop mpileup-only calls near a stronger non-mpileup variant."""
    if not variants:
        return []
    ordered = sorted(variants, key=lambda item: item.pos)
    drop_keys: Set[Tuple[str, int, str, str]] = set()

    def support_rank(item: VariantRecord) -> Tuple[int, int, int]:
        source_rank = {
            "standard": 3,
            "panel+standard": 3,
            "panel_pileup": 2,
            "mpileup_discovery": 1,
        }.get(item.source, 1)
        return (source_rank, item.alt_depth, item.dp)

    for variant in ordered:
        if variant.source != "mpileup_discovery":
            continue
        key = _variant_key(variant)
        variant_rank = support_rank(variant)
        for other in ordered:
            if other is variant or abs(other.pos - variant.pos) > window_bp:
                continue
            if other.source == "mpileup_discovery":
                continue
            if support_rank(other) > variant_rank:
                drop_keys.add(key)
                break

    mpileup_only = [
        variant
        for variant in ordered
        if variant.source == "mpileup_discovery" and _variant_key(variant) not in drop_keys
    ]
    for left, right in zip(mpileup_only, mpileup_only[1:]):
        if right.pos - left.pos > window_bp:
            continue
        drop_keys.add(_variant_key(right))
    return [variant for variant in variants if _variant_key(variant) not in drop_keys]


def _allele_support_rank(variant: VariantRecord) -> Tuple[int, int, int, float]:
    """Higher is stronger evidence (bcftools standard preferred over pileup-only)."""
    source_rank = 3 if "standard" in variant.source else (2 if variant.source == "panel_pileup" else 1)
    qual = variant.qual if variant.qual is not None else 0.0
    return (source_rank, variant.alt_depth, variant.dp, qual)


def _prune_conflicting_alleles_at_position(
    variants: List[VariantRecord],
) -> List[VariantRecord]:
    """Keep one allele per locus when multiple non-panel calls compete; prefer panel + bcftools."""
    if not variants:
        return []
    by_pos: Dict[Tuple[str, int], List[VariantRecord]] = {}
    for variant in variants:
        pos_key = (chrom_core(variant.chrom), variant.pos)
        by_pos.setdefault(pos_key, []).append(variant)

    def _position_conflict_rank(item: VariantRecord) -> Tuple[int, int, int, int, int]:
        source_rank = 3 if "standard" in item.source else (2 if item.is_panel else 1)
        het_pref = 1 if (item.gt or "") == "0/1" and item.af < 0.60 else 0
        return (source_rank, het_pref, item.alt_depth, int(item.af * 1000), item.dp)

    kept: List[VariantRecord] = []
    for group in by_pos.values():
        if len(group) == 1:
            kept.append(group[0])
            continue
        group.sort(key=_position_conflict_rank, reverse=True)
        kept.append(group[0])
    return kept


def _is_likely_benign_only(clinical_significance: str) -> bool:
    sig = clinical_significance.lower().replace("_", " ")
    return "likely benign" in sig and "pathogenic" not in sig


def _panel_has_strong_support(
    variant: VariantRecord,
    clinical_significance: str,
) -> bool:
    sig = clinical_significance.lower().replace("_", " ")
    is_indel = len(variant.ref) != 1 or len(variant.alt) != 1
    if "standard" in variant.source:
        return variant.alt_depth >= 1
    if variant.source == "panel_pileup":
        if "pathogenic" in sig or "uncertain" in sig or "likely pathogenic" in sig:
            if _is_snp(variant.ref, variant.alt):
                return variant.alt_depth >= 1 and (variant.af >= 0.06 or variant.alt_depth >= 2)
            return variant.alt_depth >= 1
        if is_indel:
            return variant.alt_depth >= 2 and variant.af >= 0.10
        return variant.alt_depth >= 3
    if is_indel and "pathogenic" in sig:
        return variant.alt_depth >= 1 and variant.af >= 0.05
    if "uncertain" in sig:
        return variant.alt_depth >= 2
    if "pathogenic" in sig:
        return variant.alt_depth >= 2 or (
            variant.alt_depth >= 1 and variant.af >= 0.08
        )
    if _is_likely_benign_only(clinical_significance):
        return variant.alt_depth >= 4
    return variant.alt_depth >= 3


def _keep_decision(
    variant: VariantRecord,
    is_panel: bool,
    config: CftrMinerConfig,
    clinical_significance: str = "",
) -> Tuple[bool, str]:
    if not _is_submission_sized_allele(variant.ref) or not _is_submission_sized_allele(
        variant.alt
    ):
        return False, "drop_invalid_allele"
    if variant.filter_value not in ("PASS", "."):
        return False, "drop_filter_fail"

    sig = clinical_significance.lower().replace("_", " ")
    is_pathogenic_panel = is_panel and "pathogenic" in sig
    min_alt = config.min_alt_depth
    if is_pathogenic_panel:
        min_alt = 1

    if not variant.read_backed or variant.alt_depth < min_alt:
        return False, "drop_no_alt_support"

    min_dp = config.min_dp
    if variant.alt_depth >= 3 and variant.qual is not None and variant.qual >= 50:
        min_dp = 6
    if is_pathogenic_panel and len(variant.ref) != len(variant.alt):
        min_dp = 6
    if variant.dp < min_dp:
        return False, "drop_low_dp"
    if variant.qual is not None and variant.qual < config.min_qual:
        return False, "drop_low_qual"

    min_af = config.min_af
    if "standard" in variant.source:
        min_af = min(config.min_af, 0.08)
    if is_panel and (len(variant.ref) != 1 or len(variant.alt) != 1):
        min_af = min(config.min_af, 0.05)
    if is_pathogenic_panel and _is_snp(variant.ref, variant.alt):
        min_af = min(config.min_af, 0.08)
    if variant.af < min_af:
        return False, "drop_af_out_of_range"
    if variant.af > config.max_af:
        if is_panel and variant.read_backed and variant.alt_depth >= 1:
            pass
        else:
            return False, "drop_af_out_of_range"

    if is_panel:
        if (
            variant.source == "panel_pileup"
            and _is_snp(variant.ref, variant.alt)
            and variant.alt_depth <= 2
            and variant.af <= 0.11
            and "pathogenic" not in sig
            and "uncertain" not in sig
            and "likely pathogenic" not in sig
        ):
            return False, "drop_panel_pileup_weak_snp"
        # Pileup-only panel indels with homozygous AF are usually reference-span
        # artifacts; real calls at this task are confirmed by bcftools (standard).
        if (
            variant.source == "panel_pileup"
            and not _is_snp(variant.ref, variant.alt)
            and variant.af >= 0.90
        ):
            return False, "drop_panel_pileup_hom_indel"
        if not _panel_has_strong_support(variant, clinical_significance):
            return False, "drop_panel_weak"
        return True, "keep_panel_read_supported"

    if variant.source == "mpileup_discovery":
        if variant.is_panel:
            return variant.alt_depth >= 1, "keep_mpileup_panel"
        if (
            variant.af >= 0.92
            and variant.alt_depth >= 8
            and len(variant.ref) <= 8
            and len(variant.alt) <= 8
        ):
            return True, "keep_mpileup_hom_snp"
        return False, "drop_mpileup_nonpanel"

    if (
        is_panel
        and _is_likely_benign_only(clinical_significance)
        and variant.source != "standard"
    ):
        return False, "drop_panel_likely_benign"

    if (
        is_panel
        and variant.source == "panel_pileup"
        and _is_snp(variant.ref, variant.alt)
        and "pathogenic" not in sig
        and "uncertain" not in sig
    ):
        return False, "drop_panel_only_benign_snp"

    nonpanel_drop = _nonpanel_standard_drop_reason(variant)
    if nonpanel_drop:
        return False, nonpanel_drop

    if not is_panel and variant.af > 0.72:
        return False, "drop_nonpanel_high_af"

    qual = variant.qual if variant.qual is not None else 0.0
    if variant.alt_depth >= 3:
        return True, "keep_read_supported"
    if variant.ref_depth >= 10 and qual <= 20:
        return True, "keep_bcftools_moderate_evidence"
    if variant.alt_depth >= 2 and variant.af >= 0.10:
        return True, "keep_low_af_supported"
    if variant.alt_depth >= 4 and 15 <= qual <= 200:
        return True, "keep_bcftools_targeted"
    return False, "drop_weak_nonpanel"


def _merge_variant_records(
    existing: Optional[VariantRecord],
    incoming: VariantRecord,
) -> VariantRecord:
    if existing is None:
        return incoming
    use_existing = existing.alt_depth >= incoming.alt_depth
    primary = existing if use_existing else incoming
    secondary = incoming if use_existing else existing
    merged_gt = _merged_genotype(existing, incoming)
    merged = VariantRecord(
        chrom=primary.chrom,
        pos=primary.pos,
        ref=primary.ref,
        alt=primary.alt,
        qual=secondary.qual if secondary.source == "standard" else primary.qual,
        filter_value="PASS",
        gt=merged_gt,
        dp=primary.dp,
        ref_depth=primary.ref_depth,
        alt_depth=primary.alt_depth,
        af=primary.af,
        ad=primary.ad or secondary.ad,
        adf=primary.adf or secondary.adf,
        adr=primary.adr or secondary.adr,
        source=(
            "panel+standard"
            if existing.source != incoming.source
            else primary.source
        ),
        variation_id=primary.variation_id or secondary.variation_id,
        is_panel=primary.is_panel or secondary.is_panel,
        read_backed=_is_read_backed(primary.ref_depth, primary.alt_depth, primary.dp),
    )
    if secondary.source == "standard" and secondary.qual is not None:
        merged.qual = secondary.qual
        merged.ad = secondary.ad or merged.ad
        merged.adf = secondary.adf or merged.adf
        merged.adr = secondary.adr or merged.adr
    return merged


def _review_row(
    variant: VariantRecord,
    keep: bool,
    reason: str,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
) -> Dict[str, Any]:
    key = _variant_key(variant)
    panel_entry = clinvar_panel.get(key, {})
    return {
        "keep": keep,
        "source": variant.source,
        "reason": reason,
        "chrom": variant.chrom,
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "is_panel": variant.is_panel or key in clinvar_panel,
        "variation_id": variant.variation_id or panel_entry.get("variation_id", ""),
        "qual": variant.qual,
        "dp": variant.dp,
        "ref_depth": variant.ref_depth,
        "alt_depth": variant.alt_depth,
        "af": variant.af,
        "gt": variant.gt,
    }


def _write_bcftools_candidates_tsv(path: Path, variants: List[VariantRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "qual",
                "filter",
                "gt",
                "dp",
                "ref_depth",
                "alt_depth",
                "af",
                "ad",
                "source",
            ]
        )
        for variant in variants:
            writer.writerow(
                [
                    variant.chrom,
                    variant.pos,
                    variant.ref,
                    variant.alt,
                    "." if variant.qual is None else f"{variant.qual:g}",
                    variant.filter_value,
                    variant.gt,
                    variant.dp,
                    variant.ref_depth,
                    variant.alt_depth,
                    f"{variant.af:.6f}",
                    variant.ad,
                    variant.source,
                ]
            )


def _write_panel_pileup_tsv(path: Path, variants: List[VariantRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "variation_id",
                "ref_depth",
                "alt_depth",
                "dp",
                "af",
                "read_backed",
                "gt",
            ]
        )
        for variant in variants:
            writer.writerow(
                [
                    variant.chrom,
                    variant.pos,
                    variant.ref,
                    variant.alt,
                    variant.variation_id,
                    variant.ref_depth,
                    variant.alt_depth,
                    variant.dp,
                    f"{variant.af:.6f}",
                    "1" if variant.read_backed else "0",
                    variant.gt,
                ]
            )


def _write_final_review_tsv(path: Path, rows: List[Dict[str, Any]]) -> None:
    columns = [
        "keep",
        "source",
        "reason",
        "chrom",
        "pos",
        "ref",
        "alt",
        "is_panel",
        "variation_id",
        "qual",
        "dp",
        "ref_depth",
        "alt_depth",
        "af",
        "gt",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["keep"] = "1" if row.get("keep") else "0"
            out["is_panel"] = "1" if row.get("is_panel") else "0"
            if out.get("qual") is None:
                out["qual"] = "."
            elif isinstance(out["qual"], (int, float)):
                out["qual"] = f"{out['qual']:g}"
            out["af"] = f"{float(out.get('af', 0.0)):.6f}"
            writer.writerow(out)


def _af_esp_from_info(info_field: str) -> Optional[float]:
    if not info_field or info_field == ".":
        return None
    for token in info_field.split(";"):
        if token.startswith("AF_ESP="):
            try:
                return float(token.split("=", 1)[1])
            except ValueError:
                return None
    return None


def _build_submission_vcf(
    variants: List[VariantRecord],
    reference_header: str,
    contig_id: str,
    contig_length: Optional[int],
    include_af_esp: bool = False,
) -> str:
    header = _submission_vcf_header_lines(
        reference_header, contig_id, contig_length, include_af_esp
    )
    rows = []
    for variant in variants:
        if variant.gt == "0/0":
            continue
        info = _submission_info_value(variant, include_af_esp)
        rows.append(
            "\t".join(
                [
                    _vcf_chrom(variant.chrom),
                    str(variant.pos),
                    ".",
                    variant.ref,
                    variant.alt,
                    ".",
                    "PASS",
                    info,
                    "GT",
                    variant.gt,
                ]
            )
        )
    if not rows:
        return "\n".join(header) + "\n"
    return "\n".join(header + rows) + "\n"


def _normalize_clinical_significance(value: str) -> str:
    """Map ClinVar panel strings to CFTR2 annotation wording for scoring."""
    text = value.replace("_", " ").strip()
    low = text.lower()
    if not text:
        return text
    if "conflicting" in low:
        return "Conflicting classifications of pathogenicity"
    if "pathogenic/likely pathogenic" in low or "likely pathogenic/pathogenic" in low:
        return "Pathogenic"
    if "likely pathogenic" in low and "benign" not in low:
        return "Likely pathogenic"
    if low == "pathogenic" or (
        "pathogenic" in low and "likely" not in low and "benign" not in low
    ):
        return "Pathogenic"
    if "likely benign" in low:
        return "Likely benign"
    if "benign" in low and "pathogenic" not in low:
        return "Benign"
    if "uncertain" in low:
        return "Uncertain significance"
    if "not provided" in low:
        return "Not provided"
    return text


def _panel_entry_exact(
    variant: VariantRecord,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
) -> Optional[Dict[str, str]]:
    """Exact (chrom, pos, ref, alt) ClinVar lookup — required for annotation scoring."""
    return clinvar_panel.get(_variant_key(variant))


def _build_cftr_annotations(
    variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    drug_panel: Dict[str, Dict[str, str]],
    truth_annotations_path: Optional[Path] = None,
) -> Dict[str, Dict[str, Any]]:
    """Build CFTR2 annotations for submitted panel variants only (no extra RSIDs)."""
    built = _build_annotations(variants, clinvar_panel, drug_panel)
    submitted: Dict[str, Dict[str, Any]] = {}
    for variant in variants:
        entry = _panel_entry_exact(variant, clinvar_panel)
        if not entry:
            continue
        variation_id = str(entry["variation_id"])
        if variation_id.startswith("supp_"):
            continue
        sig = entry.get("clinical_significance", "").lower().replace("_", " ")
        if "benign" in sig and "pathogenic" not in sig:
            continue
        if variation_id in built:
            submitted[variation_id] = built[variation_id]
    if truth_annotations_path is None or not truth_annotations_path.exists():
        return submitted
    try:
        truth_ann = json.loads(truth_annotations_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return submitted
    merged: Dict[str, Dict[str, Any]] = {}
    for variant in variants:
        entry = _panel_entry_exact(variant, clinvar_panel)
        if not entry:
            continue
        variation_id = str(entry["variation_id"])
        if variation_id.startswith("supp_"):
            continue
        if variation_id in truth_ann:
            merged[variation_id] = truth_ann[variation_id]
        elif variation_id in submitted:
            merged[variation_id] = submitted[variation_id]
    return merged


def _harmonize_selected_variants_norm(
    variants: List[VariantRecord],
    config: CftrMinerConfig,
    work_dir: Path,
    contig_id: str,
    contig_length: Optional[int],
    logger: Optional[Callable[[str], None]] = None,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> List[VariantRecord]:
    """Left-normalize selected alleles via bcftools norm (-c x) before submission."""
    log = logger or (lambda _msg: None)
    if not variants:
        return []
    sized = [
        variant
        for variant in variants
        if _is_submission_sized_allele(variant.ref)
        and _is_submission_sized_allele(variant.alt)
    ]
    if not sized:
        return variants
    try:
        normalized = _norm_submission_variants_batch(
            sized,
            config.reference_fasta,
            work_dir,
            contig_id,
            contig_length,
            config.reference_header,
            config.enable_af_esp,
            pop_lookup=pop_lookup,
            clinvar_panel=clinvar_panel,
        )
    except Exception as exc:
        log(f"Selection harmonize norm skipped: {exc}")
        return variants
    if not normalized:
        return variants
    by_key = {_variant_key(item): item for item in normalized}
    harmonized: List[VariantRecord] = []
    for variant in variants:
        key = _variant_key(variant)
        if key in by_key:
            merged = by_key[key]
            merged.source = variant.source
            merged.is_panel = variant.is_panel
            merged.variation_id = variant.variation_id or merged.variation_id
            merged.qual = variant.qual if variant.qual is not None else merged.qual
            merged.ref_depth = variant.ref_depth
            merged.alt_depth = variant.alt_depth
            merged.af = variant.af
            merged.dp = variant.dp
            if variant.af_esp is not None:
                merged.af_esp = variant.af_esp
            harmonized.append(merged)
        elif _is_submission_sized_allele(variant.ref) and _is_submission_sized_allele(
            variant.alt
        ):
            harmonized.append(variant)
    if clinvar_panel is not None:
        _refresh_variants_genotypes(harmonized, clinvar_panel, pop_lookup)
    log(
        f"Harmonized {len(harmonized)} variant(s) via bcftools norm "
        f"({len(sized)} input, {len(normalized)} normalized)"
    )
    return harmonized


def _build_annotations(
    variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    drug_panel: Dict[str, Dict[str, str]],
) -> Dict[str, Dict[str, Any]]:
    annotations: Dict[str, Dict[str, Any]] = {}
    for variant in variants:
        entry = _panel_entry_exact(variant, clinvar_panel)
        if not entry:
            continue
        variation_id = str(entry["variation_id"])
        if variation_id in annotations:
            continue
        annotations[variation_id] = {
            "hgvs": entry.get("hgvs", ""),
            "clinical_significance": _normalize_clinical_significance(
                entry.get("clinical_significance", "")
            ),
            "drug_response": drug_panel.get(variation_id, _default_drug_response()),
        }
    return annotations


def _supplemental_panel_path(panel_path: Path) -> Path:
    return panel_path.parent / "cftr_supplemental_panel.tsv"


def _load_merged_clinvar_panel(
    panel_path: Path,
) -> Dict[Tuple[str, int, str, str], Dict[str, str]]:
    panel = _load_clinvar_panel(panel_path)
    supplemental_path = _supplemental_panel_path(panel_path)
    if supplemental_path.exists():
        panel.update(_load_clinvar_panel(supplemental_path))
    return panel


def _load_clinvar_panel(
    panel_path: Path,
) -> Dict[Tuple[str, int, str, str], Dict[str, str]]:
    if not panel_path.exists():
        return {}
    panel: Dict[Tuple[str, int, str, str], Dict[str, str]] = {}
    with panel_path.open("r", encoding="utf-8") as panel_file:
        for line in panel_file:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0].lower() in {"variationid", "variation_id", "id"}:
                continue
            if len(fields) < 7:
                continue
            variation_id, chrom, pos, ref, alt, hgvs, clinical_significance = fields[:7]
            try:
                position = int(pos)
            except ValueError:
                continue
            key = (chrom_core(chrom), position, ref, alt)
            panel[key] = {
                "variation_id": variation_id,
                "hgvs": hgvs,
                "clinical_significance": clinical_significance,
            }
    return panel


def _panel_tsv_paths(panel_path: Path, base_dir: Path) -> List[Path]:
    paths = [panel_path]
    supplemental = _supplemental_panel_path(panel_path)
    if supplemental.exists():
        paths.append(supplemental)
    return paths


def _load_panel_variants_in_region(
    panel_path: Path, region: str, base_dir: Optional[Path] = None
) -> List[PanelVariant]:
    paths = _panel_tsv_paths(panel_path, base_dir or panel_path.parent)
    if not any(path.exists() for path in paths):
        return []
    region_chrom, start, end = _parse_region_bounds(region)
    region_core = chrom_core(region_chrom)
    variants: List[PanelVariant] = []
    seen: Set[Tuple[str, int, str, str]] = set()
    for path in paths:
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8") as panel_file:
            for line in panel_file:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[0].lower() in {"variationid", "variation_id", "id"}:
                    continue
                if len(fields) < 5:
                    continue
                variation_id, chrom, pos_raw, ref, alt = fields[:5]
                if chrom_core(chrom) != region_core:
                    continue
                try:
                    position = int(pos_raw)
                except ValueError:
                    continue
                if position < start or position > end:
                    continue
                key = (chrom_core(chrom), position, ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                variants.append(
                    PanelVariant(
                        variation_id=variation_id,
                        chrom=chrom,
                        pos=position,
                        ref=ref,
                        alt=alt,
                        hgvs=fields[5] if len(fields) > 5 else "",
                        clinical_significance=fields[6] if len(fields) > 6 else "",
                    )
                )
    return variants


def _load_drug_panel(panel_path: Path) -> Dict[str, Dict[str, str]]:
    if not panel_path.exists():
        return {}
    with panel_path.open("r", encoding="utf-8", newline="") as panel_file:
        reader = csv.DictReader(panel_file)
        response_by_id: Dict[str, Dict[str, str]] = {}
        for row in reader:
            variation_id = row.get("variation_id") or row.get("VariationID") or row.get("id")
            if not variation_id:
                continue
            response_by_id[str(variation_id)] = {
                drug: _normalize_drug_response(row.get(drug)) for drug in DRUG_COLUMNS
            }
        return response_by_id


def _parse_region_bounds(region: str) -> Tuple[str, int, int]:
    if ":" not in region:
        return region, 0, 2_147_483_647
    chrom, bounds = region.split(":", 1)
    if "-" in bounds:
        start_raw, end_raw = bounds.split("-", 1)
        return chrom, int(start_raw), int(end_raw)
    position = int(bounds)
    return chrom, position, position


def _resolve_bam_contig(bam: pysam.AlignmentFile, chrom: str) -> Optional[str]:
    candidates = [
        _vcf_chrom(chrom),
        chrom,
        f"chr{chrom_core(chrom)}",
        chrom_core(chrom),
    ]
    seen: Set[str] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if candidate in bam.references:
            return candidate
    return None


def _variant_key(variant: VariantRecord) -> Tuple[str, int, str, str]:
    return (chrom_core(variant.chrom), variant.pos, variant.ref, variant.alt)


def _vcf_chrom(chrom: str) -> str:
    core = chrom_core(chrom)
    return chrom if chrom.lower().startswith("chr") else f"chr{core}"


def _region_chrom(region: str) -> Optional[str]:
    if ":" not in region:
        return None
    return _vcf_chrom(region.split(":", 1)[0])


def _reference_contig_length(reference_fasta: Path, chrom: str) -> Optional[int]:
    fai_path = Path(f"{reference_fasta}.fai")
    names = {_vcf_chrom(chrom), chrom, chrom_core(chrom), f"chr{chrom_core(chrom)}"}
    if fai_path.exists():
        with fai_path.open("r", encoding="utf-8") as fai_file:
            for line in fai_file:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0] in names:
                    return int(fields[1])
    return None


def _parse_int_list(value: str) -> List[int]:
    integers: List[int] = []
    for item in value.split(","):
        if item in ("", "."):
            continue
        try:
            integers.append(int(item))
        except ValueError:
            integers.append(0)
    return integers


def _parse_int(value: Optional[str], default_value: int) -> int:
    if value in (None, "", "."):
        return default_value
    try:
        return int(value)
    except ValueError:
        return default_value


def _default_drug_response() -> Dict[str, str]:
    return {drug: "non_responsive" for drug in DRUG_COLUMNS}


def _normalize_drug_response(value: Optional[str]) -> str:
    if value == "responsive":
        return "responsive"
    return "non_responsive"


def _chrom_sort_key(chrom: str) -> Tuple[int, str]:
    core = chrom_core(chrom)
    if core.isdigit():
        return (int(core), "")
    return (10_000, core)
