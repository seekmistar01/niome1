"""CFTR variant-calling workflow used by the subnet 55 miner.

Read-evidence-first pipeline: FASTQ -> BAM -> GATK + pathogenic panel candidates ->
optional sensitive multi-caller strategy (bcftools, FreeBayes, weakscan; see
``NIOME_USE_STRATEGY_AUTO``) -> normalize -> pileup evidence ->
hard filter -> genotype from AD/AF only -> rank/select -> submission.

gnomAD is used only for optional INFO (AF_ESP) and tie-breaking during ranking;
it never creates variants or sets genotypes.
"""

from __future__ import annotations

import csv
import fcntl
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
from niome_subnet.genomics.vcf_norm import (
    normalize_vcf,
    preprocess_vcf,
    verify_vcf_for_validator_scoring,
)
from niome_subnet import strategy_auto


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
    source_gatk: bool = False
    source_bcftools: bool = False
    source_assembly: bool = False
    source_panel_probe: bool = False
    source_panel_rescue: bool = False
    evidence_score: float = 0.0
    evidence_status: str = ""
    evidence_reason: str = ""
    ref_match: bool = True
    strand_balanced: bool = False
    mean_alt_bq: float = 0.0
    mean_mapq: float = 0.0


# Evidence-first hard thresholds (override via env where noted).
EVIDENCE_MIN_DP = 3
EVIDENCE_MIN_ALT_DEPTH = 2
EVIDENCE_MIN_AF = 0.12

# Sensitive multi-caller strategy (niome_subnet/strategy_auto.py) — NIOME allele imbalance.
def _use_strategy_auto() -> bool:
    return os.environ.get("NIOME_USE_STRATEGY_AUTO", "1").strip().lower() in (
        "1",
        "true",
        "yes",
    )


def _strategy_preset() -> Dict[str, Any]:
    name = os.environ.get("NIOME_STRATEGY_SENSITIVITY", "sensitive").strip().lower()
    base = strategy_auto.SENSITIVITY_PRESETS.get(
        name, strategy_auto.SENSITIVITY_PRESETS["sensitive"]
    )
    return dict(base)


def _strategy_submit_top_n() -> int:
    return max(1, _env_int("NIOME_STRATEGY_SUBMIT_TOP_N", 30))


def _strategy_submit_allowlist() -> bool:
    return os.environ.get("NIOME_STRATEGY_SUBMIT_ALLOWLIST", "1").strip().lower() in (
        "1",
        "true",
        "yes",
    )


def _load_strategy_selected_key_set(work_dir: Path) -> Set[Tuple[str, int, str, str]]:
    path = work_dir / "strategy_auto" / "rank" / "selected.keys.tsv"
    if not path.exists():
        return set()
    keys: Set[Tuple[str, int, str, str]] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        keys.add(
            (chrom_core(parts[0]), int(parts[1]), parts[2], parts[3])
        )
    return keys


def _load_strategy_rank_table(
    work_dir: Path,
) -> Tuple[List[Tuple[str, int, str, str]], Dict[Tuple[str, int, str, str], str]]:
    """Rank-ordered keys and GT from strategy_auto ranked_candidates.tsv."""
    path = work_dir / "strategy_auto" / "rank" / "ranked_candidates.tsv"
    if not path.exists():
        return [], {}
    ordered: List[Tuple[str, int, str, str]] = []
    gt_map: Dict[Tuple[str, int, str, str], str] = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (
                chrom_core(row["CHROM"]),
                int(row["POS"]),
                row["REF"],
                row["ALT"],
            )
            ordered.append(key)
            gt = (row.get("GT") or "").strip()
            if gt and gt not in (".", "./."):
                gt_map[key] = _normalize_submission_gt(gt) or gt
    return ordered, gt_map


def _strategy_use_rank_gt() -> bool:
    return os.environ.get("NIOME_STRATEGY_USE_RANK_GT", "1").strip().lower() in (
        "1",
        "true",
        "yes",
    )


def _strategy_hybrid_hom_gt() -> bool:
    return os.environ.get("NIOME_STRATEGY_HYBRID_HOM_GT", "1").strip().lower() in (
        "1",
        "true",
        "yes",
    )


def _apply_strategy_rank_genotypes(
    variants: List[VariantRecord],
    work_dir: Path,
) -> None:
    if not _strategy_use_rank_gt():
        return
    _, gt_map = _load_strategy_rank_table(work_dir)
    for variant in variants:
        rank_gt = gt_map.get(_variant_key(variant))
        if rank_gt:
            variant.gt = rank_gt
        if variant.alt_depth >= 1:
            safe_gt = _normalize_submission_gt(variant.gt or "")
            if safe_gt is None:
                variant.gt = "0/1"
        if _strategy_hybrid_hom_gt():
            is_snp = _is_snp(variant.ref, variant.alt)
            af = float(variant.af or 0.0)
            ad = int(variant.alt_depth or 0)
            if is_snp and ad >= 4 and af >= 0.85:
                variant.gt = "1/1"
            elif not is_snp and ad >= 3 and af >= 0.75:
                variant.gt = "1/1"


EVIDENCE_SNV_MIN_AF = 0.18
EVIDENCE_INDEL_MIN_ALT_DEPTH = 3
EVIDENCE_INDEL_MIN_AF = 0.20
EVIDENCE_LOW_DP_HET_MAX_DP = 5
EVIDENCE_LOW_DP_HET_MIN_AF = 0.35
GT_HOM_AF = 0.80
GT_HET_MIN_AF = 0.20


class BamIndexError(RuntimeError):
    """Raised when a BAM exists but a usable index cannot be created or found."""


def _miner_instance_id() -> str:
    """Per-process work dir suffix (MINER_INSTANCE env or pid)."""
    return os.environ.get("MINER_INSTANCE") or str(os.getpid())


def _task_work_dir(base_dir: Path, task_id: str) -> Path:
    return base_dir / "work" / task_id / _miner_instance_id()


def _bam_index_candidate_paths(bam_path: Path) -> List[Path]:
    """Return possible index paths for a BAM (sample.bam.bai and sample.bai)."""
    bam_path = Path(bam_path).resolve()
    return [Path(f"{bam_path}.bai"), bam_path.with_suffix(".bai")]


def _find_bam_index_path(bam_path: Path) -> Optional[Path]:
    for index_path in _bam_index_candidate_paths(bam_path):
        if index_path.exists() and index_path.stat().st_size > 0:
            return index_path.resolve()
    return None


def _has_valid_bam_index(bam_path: Path) -> bool:
    return _find_bam_index_path(bam_path) is not None


def _ensure_bam_index(
    bam_path: Path,
    threads: int = 2,
    logger: Optional[Callable[[str], None]] = None,
) -> Path:
    """Ensure a non-empty BAM index exists beside the resolved BAM path."""
    log = logger or (lambda _msg: None)
    bam_path = Path(bam_path).resolve()
    if not bam_path.exists():
        raise BamIndexError(f"BAM not found: {bam_path}")

    existing = _find_bam_index_path(bam_path)
    if existing is not None:
        log(f"BAM index already present for {bam_path}: {existing}")
        return existing

    lock_path = bam_path.parent / f"{bam_path.name}.index.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        existing = _find_bam_index_path(bam_path)
        if existing is not None:
            log(f"BAM index present after lock for {bam_path}: {existing}")
            return existing

        thread_count = max(1, int(threads))
        cmd = ["samtools", "index", f"-@{thread_count}", str(bam_path)]
        log(f"Creating BAM index: {' '.join(cmd)}")
        index_proc = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
        )
        if index_proc.returncode != 0:
            detail = (index_proc.stderr or index_proc.stdout or "").strip()
            raise BamIndexError(
                f"samtools index failed (exit {index_proc.returncode}) "
                f"for {bam_path}: {detail}"
            )

        existing = _find_bam_index_path(bam_path)
        if existing is None:
            raise BamIndexError(
                f"samtools index completed but no index found for {bam_path}"
            )
        log(f"Created BAM index for {bam_path}: {existing}")
        return existing


def _ensure_aligned_bam_alias(
    work_dir: Path,
    dedup_bam: Path,
    threads: int,
    logger: Optional[Callable[[str], None]] = None,
) -> Path:
    """Keep aligned.bam (+ .bai) as a compatibility alias to the indexed dedup BAM."""
    log = logger or (lambda _msg: None)
    dedup_bam = Path(dedup_bam).resolve()
    _ensure_bam_index(dedup_bam, threads=threads, logger=log)
    dedup_index = _find_bam_index_path(dedup_bam)
    aligned_bam = work_dir / "aligned.bam"

    if not aligned_bam.exists():
        try:
            aligned_bam.symlink_to(dedup_bam)
            log(f"Created aligned.bam symlink -> {dedup_bam}")
        except OSError:
            shutil.copy2(dedup_bam, aligned_bam)
            log(f"Copied dedup BAM to aligned.bam: {aligned_bam}")
            _ensure_bam_index(aligned_bam, threads=threads, logger=log)

    if dedup_index is not None:
        aligned_index = Path(f"{aligned_bam}.bai")
        if not aligned_index.exists():
            try:
                aligned_index.symlink_to(dedup_index)
                log(f"Created aligned.bam.bai symlink -> {dedup_index}")
            except OSError:
                shutil.copy2(dedup_index, aligned_index)
                log(f"Copied BAM index to {aligned_index}")

    return aligned_bam


def build_empty_miner_response(base_dir: Path) -> Dict[str, Any]:
    """Schema-valid empty submission for miner forward error handling."""
    config = _config_from_env(base_dir)
    region_chrom = _region_chrom(DEFAULT_REGION) or "chr7"
    contig_length = _reference_contig_length(config.reference_fasta, region_chrom)
    vcf_content = _build_submission_vcf(
        [],
        config.reference_header,
        region_chrom,
        contig_length,
        config.enable_af_esp,
    )
    return {
        "vcf_content": vcf_content,
        "cftr_annotations": {},
        "elapsed_time": 0.0,
    }


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
    work_dir = _task_work_dir(config.base_dir, safe_task_id)
    output_dir = config.base_dir / "outputs" / safe_task_id / _miner_instance_id()
    for directory in (task_dir, reads_dir, work_dir, output_dir):
        directory.mkdir(parents=True, exist_ok=True)
    log(
        f"CFTR miner task {task_id}: work dir {work_dir} "
        f"output dir {output_dir} (instance={_miner_instance_id()})"
    )

    task_json_path = task_dir / "task.json"
    task_json_path.write_text(json.dumps(task_data, indent=2, sort_keys=True), encoding="utf-8")

    read1_path = reads_dir / "reads_1.fq"
    read2_path = reads_dir / "reads_2.fq"
    log(f"CFTR miner task {task_id}: downloading FASTQ inputs")
    _download_file(read1_url, read1_path)
    _download_file(read2_url, read2_path)

    for tool in ("bwa", "samtools", "bcftools", "tabix"):
        _ensure_executable(tool)
    if not config.gatk_path.exists():
        raise FileNotFoundError(f"GATK executable not found: {config.gatk_path}")
    if not config.reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {config.reference_fasta}")

    expected_variant_count = int(task_data.get("expected_variant_count") or 0)
    region_chrom = _region_chrom(region) or "chr7"
    contig_length = _reference_contig_length(config.reference_fasta, region_chrom)

    log(f"CFTR miner task {task_id}: preparing reference indexes")
    _ensure_reference_indexes(config, log)
    if _use_strategy_auto():
        log(
            f"CFTR miner task {task_id}: strategy_auto=on "
            f"sensitivity={os.environ.get('NIOME_STRATEGY_SENSITIVITY', 'sensitive')} "
            f"submit_top_n={_strategy_submit_top_n()} "
            f"allowlist={int(_strategy_submit_allowlist())} "
            f"rank_gt={int(_strategy_use_rank_gt())}"
        )

    log(
        f"CFTR miner task {task_id}: aligning reads (bwa mem) and "
        f"preprocessing BAM (GATK RG + MarkDuplicates)"
    )
    dedup_bam = _prepare_dedup_bam_from_reads(
        config, read1_path, read2_path, work_dir, safe_task_id, log
    )
    evidence_bam = Path(dedup_bam).resolve()
    _ensure_bam_index(evidence_bam, threads=config.threads, logger=log)
    log(f"CFTR miner task {task_id}: evidence BAM {evidence_bam}")
    index_path = _find_bam_index_path(evidence_bam)
    log(
        f"CFTR miner task {task_id}: evidence index "
        f"{'present' if index_path else 'missing'} "
        f"({index_path if index_path else 'none'})"
    )
    bam_path = _ensure_aligned_bam_alias(
        work_dir, evidence_bam, threads=config.threads, logger=log
    )

    clinvar_panel = _load_merged_clinvar_panel(config.clinvar_panel)
    drug_panel = _load_drug_panel(config.drug_panel)

    pop_lookup: Optional[PopulationAfLookup] = None
    if config.enable_af_esp:
        try:
            pop_lookup = _load_population_af_lookup(config)
            log(
                f"CFTR miner task {task_id}: gnomAD cache for INFO/tie-break only "
                f"({len(pop_lookup.by_allele)} alleles)"
            )
        except Exception as exc:
            log(f"CFTR miner task {task_id}: gnomAD cache skipped: {exc}")

    log(f"CFTR miner task {task_id}: collecting GATK + panel candidates")
    raw_candidates, caller_paths = _collect_all_candidates(
        config,
        evidence_bam,
        region,
        work_dir,
        clinvar_panel,
        log,
    )
    log(f"CFTR miner task {task_id}: union {len(raw_candidates)} raw candidate allele(s)")

    log(f"CFTR miner task {task_id}: normalizing candidate alleles (bcftools norm)")
    normalized_candidates = _normalize_candidate_alleles(
        raw_candidates,
        config,
        work_dir,
        region_chrom,
        contig_length,
        log,
    )
    log(
        f"CFTR miner task {task_id}: {len(normalized_candidates)} candidate(s) after norm"
    )

    log(f"CFTR miner task {task_id}: computing read evidence from BAM")
    evidence_rows = _apply_evidence_engine(
        normalized_candidates,
        evidence_bam,
        config,
        region,
        config.reference_fasta,
        clinvar_panel,
        log,
        work_dir=work_dir,
    )
    accepted = [row for row in evidence_rows if row.evidence_status == "ACCEPTED"]
    rejected = len(evidence_rows) - len(accepted)
    log(
        f"CFTR miner task {task_id}: evidence filter {len(accepted)} accepted, "
        f"{rejected} rejected"
    )
    fp_before = len(accepted)
    accepted = _filter_nonpanel_gatk_false_positives(accepted, clinvar_panel)
    if len(accepted) < fp_before:
        log(
            f"CFTR miner task {task_id}: dropped {fp_before - len(accepted)} "
            f"non-panel GATK FP(s) after evidence"
        )

    log(
        f"CFTR miner task {task_id}: ranking/selecting "
        f"(expected_variant_count={expected_variant_count})"
    )
    if _use_strategy_auto():
        selected_variants = _select_strategy_aligned_variants(
            evidence_rows,
            work_dir,
            expected_variant_count,
            pop_lookup,
            log,
        )
    else:
        selected_variants = _select_evidence_ranked_variants(
            accepted,
            expected_variant_count,
            pop_lookup,
            log,
        )

    evidence_tsv = output_dir / "candidate_evidence.tsv"
    _write_candidate_evidence_tsv(evidence_tsv, evidence_rows)

    vcf_content = _prepare_submission_vcf(
        selected_variants,
        config,
        work_dir,
        region_chrom,
        contig_length,
        log,
        pop_lookup=pop_lookup,
        clinvar_panel=clinvar_panel,
        bam_path=evidence_bam,
        region=region,
    )
    try:
        verify_vcf_for_validator_scoring(
            vcf_content,
            config.reference_fasta,
            work_dir=work_dir / "validator_vcf_check",
        )
    except Exception as exc:
        log(f"WARNING: synapse VCF failed validator norm pre-check: {exc}")

    vcf_path = output_dir / "submission.vcf"
    synapse_vcf_path = output_dir / "synapse.vcf"
    output_lock = output_dir / ".submission.write.lock"
    with open(output_lock, "w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        # Synapse wire format: exact validator-safe VCF from pipeline (no extra headers).
        synapse_vcf_path.write_text(vcf_content, encoding="utf-8")
        vcf_lines = vcf_content.splitlines()
        if vcf_lines:
            vcf_lines.insert(1, f"##niome_task_id={task_id}")
            vcf_lines.insert(2, f"##niome_instance={_miner_instance_id()}")
        vcf_path.write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")

    cftr_annotations = _build_cftr_annotations(
        selected_variants,
        clinvar_panel,
        drug_panel,
    )
    annotation_path = output_dir / "annotation.json"
    with open(output_lock, "w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        annotation_path.write_text(
            json.dumps(cftr_annotations, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    elapsed_time = time.time() - started_at
    # Synapse must return exactly what was written to disk (not a stale shared buffer).
    vcf_content = synapse_vcf_path.read_text(encoding="utf-8")
    variant_lines = [
        line for line in vcf_content.splitlines() if line and not line.startswith("#")
    ]
    log(
        f"CFTR miner task {task_id}: submitted {len(selected_variants)} variants "
        f"({len(variant_lines)} VCF rows) with {len(cftr_annotations)} annotations "
        f"in {elapsed_time:.2f}s -> {synapse_vcf_path}"
    )

    return {
        "vcf_content": vcf_content,
        "cftr_annotations": cftr_annotations,
        "elapsed_time": elapsed_time,
        "counts": {
            "raw_candidates": len(raw_candidates),
            "normalized_candidates": len(normalized_candidates),
            "evidence_accepted": len(accepted),
            "evidence_rejected": rejected,
            "submitted": len(selected_variants),
            "annotated": len(cftr_annotations),
            "expected_variant_count": expected_variant_count,
        },
        "paths": {
            "task_json": str(task_json_path),
            "read1_fastq": str(read1_path),
            "read2_fastq": str(read2_path),
            "aligned_bam": str(bam_path),
            "dedup_bam": str(evidence_bam),
            "evidence_bam": str(evidence_bam),
            **caller_paths,
            "candidate_evidence_tsv": str(evidence_tsv),
            "submission_vcf": str(vcf_path),
            "synapse_vcf": str(synapse_vcf_path),
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
    gatk_path: Path
    caller: str
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
    panel_rescue_max: int = 10
    panel_rescue_min_dp: int = 10


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


def _is_actionable_panel_significance(clinical_significance: str) -> bool:
    """Panel sites worth probing when read-backed (pathogenic / VUS / likely)."""
    sig = clinical_significance.lower().replace("_", " ")
    if "benign" in sig and "pathogenic" not in sig:
        return False
    return any(
        token in sig
        for token in ("pathogenic", "uncertain", "likely pathogenic", "vus")
    )


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


def _gt_from_gnomad(variant: VariantRecord) -> Optional[str]:
    """Infer diploid GT from gnomAD population AF (AF_ESP) and sample allele fraction."""
    pop_af = variant.af_esp
    if pop_af is None:
        return None
    af = variant.af
    if variant.alt_depth < 1:
        return None
    if pop_af >= 0.90:
        return "1/1" if af >= 0.35 else None
    if pop_af >= 0.25:
        if af >= 0.72:
            return "1/1"
        if af >= 0.15:
            return "0/1"
    if pop_af < 0.01:
        if af >= 0.82:
            return "1/1"
        if af >= 0.10:
            return "0/1"
    if 0.20 <= af <= 0.80:
        return "0/1"
    return None


def _finalize_genotype(
    variant: VariantRecord,
    clinical_significance: str = "",
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> str:
    """Finalize GT from BAM read evidence only (gnomAD does not set GT)."""
    _ = clinical_significance
    _ = pop_lookup
    is_snp = _is_snp(variant.ref, variant.alt)
    gt, _ = _gt_from_read_evidence_only(
        variant.ref_depth,
        variant.alt_depth,
        variant.dp,
        variant.af,
        is_snp,
        not is_snp,
        ref=variant.ref,
        alt=variant.alt,
        clinical_significance=clinical_significance,
        is_panel=variant.is_panel,
    )
    if gt:
        return gt
    if variant.gt and variant.gt not in ("0/0", "./.", "."):
        return _normalize_submission_gt(variant.gt) or "0/1"
    return "0/1"


def _is_primary_caller_source(source: str) -> bool:
    """True for GATK HaplotypeCaller or legacy bcftools mpileup/call variants."""
    return "gatk" in source or "standard" in source


def _use_gatk_caller(config: CftrMinerConfig) -> bool:
    return config.caller.strip().lower() != "bcftools"


def _nonpanel_standard_drop_reason(variant: VariantRecord) -> Optional[str]:
    """Drop high-confidence false-positive caller-only calls (not in ClinVar panel)."""
    if variant.is_panel or not _is_primary_caller_source(variant.source):
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
        if _is_primary_caller_source(record.source):
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


def _compute_detailed_pileup_evidence(
    bam: pysam.AlignmentFile,
    contig: str,
    pos: int,
    ref: str,
    alt: str,
    config: CftrMinerConfig,
) -> Tuple[int, int, int, float, bool, float, float]:
    """Return ref/alt depths, AF, strand_balanced, mean_alt_bq, mean_mapq."""
    if _is_snp(ref, alt):
        return _compute_snp_pileup_metrics(bam, contig, pos, ref, alt, config)
    ref_d, alt_d, dp, af = count_allele_support_with_pysam(
        bam, contig, pos, ref, alt
    )
    strand_balanced = True
    mean_bq = float(config.min_baseq)
    mean_mq = float(config.min_mapq)
    if alt_d >= 2:
        strand_balanced = True
    return ref_d, alt_d, dp, af, strand_balanced, mean_bq, mean_mq


def _compute_snp_pileup_metrics(
    bam: pysam.AlignmentFile,
    contig: str,
    pos: int,
    ref: str,
    alt: str,
    config: CftrMinerConfig,
) -> Tuple[int, int, int, float, bool, float, float]:
    ref_upper = ref.upper()
    alt_upper = alt.upper()
    ref_depth = 0
    alt_depth = 0
    alt_forward = 0
    alt_reverse = 0
    bq_sum = 0.0
    bq_count = 0
    mq_sum = 0.0
    mq_count = 0
    min_bq = config.min_baseq
    min_mq = config.min_mapq

    for column in bam.pileup(
        contig,
        pos - 1,
        pos,
        truncate=True,
        stepper="all",
        min_base_quality=min_bq,
        min_mapping_quality=min_mq,
    ):
        if column.reference_pos != pos - 1:
            continue
        for pileupread in column.pileups:
            if pileupread.is_refskip or pileupread.is_del:
                continue
            query_position = pileupread.query_position
            if query_position is None:
                continue
            aln = pileupread.alignment
            query_sequence = aln.query_sequence
            if query_sequence is None or query_position >= len(query_sequence):
                continue
            base = query_sequence[query_position].upper()
            if aln.mapping_quality < min_mq:
                continue
            if base == ref_upper:
                ref_depth += 1
            elif base == alt_upper:
                alt_depth += 1
                if aln.is_reverse:
                    alt_reverse += 1
                else:
                    alt_forward += 1
                if aln.query_qualities and query_position < len(aln.query_qualities):
                    bq_sum += aln.query_qualities[query_position]
                    bq_count += 1
                mq_sum += aln.mapping_quality
                mq_count += 1

    dp = ref_depth + alt_depth
    af = (alt_depth / dp) if dp > 0 else 0.0
    strand_balanced = True
    if alt_depth >= 2:
        minor = min(alt_forward, alt_reverse)
        major = max(alt_forward, alt_reverse)
        strand_balanced = minor >= 1 and (minor / major) >= 0.25
    mean_bq = (bq_sum / bq_count) if bq_count else float(min_bq)
    mean_mq = (mq_sum / mq_count) if mq_count else float(min_mq)
    return ref_depth, alt_depth, dp, af, strand_balanced, mean_bq, mean_mq


def _variant_in_region(variant: VariantRecord, region: str) -> bool:
    chrom, start, end = _parse_region_bounds(region)
    if chrom_core(variant.chrom) != chrom_core(chrom):
        return False
    return start <= variant.pos <= end


def _check_ref_match(
    variant: VariantRecord,
    fasta: pysam.FastaFile,
) -> bool:
    contig = _resolve_reference_contig(fasta, variant.chrom)
    if contig is None:
        return False
    ref_seq = _reference_bases_at(fasta, contig, variant.pos, len(variant.ref))
    return ref_seq is not None and ref_seq == variant.ref


def _caller_support_count(variant: VariantRecord) -> int:
    return sum(
        (
            bool(variant.source_gatk),
            bool(variant.source_bcftools),
            bool(variant.source_assembly),
            bool(variant.source_panel_probe),
        )
    )


def _primary_caller_count(variant: VariantRecord) -> int:
    """GATK + bcftools only (used for SNV AF relax, not assembly/panel)."""
    return int(bool(variant.source_gatk)) + int(bool(variant.source_bcftools))


def _is_relaxed_indel_evidence(variant: VariantRecord, is_indel: bool) -> bool:
    if not is_indel:
        return False
    if variant.source == "mpileup_indel_discovery":
        return True
    return bool(
        variant.source_assembly
        and not variant.source_gatk
        and not variant.source_bcftools
        and variant.alt_depth >= 2
        and len(variant.alt) - len(variant.ref) <= 3
    )


def _gt_from_read_evidence_only(
    ref_depth: int,
    alt_depth: int,
    dp: int,
    af: float,
    is_snp: bool,
    is_indel: bool,
    caller_count: int = 0,
    primary_caller_count: int = 0,
    panel_backed: bool = False,
    relaxed_indel: bool = False,
    ref: str = "",
    alt: str = "",
    clinical_significance: str = "",
    is_panel: bool = False,
) -> Tuple[Optional[str], Optional[str]]:
    """Evidence filter + GT from BAM AD/AF (GT uses validator-aligned infer_gt_from_depth)."""
    min_alt_depth = EVIDENCE_MIN_ALT_DEPTH
    if _use_strategy_auto():
        min_alt_depth = max(1, _env_int("NIOME_EVIDENCE_MIN_ALT_DEPTH", 1))
    if alt_depth < min_alt_depth:
        return None, "reject_alt_depth_lt_2"
    if dp < EVIDENCE_MIN_DP:
        return None, "reject_dp_lt_3"
    min_af = EVIDENCE_MIN_AF
    if _use_strategy_auto():
        min_af = min(min_af, _env_float("NIOME_EVIDENCE_MIN_AF", 0.08))
    if relaxed_indel and is_indel and alt_depth >= 2:
        min_af = min(min_af, 0.10)
    if panel_backed and is_indel and alt_depth >= 2:
        min_af = min(min_af, 0.10)
    if _use_strategy_auto() and panel_backed and alt_depth >= min_alt_depth:
        min_af = min(min_af, 0.05)
    if af < min_af:
        return None, "reject_af_lt_min"

    snv_min_af = EVIDENCE_SNV_MIN_AF
    if _use_strategy_auto():
        snv_min_af = min(snv_min_af, _env_float("NIOME_EVIDENCE_SNV_MIN_AF", 0.10))
    if primary_caller_count >= 2 and alt_depth >= 3:
        snv_min_af = min(snv_min_af, 0.15)
    if primary_caller_count >= 2 and alt_depth >= 4:
        snv_min_af = min(snv_min_af, 0.12)
    if panel_backed and primary_caller_count >= 1 and alt_depth >= 3:
        snv_min_af = min(snv_min_af, 0.15)

    if (
        is_snp
        and af < snv_min_af
        and dp > EVIDENCE_LOW_DP_HET_MAX_DP
        and not (_use_strategy_auto() and panel_backed and alt_depth >= min_alt_depth)
    ):
        return None, "reject_snv_af_lt_0.18"
    if is_indel:
        if relaxed_indel and alt_depth >= 2 and af >= min_af:
            pass
        else:
            multi_caller_indel = (
                caller_count >= 2 and alt_depth >= 2 and af >= min_af
            )
            if not multi_caller_indel:
                if alt_depth < EVIDENCE_INDEL_MIN_ALT_DEPTH and af < EVIDENCE_INDEL_MIN_AF:
                    return None, "reject_indel_weak_support"
                if af < EVIDENCE_INDEL_MIN_AF and dp > EVIDENCE_LOW_DP_HET_MAX_DP:
                    return None, "reject_indel_af_lt_0.20"

    if dp <= EVIDENCE_LOW_DP_HET_MAX_DP:
        if alt_depth >= EVIDENCE_MIN_ALT_DEPTH and af >= EVIDENCE_LOW_DP_HET_MIN_AF:
            return (
                infer_gt_from_depth(
                    ref_depth,
                    alt_depth,
                    clinical_significance,
                    ref=ref,
                    alt=alt,
                    is_panel=is_panel or panel_backed,
                ),
                None,
            )
        return None, "reject_low_dp_insufficient_af"

    het_min_af = GT_HET_MIN_AF
    if is_snp and alt_depth >= 7:
        het_min_af = min(het_min_af, 0.18)
    if is_snp and caller_count >= 2 and alt_depth >= 4:
        het_min_af = min(het_min_af, 0.15)
    if is_snp and primary_caller_count >= 2 and alt_depth >= 3:
        het_min_af = min(het_min_af, 0.15)
    if is_indel and primary_caller_count >= 2 and alt_depth >= 2:
        het_min_af = min(het_min_af, 0.15)
    if panel_backed and alt_depth >= 3:
        het_min_af = min(het_min_af, 0.15)
    if panel_backed and is_indel and alt_depth >= 2:
        het_min_af = min(het_min_af, 0.10)
    if relaxed_indel and alt_depth >= 2:
        het_min_af = min(het_min_af, 0.10)

    if af < het_min_af:
        return None, "reject_af_below_het_threshold"
    gt = infer_gt_from_depth(
        ref_depth,
        alt_depth,
        clinical_significance,
        ref=ref,
        alt=alt,
        is_panel=is_panel or panel_backed,
    )
    # Top leaderboard miners keep more hets under allele imbalance; avoid auto-hom upgrades.
    if not _use_strategy_auto():
        if panel_backed and is_snp and gt == "0/1":
            sig = clinical_significance.lower().replace("_", " ")
            if "pathogenic" in sig and "benign" not in sig:
                if af >= 0.58 and alt_depth >= 5:
                    gt = "1/1"
                elif af >= 0.52 and alt_depth >= 8:
                    gt = "1/1"
    if panel_backed and is_indel and gt == "0/1":
        sig = clinical_significance.lower().replace("_", " ")
        if "pathogenic" in sig and af >= 0.50 and alt_depth >= 4:
            gt = "1/1"
    return gt, None


def _hard_reject_evidence(
    variant: VariantRecord,
    region: str,
    ref_match: bool,
) -> Optional[str]:
    if not _variant_in_region(variant, region):
        return "reject_outside_region"
    if not ref_match:
        return "reject_ref_mismatch"
    if not _is_submission_sized_allele(variant.ref) or not _is_submission_sized_allele(
        variant.alt
    ):
        return "reject_invalid_allele"
    if variant.ref_depth < 2 and not variant.source_panel_rescue:
        return "reject_insufficient_ref_depth"
    if variant.alt_depth == 0 and not variant.source_panel_rescue:
        return "reject_alt_depth_zero"
    if variant.af <= 0.0 and not variant.source_panel_rescue:
        return "reject_af_zero"
    if variant.source_panel_probe and variant.alt_depth < EVIDENCE_MIN_ALT_DEPTH:
        return "reject_panel_no_alt_reads"
    return None


def _compute_evidence_confidence_score(variant: VariantRecord) -> float:
    score = 0.0
    if variant.source_panel_rescue:
        score += 60.0
        if variant.dp >= 15:
            score += 10.0
    score += min(40.0, variant.alt_depth * 5)
    score += min(25.0, variant.af * 25.0)
    if variant.source_gatk:
        score += 10.0
    if variant.source_bcftools:
        score += 10.0
    if variant.source_assembly:
        score += 15.0
    if (
        variant.source_assembly
        and not _is_snp(variant.ref, variant.alt)
        and variant.alt_depth >= 2
        and variant.af >= 0.10
    ):
        score += 20.0
    if variant.source_panel_probe or variant.is_panel:
        score += 15.0
    score += max(0, _caller_support_count(variant) - 1) * 8.0
    if variant.strand_balanced:
        score += 5.0
    if variant.mean_mapq >= 20:
        score += 5.0
    if variant.mean_alt_bq >= 20:
        score += 5.0
    return score


def _merge_candidate_into_pool(
    pool: Dict[Tuple[str, int, str, str], VariantRecord],
    incoming: VariantRecord,
) -> None:
    key = _variant_key(incoming)
    existing = pool.get(key)
    if existing is None:
        pool[key] = incoming
        return
    existing.source_gatk = existing.source_gatk or incoming.source_gatk
    existing.source_bcftools = existing.source_bcftools or incoming.source_bcftools
    existing.source_assembly = existing.source_assembly or incoming.source_assembly
    existing.source_panel_probe = (
        existing.source_panel_probe
        or incoming.source_panel_probe
        or incoming.source in ("panel_pileup", "panel_probe")
    )
    existing.source_panel_rescue = (
        existing.source_panel_rescue or incoming.source_panel_rescue
    )
    existing.is_panel = existing.is_panel or incoming.is_panel
    if incoming.variation_id and not existing.variation_id:
        existing.variation_id = incoming.variation_id
    if incoming.qual is not None and (
        existing.qual is None or incoming.qual > existing.qual
    ):
        existing.qual = incoming.qual
    if incoming.alt_depth > existing.alt_depth:
        existing.ref_depth = incoming.ref_depth
        existing.alt_depth = incoming.alt_depth
        existing.dp = incoming.dp
        existing.af = incoming.af
        existing.gt = incoming.gt
    if incoming.evidence_score > existing.evidence_score:
        existing.evidence_score = incoming.evidence_score
        existing.evidence_status = incoming.evidence_status
        existing.evidence_reason = incoming.evidence_reason


def _dedupe_variants_by_key(variants: List[VariantRecord]) -> List[VariantRecord]:
    pool: Dict[Tuple[str, int, str, str], VariantRecord] = {}
    for variant in variants:
        _merge_candidate_into_pool(pool, variant)
    return sorted(
        pool.values(),
        key=lambda item: (_chrom_sort_key(item.chrom), item.pos, item.ref, item.alt),
    )


def _variants_from_vcf_path(vcf_path: Path, source_label: str) -> List[VariantRecord]:
    records: List[VariantRecord] = []
    for variant in _parse_vcf_records(vcf_path, source_label):
        if not _is_valid_allele(variant.ref) or not _is_valid_allele(variant.alt):
            continue
        if source_label == "gatk":
            variant.source_gatk = True
        elif source_label in ("bcftools", "standard"):
            variant.source_bcftools = True
        elif source_label == "assembly":
            variant.source_assembly = True
        records.append(variant)
    return records


def _pathogenic_panel_entries_in_region(
    config: CftrMinerConfig,
    region: str,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
) -> List[PanelVariant]:
    """Pathogenic/VUS panel alleles in the task region (supplemental union)."""
    entries = _load_panel_variants_in_region(
        config.clinvar_panel, region, config.base_dir
    )
    filtered: List[PanelVariant] = []
    for entry in entries:
        key = (chrom_core(entry.chrom), entry.pos, entry.ref, entry.alt)
        panel_meta = clinvar_panel.get(key, {})
        clin_sig = panel_meta.get("clinical_significance", "") or entry.clinical_significance
        if not _is_actionable_panel_significance(clin_sig):
            continue
        if panel_meta.get("variation_id") and not entry.variation_id:
            entry = PanelVariant(
                variation_id=str(panel_meta["variation_id"]),
                chrom=entry.chrom,
                pos=entry.pos,
                ref=entry.ref,
                alt=entry.alt,
                hgvs=panel_meta.get("hgvs", entry.hgvs),
                clinical_significance=clin_sig,
            )
        filtered.append(entry)
    return filtered


def _variant_from_strategy_selection(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    agg: Dict[str, Any],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
) -> VariantRecord:
    """Map strategy_auto ranked entry into a miner VariantRecord."""
    dp = max(int(agg.get("max_dp") or 0), 1)
    alt_d = max(int(agg.get("max_alt") or 0), 0)
    ref_d = max(int(agg.get("max_ref") or 0), dp - alt_d)
    af = alt_d / dp if dp else 0.0
    families = agg.get("families") or set()

    variant = VariantRecord(
        chrom=_vcf_chrom(chrom),
        pos=pos,
        ref=ref,
        alt=alt,
        qual=float(agg.get("max_qual") or 0.0) or None,
        filter_value="PASS",
        gt=str(agg.get("best_gt") or "0/1"),
        dp=dp,
        ref_depth=ref_d,
        alt_depth=alt_d,
        af=af,
        ad=f"{ref_d},{alt_d}",
        source="strategy_auto",
        read_backed=alt_d >= 1,
    )
    if "deepvariant" in families or "gatk" in families:
        variant.source_gatk = True
    if "bcftools" in families:
        variant.source_bcftools = True
    if "freebayes" in families or "weakscan" in families:
        variant.source_assembly = True

    panel_key = (chrom_core(chrom), pos, ref, alt)
    panel_meta = clinvar_panel.get(panel_key)
    if panel_meta or agg.get("clinvar"):
        variant.is_panel = bool(panel_meta)
        variant.source_panel_probe = bool(panel_meta)
        if panel_meta:
            variant.variation_id = str(panel_meta.get("variation_id", ""))
    return variant


def _merge_strategy_auto_candidates(
    config: CftrMinerConfig,
    bam_path: Path,
    region: str,
    work_dir: Path,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    pool: Dict[Tuple[str, int, str, str], VariantRecord],
    paths: Dict[str, str],
    gatk_plain: Path,
    logger: Callable[[str], None],
) -> None:
    """Run sensitive multi-caller strategy and merge selected alleles into the pool."""
    log = logger
    preset = _strategy_preset()
    run_dv = os.environ.get("NIOME_STRATEGY_DEEPVARIANT", "0").strip().lower() in (
        "1",
        "true",
        "yes",
    )
    log(
        f"Strategy auto: preset={os.environ.get('NIOME_STRATEGY_SENSITIVITY', 'sensitive')} "
        f"deepvariant={run_dv}"
    )
    _ranked, selected, strat_paths = strategy_auto.run_for_miner_bam(
        ref=str(config.reference_fasta),
        bam=str(bam_path),
        region=region,
        work_dir=str(work_dir),
        clinvar_panel=clinvar_panel,
        preset=preset,
        threads=config.threads,
        logger=log,
        gatk_plain_vcf=str(gatk_plain) if gatk_plain.exists() else None,
        enable_deepvariant=run_dv,
        base_dir=str(config.base_dir),
    )
    paths.update({f"strategy_{k}": v for k, v in strat_paths.items()})
    merged = 0
    for _score, key, agg in selected:
        chrom, pos, ref, alt = key
        variant = _variant_from_strategy_selection(
            chrom, pos, ref, alt, agg, clinvar_panel
        )
        _merge_candidate_into_pool(pool, variant)
        merged += 1
    log(f"Strategy auto merged {merged} selected allele(s) into candidate pool")


def _collect_all_candidates(
    config: CftrMinerConfig,
    bam_path: Path,
    region: str,
    work_dir: Path,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    logger: Callable[[str], None],
) -> Tuple[List[VariantRecord], Dict[str, str]]:
    """Collect candidates from GATK plus read-backed pathogenic panel probes."""
    log = logger
    paths: Dict[str, str] = {}
    pool: Dict[Tuple[str, int, str, str], VariantRecord] = {}

    gatk_raw = work_dir / "calls.gatk.raw.vcf.gz"
    gatk_plain = work_dir / "calls.gatk.vcf"
    log("Running GATK HaplotypeCaller")
    _call_gatk_haplotypecaller(config, bam_path, region, gatk_raw)
    _run_command(
        ["bcftools", "view", "-Ov", "-o", str(gatk_plain), str(gatk_raw)],
        "convert GATK VCF to plain text",
    )
    paths["gatk_raw_vcf"] = str(gatk_raw)
    paths["gatk_plain_vcf"] = str(gatk_plain)
    gatk_variants = _variants_from_vcf_path(gatk_plain, "gatk")
    for variant in gatk_variants:
        _merge_candidate_into_pool(pool, variant)
    log(f"GATK contributed {len(gatk_variants)} allele(s)")

    panel_entries = _pathogenic_panel_entries_in_region(
        config, region, clinvar_panel
    )
    if panel_entries:
        log(
            f"Panel pileup scan for {len(panel_entries)} pathogenic/VUS site(s) "
            f"in region"
        )
        pileup_variants = _panel_pileup_scan(bam_path, panel_entries, config)
        pileup_added = 0
        for variant in pileup_variants:
            if variant.alt_depth < EVIDENCE_MIN_ALT_DEPTH or variant.af <= 0.0:
                continue
            key = _variant_key(variant)
            panel_meta = clinvar_panel.get(key, {})
            if panel_meta and not variant.variation_id:
                variant.variation_id = str(panel_meta.get("variation_id", ""))
            variant.source_panel_probe = True
            _merge_candidate_into_pool(pool, variant)
            pileup_added += 1
        log(f"Panel pileup contributed {pileup_added} read-backed allele(s)")

        missing_probe = [
            entry
            for entry in panel_entries
            if (chrom_core(entry.chrom), entry.pos, entry.ref, entry.alt)
            not in pool
        ]
        if missing_probe:
            probed = _probe_panel_candidates_with_reads(
                bam_path,
                missing_probe,
                config,
                clinvar_panel,
                config.reference_fasta,
            )
            probe_added = 0
            for variant in probed:
                variant.source_panel_probe = True
                _merge_candidate_into_pool(pool, variant)
                probe_added += 1
            log(f"Panel probe contributed {probe_added} additional allele(s)")

        if config.panel_rescue_max > 0:
            already_keys = set(pool.keys())
            rescued = _collect_panel_rescue_candidates(
                bam_path,
                panel_entries,
                config,
                clinvar_panel,
                already_keys,
                pileup_variants=pileup_variants,
            )
            rescue_added = 0
            for variant in rescued:
                _merge_candidate_into_pool(pool, variant)
                rescue_added += 1
                log(
                    f"  panel_rescue picked chr7:{variant.pos} {variant.ref}>{variant.alt} "
                    f"DP={variant.dp} (vid={variant.variation_id})"
                )
            log(
                f"Panel zero-alt rescue contributed {rescue_added} allele(s) "
                f"(cap {config.panel_rescue_max}, min_dp {config.panel_rescue_min_dp})"
            )

    if _use_strategy_auto():
        _merge_strategy_auto_candidates(
            config,
            bam_path,
            region,
            work_dir,
            clinvar_panel,
            pool,
            paths,
            gatk_plain,
            log,
        )

    return list(pool.values()), paths


def _filter_nonpanel_gatk_false_positives(
    variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
) -> List[VariantRecord]:
    """Drop GATK-only calls with weak het-band support and no panel match."""
    kept: List[VariantRecord] = []
    for variant in variants:
        key = _variant_key(variant)
        if variant.is_panel or key in clinvar_panel or variant.source_panel_probe:
            kept.append(variant)
            continue
        drop_reason = _nonpanel_standard_drop_reason(variant)
        if drop_reason:
            continue
        if not variant.source_gatk or variant.source_bcftools or variant.source_assembly:
            kept.append(variant)
            continue
        is_snp = _is_snp(variant.ref, variant.alt)
        if (
            is_snp
            and 0.50 <= variant.af <= 0.56
            and variant.alt_depth < 9
            and variant.dp < 24
        ):
            continue
        if (
            is_snp
            and variant.af >= 0.72
            and variant.alt_depth < 6
            and not variant.is_panel
        ):
            continue
        kept.append(variant)
    return kept


def _discover_mpileup_snaps_for_evidence(
    bam_path: Path,
    reference_fasta: Path,
    region: str,
    config: CftrMinerConfig,
) -> List[VariantRecord]:
    """Local assembly-style discovery (mpileup SNPs + indels + pileup confirm)."""
    discovered = _discover_mpileup_snps(bam_path, reference_fasta, region, config)
    for variant in discovered:
        variant.source = "assembly"
        variant.source_assembly = True
        if variant.alt_depth >= EVIDENCE_MIN_ALT_DEPTH and variant.af > 0:
            variant.read_backed = True
    return [v for v in discovered if v.alt_depth >= 1]


def _probe_panel_candidates_with_reads(
    bam_path: Path,
    panel_entries: List[PanelVariant],
    config: CftrMinerConfig,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    reference_fasta: Path,
) -> List[VariantRecord]:
    """Panel is probe-only: only enter pool when ALT reads exist (>=2, AF>0)."""
    probed: List[VariantRecord] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
        str(reference_fasta)
    ) as fasta:
        for entry in panel_entries:
            if not _is_submission_sized_allele(entry.ref) or not _is_submission_sized_allele(
                entry.alt
            ):
                continue
            contig = _resolve_bam_contig(bam, entry.chrom)
            if contig is None:
                continue
            probe_variant = VariantRecord(
                chrom=_vcf_chrom(contig),
                pos=entry.pos,
                ref=entry.ref,
                alt=entry.alt,
                filter_value="PASS",
                gt="0/1",
                source="panel_probe",
            )
            if not _check_ref_match(probe_variant, fasta):
                continue
            (
                ref_d,
                alt_d,
                dp,
                af,
                _strand,
                _bq,
                _mq,
            ) = _compute_detailed_pileup_evidence(
                bam, contig, entry.pos, entry.ref, entry.alt, config
            )
            if ref_d < 2 or alt_d < EVIDENCE_MIN_ALT_DEPTH or af <= 0.0:
                continue
            panel_entry = clinvar_panel.get(
                (chrom_core(entry.chrom), entry.pos, entry.ref, entry.alt), {}
            )
            clin_sig = str(
                panel_entry.get("clinical_significance", entry.clinical_significance)
            )
            gt = infer_gt_from_depth(
                ref_d,
                alt_d,
                clin_sig,
                ref=entry.ref,
                alt=entry.alt,
                is_panel=True,
            )
            probed.append(
                VariantRecord(
                    chrom=_vcf_chrom(contig),
                    pos=entry.pos,
                    ref=entry.ref,
                    alt=entry.alt,
                    filter_value="PASS",
                    gt=gt,
                    dp=dp,
                    ref_depth=ref_d,
                    alt_depth=alt_d,
                    af=af,
                    source="panel_probe",
                    variation_id=entry.variation_id,
                    is_panel=True,
                    read_backed=True,
                    source_panel_probe=True,
                )
            )
            if panel_entry and not probed[-1].variation_id:
                probed[-1].variation_id = str(panel_entry.get("variation_id", ""))
    return probed


_PANEL_RESCUE_SIG_PRIORITY = {
    "pathogenic": 0,
    "likely_pathogenic": 1,
    "uncertain": 2,
    "vus": 2,
    "conflicting": 3,
}


def _panel_rescue_significance_rank(clinical_significance: str) -> int:
    sig = (clinical_significance or "").lower().replace("_", " ")
    if "pathogenic" in sig and "likely" not in sig and "benign" not in sig:
        return 0
    if "likely pathogenic" in sig:
        return 1
    if "uncertain" in sig or "vus" in sig:
        return 2
    if "conflicting" in sig:
        return 3
    return 9


def _collect_panel_rescue_candidates(
    bam_path: Path,
    panel_entries: List[PanelVariant],
    config: CftrMinerConfig,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    already_in_pool: Set[Tuple[str, int, str, str]],
    pileup_variants: Optional[List[VariantRecord]] = None,
) -> List[VariantRecord]:
    """Zero-alt-read panel rescue: include actionable panel variants where the BAM
    is well covered (dp >= panel_rescue_min_dp) but shows no alt reads.

    The Niome simulator deliberately introduces 15-85% haplotype imbalance with low
    coverage; a het variant on the minority haplotype can yield zero alt reads in
    the BAM at a covered position. Standard panel pileup drops these; this function
    rescues them as conservative submission candidates with GT 0/1.

    When `pileup_variants` is supplied (from `_panel_pileup_scan`), this function
    reuses those depths instead of re-pileup-ing the BAM — saves ~25-30s per task.
    """
    if config.panel_rescue_max <= 0 or not panel_entries:
        return []
    rescued: List[VariantRecord] = []
    if pileup_variants is not None:
        sig_by_key = {
            (chrom_core(entry.chrom), entry.pos, entry.ref, entry.alt): entry
            for entry in panel_entries
        }
        for variant in pileup_variants:
            key = (chrom_core(variant.chrom), variant.pos, variant.ref, variant.alt)
            if key in already_in_pool:
                continue
            entry = sig_by_key.get(key)
            if entry is None:
                continue
            if not _is_actionable_panel_significance(entry.clinical_significance):
                continue
            if not _is_submission_sized_allele(entry.ref) or not _is_submission_sized_allele(
                entry.alt
            ):
                continue
            if variant.dp < config.panel_rescue_min_dp:
                continue
            if variant.alt_depth >= EVIDENCE_MIN_ALT_DEPTH:
                continue
            panel_meta = clinvar_panel.get(key, {})
            rescued.append(
                VariantRecord(
                    chrom=variant.chrom,
                    pos=variant.pos,
                    ref=variant.ref,
                    alt=variant.alt,
                    filter_value="PASS",
                    gt="0/1",
                    dp=variant.dp,
                    ref_depth=variant.ref_depth,
                    alt_depth=variant.alt_depth,
                    af=0.0,
                    source="panel_rescue",
                    variation_id=panel_meta.get("variation_id", entry.variation_id),
                    is_panel=True,
                    read_backed=False,
                    source_panel_rescue=True,
                )
            )
    else:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for entry in panel_entries:
                if not _is_actionable_panel_significance(entry.clinical_significance):
                    continue
                if not _is_submission_sized_allele(entry.ref) or not _is_submission_sized_allele(
                    entry.alt
                ):
                    continue
                contig = _resolve_bam_contig(bam, entry.chrom)
                if contig is None:
                    continue
                key = (chrom_core(entry.chrom), entry.pos, entry.ref, entry.alt)
                if key in already_in_pool:
                    continue
                if _is_snp(entry.ref, entry.alt):
                    ref_d, alt_d, dp, _af = count_snp_support_with_pysam(
                        bam, contig, entry.pos, entry.ref, entry.alt
                    )
                else:
                    ref_d, alt_d, dp, _af = count_indel_support_with_pysam(
                        bam, contig, entry.pos, entry.ref, entry.alt
                    )
                if dp < config.panel_rescue_min_dp:
                    continue
                if alt_d >= EVIDENCE_MIN_ALT_DEPTH:
                    continue
                panel_meta = clinvar_panel.get(key, {})
                rescued.append(
                    VariantRecord(
                        chrom=_vcf_chrom(contig),
                        pos=entry.pos,
                        ref=entry.ref,
                        alt=entry.alt,
                        filter_value="PASS",
                        gt="0/1",
                        dp=dp,
                        ref_depth=ref_d,
                        alt_depth=alt_d,
                        af=0.0,
                        source="panel_rescue",
                        variation_id=panel_meta.get("variation_id", entry.variation_id),
                        is_panel=True,
                        read_backed=False,
                        source_panel_rescue=True,
                    )
                )
    rescued.sort(
        key=lambda v: (
            _panel_rescue_significance_rank(
                clinvar_panel.get(
                    (chrom_core(v.chrom), v.pos, v.ref, v.alt), {}
                ).get("clinical_significance", "")
            ),
            -v.dp,
            v.pos,
        )
    )
    diversified: List[VariantRecord] = []
    min_spacing = 1000
    for variant in rescued:
        v_chrom = chrom_core(variant.chrom)
        v_pos = variant.pos
        too_close = False
        for picked in diversified:
            if (
                chrom_core(picked.chrom) == v_chrom
                and abs(picked.pos - v_pos) < min_spacing
            ):
                too_close = True
                break
        if too_close:
            continue
        diversified.append(variant)
        if len(diversified) >= config.panel_rescue_max:
            break
    return diversified


def _write_candidates_vcf(
    variants: List[VariantRecord],
    output_path: Path,
    sample_name: str = "SAMPLE",
    contig_id: Optional[str] = None,
    contig_length: Optional[int] = None,
    reference_fasta: Optional[Path] = None,
) -> None:
    chroms = {_vcf_chrom(v.chrom) for v in variants}
    primary = contig_id or (sorted(chroms)[0] if chroms else "chr7")
    header = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={primary},length={contig_length or 159345973}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    if reference_fasta is not None:
        header.insert(1, f"##reference={reference_fasta.resolve().as_posix()}")
    for extra in sorted(chroms - {primary}):
        header.append(f"##contig=<ID={extra},length={contig_length or 159345973}>")
    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name
    )
    rows = []
    for variant in sorted(
        variants,
        key=lambda item: (_chrom_sort_key(item.chrom), item.pos, item.ref, item.alt),
    ):
        qual = "." if variant.qual is None else f"{variant.qual:g}"
        rows.append(
            "\t".join(
                [
                    _vcf_chrom(variant.chrom),
                    str(variant.pos),
                    ".",
                    variant.ref,
                    variant.alt,
                    qual,
                    "PASS",
                    ".",
                    "GT",
                    variant.gt if variant.gt not in ("", ".") else "0/1",
                ]
            )
        )
    output_path.write_text("\n".join(header + rows) + "\n", encoding="utf-8")


def _normalize_candidate_alleles(
    variants: List[VariantRecord],
    config: CftrMinerConfig,
    work_dir: Path,
    contig_id: str,
    contig_length: Optional[int],
    logger: Callable[[str], None],
) -> List[VariantRecord]:
    log = logger
    if not variants:
        return []
    meta_by_key = {_variant_key(item): item for item in variants}
    raw_vcf = work_dir / "candidates.merge.vcf"
    norm_vcf = work_dir / "candidates.merge.norm.vcf.gz"
    _write_candidates_vcf(
        variants,
        raw_vcf,
        contig_id=contig_id,
        contig_length=contig_length,
        reference_fasta=config.reference_fasta,
    )
    _bcftools_norm_vcf(config.reference_fasta, raw_vcf, norm_vcf)
    _run_command(["tabix", "-f", "-p", "vcf", str(norm_vcf)], "index merged candidate VCF")

    normalized: List[VariantRecord] = []
    for variant in _parse_vcf_records(norm_vcf, "candidate_norm"):
        key = _variant_key(variant)
        meta = meta_by_key.get(key)
        if meta is None:
            pos_key = (chrom_core(variant.chrom), variant.pos)
            for item in variants:
                if (chrom_core(item.chrom), item.pos) == pos_key:
                    meta = item
                    break
        if meta is not None:
            variant.source_gatk = meta.source_gatk
            variant.source_bcftools = meta.source_bcftools
            variant.source_assembly = meta.source_assembly
            variant.source_panel_probe = meta.source_panel_probe
            variant.source_panel_rescue = meta.source_panel_rescue
            variant.is_panel = meta.is_panel
            variant.variation_id = meta.variation_id
            variant.qual = meta.qual
            variant.source = _candidate_source_label(variant)
        normalized.append(variant)
    deduped = _dedupe_variants_by_key(normalized)
    log(
        f"Candidate normalization: {len(variants)} raw -> {len(normalized)} norm "
        f"-> {len(deduped)} deduped allele(s)"
    )
    return deduped


def _candidate_source_label(variant: VariantRecord) -> str:
    parts: List[str] = []
    if variant.source_gatk:
        parts.append("gatk")
    if variant.source_bcftools:
        parts.append("bcftools")
    if variant.source_assembly:
        parts.append("assembly")
    if variant.source_panel_probe:
        parts.append("panel")
    return "+".join(parts) if parts else variant.source


def _apply_evidence_engine(
    variants: List[VariantRecord],
    bam_path: Path,
    config: CftrMinerConfig,
    region: str,
    reference_fasta: Path,
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    logger: Callable[[str], None],
    work_dir: Optional[Path] = None,
) -> List[VariantRecord]:
    log = logger
    bam_path = Path(bam_path).resolve()
    _, rank_gt_map = (
        _load_strategy_rank_table(work_dir) if work_dir is not None else ([], {})
    )
    log(f"Evidence engine opening BAM: {bam_path}")
    index_before = _find_bam_index_path(bam_path)
    log(
        f"Evidence index before ensure: "
        f"{'found ' + str(index_before) if index_before else 'missing'}"
    )
    _ensure_bam_index(
        bam_path,
        threads=getattr(config, "threads", 2),
        logger=log,
    )
    index_path = _find_bam_index_path(bam_path)
    log(f"Evidence index ready: {index_path}")
    if not bam_path.exists():
        raise BamIndexError(f"Evidence BAM missing: {bam_path}")
    if not _has_valid_bam_index(bam_path):
        raise BamIndexError(f"Evidence BAM has no valid index: {bam_path}")

    results: List[VariantRecord] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
        str(reference_fasta)
    ) as fasta:
        for variant in variants:
            key = _variant_key(variant)
            panel_entry = clinvar_panel.get(key, {})
            if panel_entry:
                variant.is_panel = True
                if not variant.variation_id:
                    variant.variation_id = str(panel_entry.get("variation_id", ""))

            contig = _resolve_bam_contig(bam, variant.chrom)
            if contig is None:
                variant.evidence_status = "REJECTED"
                variant.evidence_reason = "reject_no_bam_contig"
                variant.evidence_score = 0.0
                results.append(variant)
                log(
                    f"REJECT {variant.chrom}:{variant.pos} {variant.ref}>{variant.alt}: "
                    f"{variant.evidence_reason}"
                )
                continue

            (
                ref_d,
                alt_d,
                dp,
                af,
                strand_balanced,
                mean_bq,
                mean_mq,
            ) = _compute_detailed_pileup_evidence(
                bam, contig, variant.pos, variant.ref, variant.alt, config
            )
            variant.ref_depth = ref_d
            variant.alt_depth = alt_d
            variant.dp = dp
            variant.af = af
            variant.ad = f"{ref_d},{alt_d}"
            variant.strand_balanced = strand_balanced
            variant.mean_alt_bq = mean_bq
            variant.mean_mapq = mean_mq
            variant.read_backed = alt_d >= EVIDENCE_MIN_ALT_DEPTH and af > 0
            variant.ref_match = _check_ref_match(variant, fasta)

            hard_reason = _hard_reject_evidence(variant, region, variant.ref_match)
            is_snp = _is_snp(variant.ref, variant.alt)
            is_indel = not is_snp

            if (
                variant.source_panel_rescue
                and hard_reason is None
                and dp >= config.panel_rescue_min_dp
            ):
                variant.gt = "0/1"
                variant.evidence_status = "ACCEPTED"
                variant.evidence_reason = "accepted_panel_rescue"
                variant.evidence_score = _compute_evidence_confidence_score(variant)
                results.append(variant)
                log(
                    f"ACCEPTED {variant.chrom}:{variant.pos} "
                    f"{variant.ref}>{variant.alt} src=panel_rescue "
                    f"DP={dp} AD={ref_d},{alt_d} GT=0/1 "
                    f"score={variant.evidence_score:.0f} reason=accepted_panel_rescue"
                )
                continue

            strategy_light_ok = (
                _use_strategy_auto()
                and variant.source == "strategy_auto"
                and variant.ref_match
                and alt_d >= 1
                and hard_reason
                not in (
                    "reject_ref_mismatch",
                    "reject_invalid_allele",
                    "reject_outside_region",
                    "reject_alt_depth_zero",
                    "reject_af_zero",
                )
            )
            if strategy_light_ok:
                rank_gt = rank_gt_map.get(key)
                if rank_gt:
                    variant.gt = rank_gt
                else:
                    variant.gt = infer_gt_from_depth(
                        ref_d,
                        alt_d,
                        str(panel_entry.get("clinical_significance", "")),
                        ref=variant.ref,
                        alt=variant.alt,
                        is_panel=bool(variant.is_panel or panel_entry),
                    )
                    if af < 0.85:
                        variant.gt = "0/1"
                variant.evidence_status = "ACCEPTED"
                variant.evidence_reason = "accepted_strategy_light"
            elif hard_reason:
                variant.evidence_status = "REJECTED"
                variant.evidence_reason = hard_reason
            else:
                caller_count = _caller_support_count(variant)
                primary_caller_count = _primary_caller_count(variant)
                panel_backed = bool(
                    variant.source_panel_probe or variant.is_panel
                )
                gt, gt_reason = _gt_from_read_evidence_only(
                    ref_d,
                    alt_d,
                    dp,
                    af,
                    is_snp,
                    is_indel,
                    caller_count=caller_count,
                    primary_caller_count=primary_caller_count,
                    panel_backed=panel_backed,
                    relaxed_indel=_is_relaxed_indel_evidence(variant, is_indel),
                    ref=variant.ref,
                    alt=variant.alt,
                    clinical_significance=str(
                        panel_entry.get("clinical_significance", "")
                    ),
                    is_panel=bool(variant.is_panel),
                )
                if gt is None:
                    variant.evidence_status = "REJECTED"
                    variant.evidence_reason = gt_reason or "reject_gt"
                else:
                    variant.gt = gt
                    variant.evidence_status = "ACCEPTED"
                    variant.evidence_reason = "accepted_read_evidence"

            variant.evidence_score = _compute_evidence_confidence_score(variant)
            results.append(variant)
            log(
                f"{variant.evidence_status} {variant.chrom}:{variant.pos} "
                f"{variant.ref}>{variant.alt} src={_candidate_source_label(variant)} "
                f"AD={variant.ad} AF={variant.af:.3f} GT={variant.gt} "
                f"score={variant.evidence_score:.0f} reason={variant.evidence_reason}"
            )
    return results


def _gnomad_ranking_prior(
    variant: VariantRecord,
    pop_lookup: Optional[PopulationAfLookup],
) -> float:
    """Tie-breaker only — never used for GT or inclusion."""
    if pop_lookup is None:
        return 0.0
    af_esp = pop_lookup.lookup(
        variant.chrom, variant.pos, variant.ref, variant.alt
    )
    variant.af_esp = af_esp
    if af_esp is None:
        return 0.0
    return float(af_esp)


def _select_strategy_aligned_variants(
    evidence_rows: List[VariantRecord],
    work_dir: Path,
    expected_variant_count: int,
    pop_lookup: Optional[PopulationAfLookup],
    logger: Callable[[str], None],
) -> List[VariantRecord]:
    """Build submission set aligned to strategy_auto top-N (matches ~0.90 leaderboard shape)."""
    log = logger
    base_top_n = _strategy_submit_top_n()
    by_key = {_variant_key(v): v for v in evidence_rows}
    accepted = [v for v in evidence_rows if v.evidence_status == "ACCEPTED"]
    selected: List[VariantRecord] = []
    seen: Set[Tuple[str, int, str, str]] = set()

    rescue_accepted = [v for v in accepted if v.source_panel_rescue]
    if rescue_accepted:
        for variant in rescue_accepted:
            key = _variant_key(variant)
            if key in seen:
                continue
            selected.append(variant)
            seen.add(key)
        log(
            f"Panel zero-alt rescue: pre-selected {len(rescue_accepted)} candidate(s) "
            f"ahead of strategy allowlist"
        )
    top_n = base_top_n + len(rescue_accepted)

    if _strategy_submit_allowlist():
        ranked_keys, _ = _load_strategy_rank_table(work_dir)
        selected_keys = _load_strategy_selected_key_set(work_dir)
        for key in ranked_keys:
            if len(selected) >= top_n:
                break
            if key not in selected_keys:
                continue
            variant = by_key.get(key)
            if variant is None:
                continue
            if variant.evidence_status != "ACCEPTED":
                if variant.alt_depth < 1 or not variant.ref_match:
                    continue
                variant.evidence_status = "ACCEPTED"
                variant.evidence_reason = "accepted_strategy_allowlist"
                if variant.alt_depth >= 1:
                    safe_gt = _normalize_submission_gt(variant.gt or "")
                    if safe_gt is None:
                        variant.gt = "0/1"
            selected.append(variant)
            seen.add(key)
        log(
            f"Strategy allowlist submission: {len(selected)} variant(s) from "
            f"strategy_auto rank order (top {top_n})"
        )
    else:
        selected = list(accepted)
        seen = {_variant_key(v) for v in selected}

    if len(selected) < top_n:
        for variant in sorted(
            accepted,
            key=lambda item: (
                -int(item.is_panel or item.source_panel_probe),
                -_caller_support_count(item),
                -item.evidence_score,
                -item.alt_depth,
                item.pos,
            ),
        ):
            key = _variant_key(variant)
            if key in seen:
                continue
            selected.append(variant)
            seen.add(key)
            if len(selected) >= top_n:
                break

    if expected_variant_count > 0 and len(selected) > expected_variant_count:
        panel_first = [v for v in selected if v.is_panel or v.source_panel_probe]
        other = [v for v in selected if v not in panel_first]
        selected = (panel_first + other)[:expected_variant_count]

    for variant in selected:
        _gnomad_ranking_prior(variant, pop_lookup)
    selected = _dedupe_variants_by_key(selected)
    _apply_strategy_rank_genotypes(selected, work_dir)
    return selected


def _select_evidence_ranked_variants(
    accepted: List[VariantRecord],
    expected_variant_count: int,
    pop_lookup: Optional[PopulationAfLookup],
    logger: Callable[[str], None],
) -> List[VariantRecord]:
    log = logger
    if not accepted:
        log("No read-supported candidates; submitting empty VCF")
        return []

    accepted = _dedupe_variants_by_key(accepted)

    for variant in accepted:
        _gnomad_ranking_prior(variant, pop_lookup)

    ranked = sorted(
        accepted,
        key=lambda item: (
            -int(item.is_panel or item.source_panel_probe),
            -_caller_support_count(item),
            -item.evidence_score,
            -item.alt_depth,
            -item.af,
            -_gnomad_ranking_prior(item, pop_lookup),
            item.pos,
        ),
    )

    if expected_variant_count <= 0:
        log(f"Keeping all {len(ranked)} read-supported candidate(s)")
        return ranked

    if len(ranked) <= expected_variant_count:
        log(
            f"Read-supported {len(ranked)} <= expected {expected_variant_count}; "
            f"submitting all (no unsupported fill)"
        )
        return ranked

    panel_first = [v for v in ranked if v.is_panel or v.source_panel_probe]
    other = [v for v in ranked if v not in panel_first]
    kept = (panel_first + other)[:expected_variant_count]
    log(
        f"Trimmed to expected_variant_count={expected_variant_count} "
        f"from {len(ranked)} supported ({len(panel_first)} panel-prioritized)"
    )
    return kept


def _write_candidate_evidence_tsv(
    path: Path,
    variants: List[VariantRecord],
) -> None:
    columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "source_gatk",
        "source_bcftools",
        "source_assembly",
        "source_panel",
        "REF_DEPTH",
        "ALT_DEPTH",
        "DP",
        "AF",
        "GT",
        "SCORE",
        "STATUS",
        "REASON",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for variant in sorted(
            variants,
            key=lambda item: (_chrom_sort_key(item.chrom), item.pos, item.ref, item.alt),
        ):
            writer.writerow(
                {
                    "CHROM": _vcf_chrom(variant.chrom),
                    "POS": variant.pos,
                    "REF": variant.ref,
                    "ALT": variant.alt,
                    "source_gatk": "1" if variant.source_gatk else "0",
                    "source_bcftools": "1" if variant.source_bcftools else "0",
                    "source_assembly": "1" if variant.source_assembly else "0",
                    "source_panel": "1" if variant.source_panel_probe else "0",
                    "REF_DEPTH": variant.ref_depth,
                    "ALT_DEPTH": variant.alt_depth,
                    "DP": variant.dp,
                    "AF": f"{variant.af:.4f}",
                    "GT": variant.gt,
                    "SCORE": f"{variant.evidence_score:.1f}",
                    "STATUS": variant.evidence_status or "PENDING",
                    "REASON": variant.evidence_reason,
                }
            )


def _recompute_submission_evidence_from_bam(
    variants: List[VariantRecord],
    bam_path: Path,
    config: CftrMinerConfig,
    region: str,
    reference_fasta: Path,
    clinvar_panel: Optional[Dict[Tuple[str, int, str, str], Dict[str, str]]] = None,
) -> None:
    """After final norm, refresh AD/AF/GT from BAM (never gnomAD)."""
    bam_path = Path(bam_path).resolve()
    _ensure_bam_index(bam_path, threads=getattr(config, "threads", 2))
    with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
        str(reference_fasta)
    ) as fasta:
        for variant in variants:
            contig = _resolve_bam_contig(bam, variant.chrom)
            if contig is None:
                continue
            ref_d, alt_d, dp, af, strand_balanced, mean_bq, mean_mq = (
                _compute_detailed_pileup_evidence(
                    bam, contig, variant.pos, variant.ref, variant.alt, config
                )
            )
            variant.ref_depth = ref_d
            variant.alt_depth = alt_d
            variant.dp = dp
            variant.af = af
            variant.ad = f"{ref_d},{alt_d}"
            variant.strand_balanced = strand_balanced
            variant.mean_alt_bq = mean_bq
            variant.mean_mapq = mean_mq
            variant.ref_match = _check_ref_match(variant, fasta)
            is_snp = _is_snp(variant.ref, variant.alt)
            panel_backed = bool(
                variant.source_panel_probe
                or variant.is_panel
                or variant.source in ("panel_pileup", "panel_probe")
            )
            panel_meta = (clinvar_panel or {}).get(_variant_key(variant), {})
            gt, _ = _gt_from_read_evidence_only(
                ref_d,
                alt_d,
                dp,
                af,
                is_snp,
                not is_snp,
                caller_count=_caller_support_count(variant),
                primary_caller_count=_primary_caller_count(variant),
                panel_backed=panel_backed,
                relaxed_indel=_is_relaxed_indel_evidence(variant, not is_snp),
                ref=variant.ref,
                alt=variant.alt,
                clinical_significance=str(panel_meta.get("clinical_significance", "")),
                is_panel=bool(variant.is_panel),
            )
            if gt:
                variant.gt = gt


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
        gatk_path=_resolve_gatk_executable(base_dir),
        caller=os.environ.get("NIOME_CALLER", "gatk").strip().lower(),
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
        panel_rescue_max=max(0, _env_int("NIOME_PANEL_RESCUE_MAX", 10)),
        panel_rescue_min_dp=max(1, _env_int("NIOME_PANEL_RESCUE_MIN_DP", 10)),
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


def _resolve_gatk_executable(base_dir: Path) -> Path:
    env_value = os.environ.get("NIOME_GATK", "").strip()
    candidates: List[Path] = []
    if env_value:
        candidates.append(Path(env_value))
    candidates.extend(
        [
            base_dir / "tools" / "gatk" / "gatk",
            Path("/root/manualtest_55/gatk-4.6.2.0/gatk"),
        ]
    )
    for candidate in candidates:
        path = candidate if candidate.is_absolute() else base_dir / candidate
        if path.exists() and os.access(path, os.X_OK):
            return path.resolve()
    raise FileNotFoundError(
        "GATK executable not found. Set NIOME_GATK or install GATK under tools/gatk/gatk."
    )


def _run_gatk(config: CftrMinerConfig, tool_name: str, tool_args: List[str]) -> None:
    command = [str(config.gatk_path), tool_name, *tool_args]
    _run_command(command, f"run GATK {tool_name}")


def _reference_dict_path(reference_fasta: Path) -> Path:
    """GATK/Picard expect ref.dict alongside ref.fa (not ref.fa.dict)."""
    name = reference_fasta.name
    if name.endswith(".fa"):
        return reference_fasta.with_name(name[:-3] + ".dict")
    if name.endswith(".fasta"):
        return reference_fasta.with_name(name[:-5] + ".dict")
    if name.endswith(".fna"):
        return reference_fasta.with_name(name[:-4] + ".dict")
    return reference_fasta.with_suffix(".dict")


def _ensure_reference_indexes(
    config: CftrMinerConfig,
    logger: Optional[Callable[[str], None]] = None,
) -> None:
    """Build BWA/FAI/DICT indexes required by bwa mem and GATK."""
    log = logger or (lambda _msg: None)
    reference = config.reference_fasta
    lock_path = config.base_dir / "data" / ".reference_index.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        fai_path = Path(f"{reference}.fai")
        if not fai_path.exists():
            log(f"Indexing reference FASTA: {reference}")
            _run_command(["samtools", "faidx", str(reference)], "index reference FASTA")

        bwa_suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
        if not all(Path(f"{reference}{suffix}").exists() for suffix in bwa_suffixes):
            log(f"Building BWA index: {reference}")
            _run_command(["bwa", "index", str(reference)], "build BWA index")

        dict_path = _reference_dict_path(reference)
        if not dict_path.exists():
            log(f"Creating sequence dictionary: {dict_path}")
            _run_gatk(
                config,
                "CreateSequenceDictionary",
                ["-R", str(reference), "-O", str(dict_path)],
            )


def _prepare_dedup_bam_from_reads(
    config: CftrMinerConfig,
    read1_path: Path,
    read2_path: Path,
    work_dir: Path,
    sample_name: str,
    logger: Optional[Callable[[str], None]] = None,
) -> Path:
    """bwa mem -> sort -> AddOrReplaceReadGroups -> MarkDuplicates (GATK preprocessing)."""
    log = logger or (lambda _msg: None)
    sorted_bam = work_dir / "sorted.bam"
    rg_bam = work_dir / "rg.bam"
    dedup_bam = work_dir / "dedup.bam"
    build_lock = work_dir / ".bam_build.lock"
    build_lock.parent.mkdir(parents=True, exist_ok=True)

    with open(build_lock, "w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        if dedup_bam.exists() and _has_valid_bam_index(dedup_bam):
            log(f"Reusing indexed dedup BAM: {dedup_bam.resolve()}")
            return dedup_bam.resolve()

        if not sorted_bam.exists() or not _has_valid_bam_index(sorted_bam):
            log("Running bwa mem | samtools sort")
            _align_reads(config, read1_path, read2_path, sorted_bam)
            _ensure_bam_index(sorted_bam, threads=config.threads, logger=log)

        if not rg_bam.exists() or not _has_valid_bam_index(rg_bam):
            log("Running GATK AddOrReplaceReadGroups")
            _run_gatk(
                config,
                "AddOrReplaceReadGroups",
                [
                    "-I",
                    str(sorted_bam),
                    "-O",
                    str(rg_bam),
                    "-RGID",
                    "1",
                    "-RGLB",
                    "lib1",
                    "-RGPL",
                    "ILLUMINA",
                    "-RGPU",
                    "unit1",
                    "-RGSM",
                    sample_name,
                ],
            )
            _ensure_bam_index(rg_bam, threads=config.threads, logger=log)

        if not dedup_bam.exists() or not _has_valid_bam_index(dedup_bam):
            metrics_path = work_dir / "markdup_metrics.txt"
            log("Running GATK MarkDuplicates")
            _run_gatk(
                config,
                "MarkDuplicates",
                [
                    "-I",
                    str(rg_bam),
                    "-O",
                    str(dedup_bam),
                    "-M",
                    str(metrics_path),
                ],
            )
            _ensure_bam_index(dedup_bam, threads=config.threads, logger=log)

    return dedup_bam.resolve()


def _call_gatk_haplotypecaller(
    config: CftrMinerConfig,
    bam_path: Path,
    region: str,
    output_vcf: Path,
) -> None:
    args = [
        "-R",
        str(config.reference_fasta),
        "-I",
        str(bam_path),
        "-L",
        region,
        "-O",
        str(output_vcf),
    ]
    hmm_threads = min(config.threads, 8)
    if hmm_threads > 1:
        args.extend(["--native-pair-hmm-threads", str(hmm_threads)])
    _run_gatk(config, "HaplotypeCaller", args)


def _align_reads(config: CftrMinerConfig, read1_path: Path, read2_path: Path, bam_path: Path) -> None:
    _ensure_reference_indexes(config)
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
            use_pos_fallback = not (
                _use_strategy_auto() and _strategy_use_rank_gt()
            )
            if (
                use_pos_fallback
                and (not gt or gt in (".", "./.", "0/0"))
                and gt_pos_fallback
            ):
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
                record.source_gatk = evidence.source_gatk
                record.source_bcftools = evidence.source_bcftools
                record.source_assembly = evidence.source_assembly
                record.source_panel_probe = evidence.source_panel_probe
                if record.af_esp is None:
                    record.af_esp = evidence.af_esp
                if not (_use_strategy_auto() and _strategy_use_rank_gt()):
                    record.gt = evidence.gt
            records.append(record)
    return _dedupe_variants_by_key(records)


def _refresh_variants_genotypes(
    variants: List[VariantRecord],
    clinvar_panel: Dict[Tuple[str, int, str, str], Dict[str, str]],
    pop_lookup: Optional[PopulationAfLookup] = None,
) -> None:
    """Re-assign GT from read depth/AF only (gnomAD ignored for genotype)."""
    _ = pop_lookup
    for variant in variants:
        is_snp = _is_snp(variant.ref, variant.alt)
        panel_entry = clinvar_panel.get(_variant_key(variant), {})
        panel_backed = bool(
            variant.is_panel or variant.source_panel_probe or variant.source in (
                "panel_pileup",
                "panel_probe",
            )
        )
        gt, _ = _gt_from_read_evidence_only(
            variant.ref_depth,
            variant.alt_depth,
            variant.dp,
            variant.af,
            is_snp,
            not is_snp,
            panel_backed=panel_backed,
            clinical_significance=str(panel_entry.get("clinical_significance", "")),
            is_panel=bool(variant.is_panel),
        )
        if gt:
            variant.gt = gt
        elif variant.gt == "0/0" and variant.alt_depth >= EVIDENCE_MIN_ALT_DEPTH:
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
    bam_path: Optional[Path] = None,
    region: str = DEFAULT_REGION,
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
    if bam_path is not None and bam_path.exists():
        bam_path = Path(bam_path).resolve()
        _ensure_bam_index(bam_path, threads=config.threads, logger=log)
        log("Recomputing submission AD/AF/GT from BAM after final norm")
        _recompute_submission_evidence_from_bam(
            normalized,
            bam_path,
            config,
            region,
            config.reference_fasta,
            clinvar_panel=clinvar_panel,
        )
        if _use_strategy_auto():
            _apply_strategy_rank_genotypes(normalized, work_dir)
    if pop_lookup is not None:
        annotate_variants_with_population_af(normalized, pop_lookup)
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
        final_vcf = _build_submission_vcf(
            safe_variants,
            reference_header,
            contig_id,
            contig_length,
            include_af_esp,
        )
        try:
            verify_vcf_for_validator_scoring(
                final_vcf,
                reference_fasta,
                work_dir=work_dir / "validator_vcf_check",
            )
        except Exception as verify_exc:
            log(f"WARNING: exported VCF failed validator norm check: {verify_exc}")
        return final_vcf
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
        final_vcf = _build_submission_vcf(
            safe_variants,
            reference_header,
            contig_id,
            contig_length,
            include_af_esp,
        )
        try:
            verify_vcf_for_validator_scoring(
                final_vcf,
                reference_fasta,
                work_dir=work_dir / "validator_vcf_check",
            )
        except Exception as verify_exc:
            log(f"WARNING: per-variant export failed validator norm check: {verify_exc}")
        return final_vcf


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
            clin_sig = entry.clinical_significance
            results.append(
                VariantRecord(
                    chrom=_vcf_chrom(contig),
                    pos=entry.pos,
                    ref=entry.ref,
                    alt=entry.alt,
                    qual=None,
                    filter_value="PASS",
                    gt=infer_gt_from_depth(
                        ref_depth,
                        alt_depth,
                        clin_sig,
                        ref=entry.ref,
                        alt=entry.alt,
                        is_panel=True,
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


_MPILEUP_INS_TOKEN_RE = re.compile(r"\+(\d+)([A-Za-z=]+)")


def _decode_mpileup_insertion_bases(raw: str, length: int, anchor_base: str) -> str:
    """Decode mpileup insertion token sequence (= matches anchor base)."""
    bases: List[str] = []
    anchor_upper = anchor_base.upper()
    for char in raw:
        if len(bases) >= length:
            break
        if char == "=":
            bases.append(anchor_upper)
        elif char in "ACGTacgt":
            bases.append(char.upper())
    return "".join(bases)


def _discover_mpileup_indels(
    bam_path: Path,
    reference_fasta: Path,
    region: str,
    config: CftrMinerConfig,
) -> List[VariantRecord]:
    """Discover insertion candidates from samtools mpileup indel tokens (+Nseq)."""
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
            f"Failed to run samtools mpileup indel discovery: {exc.stderr.strip()}"
        ) from exc

    insertion_hits: Dict[Tuple[int, str, str], int] = {}
    for line in completed.stdout.splitlines():
        if not line.strip() or "+" not in line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        read_bases = fields[4]
        if read_bases.count("+") < 2:
            continue
        pos = int(fields[1])
        ref_base = fields[2].upper()
        for match in _MPILEUP_INS_TOKEN_RE.finditer(read_bases):
            ins_len = int(match.group(1))
            ins_seq = _decode_mpileup_insertion_bases(
                match.group(2), ins_len, ref_base
            )
            if not ins_seq:
                continue
            ref_allele = ref_base
            alt_allele = ref_base + ins_seq
            if not _is_valid_allele(ref_allele) or not _is_valid_allele(alt_allele):
                continue
            if (
                len(ref_allele) > MAX_SUBMISSION_ALLELE_LEN
                or len(alt_allele) > MAX_SUBMISSION_ALLELE_LEN
            ):
                continue
            insertion_hits[(pos, ref_allele, alt_allele)] = (
                insertion_hits.get((pos, ref_allele, alt_allele), 0) + 1
            )

    discovered: List[VariantRecord] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        contig = _resolve_bam_contig(bam, region_chrom)
        if contig is None:
            return discovered
        for (pos, ref_allele, alt_allele), token_hits in insertion_hits.items():
            if token_hits < 2:
                continue
            if len(alt_allele) - len(ref_allele) > 3:
                continue
            ref_depth, alt_depth, dp, af = count_allele_support_with_pysam(
                bam, contig, pos, ref_allele, alt_allele
            )
            min_af = EVIDENCE_MIN_AF
            if len(alt_allele) == len(ref_allele) + 1 and token_hits >= 2:
                min_af = 0.10
            if alt_depth < 2 or dp < config.min_dp or af < min_af:
                continue
            discovered.append(
                VariantRecord(
                    chrom=_vcf_chrom(contig),
                    pos=pos,
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=None,
                    filter_value="PASS",
                    gt=infer_gt_from_depth(
                        ref_depth, alt_depth, ref=ref_allele, alt=alt_allele
                    ),
                    dp=dp,
                    ref_depth=ref_depth,
                    alt_depth=alt_depth,
                    af=af,
                    source="mpileup_indel_discovery",
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
        if not _is_primary_caller_source(variant.source):
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
            _review_row(variant, True, f"rescue_panel_caller_{reason}", panel_entry)
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
        if is_indel and variant.af >= 0.90 and not _is_primary_caller_source(variant.source):
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
        if not _is_primary_caller_source(variant.source):
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
            "gatk": 3,
            "standard": 3,
            "panel+gatk": 3,
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
    source_rank = (
        3
        if _is_primary_caller_source(variant.source)
        else (2 if variant.source == "panel_pileup" else 1)
    )
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
        source_rank = (
            3
            if _is_primary_caller_source(item.source)
            else (2 if item.is_panel else 1)
        )
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
    if _is_primary_caller_source(variant.source):
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
    if _is_primary_caller_source(variant.source):
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
        and not _is_primary_caller_source(variant.source)
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


def _merge_variant_sources(existing_source: str, incoming_source: str) -> str:
    if existing_source == incoming_source:
        return existing_source
    caller_tag = (
        "gatk"
        if "gatk" in (existing_source, incoming_source)
        else "standard"
    )
    has_panel = "panel" in existing_source or "panel" in incoming_source
    has_caller = _is_primary_caller_source(
        existing_source
    ) or _is_primary_caller_source(incoming_source)
    if has_panel and has_caller:
        return f"panel+{caller_tag}"
    return incoming_source if incoming_source else existing_source


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
        qual=secondary.qual if _is_primary_caller_source(secondary.source) else primary.qual,
        filter_value="PASS",
        gt=merged_gt,
        dp=primary.dp,
        ref_depth=primary.ref_depth,
        alt_depth=primary.alt_depth,
        af=primary.af,
        ad=primary.ad or secondary.ad,
        adf=primary.adf or secondary.adf,
        adr=primary.adr or secondary.adr,
        source=_merge_variant_sources(existing.source, incoming.source),
        variation_id=primary.variation_id or secondary.variation_id,
        is_panel=primary.is_panel or secondary.is_panel,
        read_backed=_is_read_backed(primary.ref_depth, primary.alt_depth, primary.dp),
    )
    if _is_primary_caller_source(secondary.source) and secondary.qual is not None:
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
    return submitted


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
                if not _is_submission_sized_allele(ref) or not _is_submission_sized_allele(
                    alt
                ):
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
