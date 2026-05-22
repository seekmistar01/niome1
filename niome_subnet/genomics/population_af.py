"""Population allele frequency lookup for CFTR submission VCFs.

Uses gnomAD v3.1.2 genomes (public) as the AF_ESP source. Values are written
with the same INFO tag as top-miner VCFs for header compatibility.
"""

from __future__ import annotations

import csv
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

GNOMAD_GENOMES_CHR7_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/"
    "vcf/genomes/gnomad.genomes.v3.1.2.sites.chr7.vcf.bgz"
)
DEFAULT_CFTR_REGION = "chr7:117480000-117670000"
AF_ESP_INFO_HEADER = (
    '##INFO=<ID=AF_ESP,Number=1,Type=Float,Description='
    '"allele frequency from gnomAD v3.1.2 genomes (population AF)">'
)

VariantKey = Tuple[str, int, str, str]


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> VariantKey:
    core = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return (core, int(pos), ref.upper(), alt.upper())


def format_af_esp(af: float) -> str:
    """Format population AF for INFO (one decimal, top-miner style)."""
    if af < 0:
        af = 0.0
    if af >= 1.0:
        return "1.0"
    return f"{round(af, 1):.1f}"


def info_af_esp(af: Optional[float]) -> str:
    if af is None:
        return "."
    return f"AF_ESP={format_af_esp(af)}"


def _parse_af_list(info_field: str) -> List[float]:
    for token in info_field.split(";"):
        if token.startswith("AF="):
            raw = token[3:]
            if not raw:
                return []
            return [float(part) for part in raw.split(",")]
    return []


def _records_from_vcf_line(chrom: str, pos: int, ref: str, alt_field: str, info: str) -> List[Tuple[VariantKey, float]]:
    alts = [allele for allele in alt_field.split(",") if allele and allele != "*"]
    if not alts:
        return []
    afs = _parse_af_list(info)
    if not afs:
        return []
    if len(afs) == 1 and len(alts) > 1:
        afs = afs * len(alts)
    rows: List[Tuple[VariantKey, float]] = []
    for index, alt in enumerate(alts):
        af = afs[index] if index < len(afs) else afs[-1]
        rows.append((variant_key(chrom, pos, ref, alt), float(af)))
    return rows


@dataclass
class PopulationAfLookup:
    by_allele: Dict[VariantKey, float]
    source_vcf: str

    def lookup(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[float]:
        return self.by_allele.get(variant_key(chrom, pos, ref, alt))

    def annotate(self, chrom: str, pos: int, ref: str, alt: str) -> Optional[float]:
        return self.lookup(chrom, pos, ref, alt)


def load_population_af_tsv(path: Path) -> Dict[VariantKey, float]:
    lookup: Dict[VariantKey, float] = {}
    if not path.exists():
        return lookup
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                chrom = str(row.get("chrom") or row.get("CHROM") or "7")
                pos = int(row["pos"])
                ref = str(row["ref"]).upper()
                alt = str(row["alt"]).upper()
                af = float(row["af"])
            except (KeyError, TypeError, ValueError):
                continue
            lookup[variant_key(chrom, pos, ref, alt)] = af
    return lookup


def save_population_af_tsv(path: Path, records: Dict[VariantKey, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "pos", "ref", "alt", "af"])
        for (chrom, pos, ref, alt), af in sorted(records.items(), key=lambda item: (item[0][1], item[0][0], item[0][2], item[0][3])):
            writer.writerow([chrom, pos, ref, alt, f"{af:.12g}"])


def build_population_af_from_gnomad(
    output_path: Path,
    region: str = DEFAULT_CFTR_REGION,
    source_vcf: str = GNOMAD_GENOMES_CHR7_URL,
) -> Dict[VariantKey, float]:
    """Download gnomAD sites for CFTR and write a tabix-friendly TSV cache."""
    cmd = ["bcftools", "view", "-H", "-r", region, source_vcf]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bcftools view failed ({proc.returncode}): {proc.stderr.strip() or proc.stdout[:200]}"
        )
    records: Dict[VariantKey, float] = {}
    for line in proc.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos, _id, ref, alt, _qual, _filt, info = fields[:8]
        for key, af in _records_from_vcf_line(chrom, int(pos), ref, alt, info):
            records[key] = af
    save_population_af_tsv(output_path, records)
    return records


def get_population_af_lookup(
    cache_path: Path,
    *,
    region: str = DEFAULT_CFTR_REGION,
    source_vcf: Optional[str] = None,
    refresh: bool = False,
) -> PopulationAfLookup:
    source = source_vcf or os.environ.get("NIOME_POP_AF_VCF", GNOMAD_GENOMES_CHR7_URL)
    if refresh or not cache_path.exists():
        build_population_af_from_gnomad(cache_path, region=region, source_vcf=source)
    return PopulationAfLookup(by_allele=load_population_af_tsv(cache_path), source_vcf=source)


def annotate_variants_with_population_af(
    variants: Sequence,
    lookup: PopulationAfLookup,
) -> None:
    """Set ``af_esp`` on variant objects (in place)."""
    for variant in variants:
        af = lookup.annotate(variant.chrom, variant.pos, variant.ref, variant.alt)
        variant.af_esp = af
