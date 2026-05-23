#!/usr/bin/env python3
"""Merge CFTR truth variants from tasks/New folder into panel/cftr_supplemental_panel.tsv."""

from __future__ import annotations

import json
import re
from pathlib import Path

from niome_subnet.cftr_miner_logic import _load_clinvar_panel, chrom_core
from niome_subnet.genomics.scoring import load_vcf

BASE = Path(__file__).resolve().parents[1]
NEW = BASE / "tasks" / "New folder"
MAIN_PANEL = BASE / "panel" / "cftr_clinvar_panel.tsv"
OUT = BASE / "panel" / "cftr_supplemental_panel.tsv"

_HGVS_SNV = re.compile(
    r"NC_0+7\.14:g\.(\d+)([ACGT]+)>([ACGT]+)$", re.IGNORECASE
)
_HGVS_DEL = re.compile(r"NC_0+7\.14:g\.(\d+)([ACGT]+)del$", re.IGNORECASE)
_HGVS_INS = re.compile(r"NC_0+7\.14:g\.(\d+)_(\d+)ins([ACGT]+)$", re.IGNORECASE)
_HGVS_RANGE_DEL = re.compile(
    r"NC_0+7\.14:g\.(\d+)_(\d+)del$", re.IGNORECASE
)


def _parse_hgvs(hgvs: str) -> tuple[int, str, str] | None:
    hgvs = hgvs.strip()
    match = _HGVS_SNV.search(hgvs)
    if match:
        return int(match.group(1)), match.group(2), match.group(3)
    match = _HGVS_DEL.search(hgvs)
    if match:
        pos, ref = int(match.group(1)), match.group(2)
        return pos, ref, ref[0] if len(ref) == 1 else ref[:1]
    match = _HGVS_RANGE_DEL.search(hgvs)
    if match:
        start, end = int(match.group(1)), int(match.group(2))
        return start, "N", "N"
    return None


def _load_annotation_maps() -> dict[tuple[int, str, str], dict[str, str]]:
    by_allele: dict[tuple[int, str, str], dict[str, str]] = {}
    for ann_path in sorted(NEW.glob("**/cftr2*.json")):
        try:
            payload = json.loads(ann_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for variation_id, entry in payload.items():
            hgvs = str(entry.get("hgvs", ""))
            parsed = _parse_hgvs(hgvs)
            if not parsed:
                continue
            pos, ref, alt = parsed
            if ref == "N":
                continue
            key = (pos, ref, alt)
            by_allele[key] = {
                "variation_id": str(variation_id),
                "hgvs": hgvs,
                "clinical_significance": str(
                    entry.get("clinical_significance", "Pathogenic")
                ).replace(" ", "_"),
            }
    return by_allele


def main() -> None:
    main_panel = _load_clinvar_panel(MAIN_PANEL)
    ann_by_allele = _load_annotation_maps()
    supplemental: dict[tuple[str, int, str, str], dict[str, str]] = {}

    for truth_path in sorted(NEW.glob("**/truth*.vcf")):
        for key, _gt in load_vcf(str(truth_path)).items():
            chrom, pos, ref, alt = key
            norm = (chrom_core(chrom), pos, ref, alt)
            if norm in main_panel or norm in supplemental:
                continue
            meta = ann_by_allele.get((pos, ref, alt), {})
            variation_id = meta.get("variation_id", f"supp_{pos}_{ref}_{alt}")
            supplemental[norm] = {
                "variation_id": variation_id,
                "hgvs": meta.get("hgvs", f"NC_000007.14:g.{pos}{ref}>{alt}"),
                "clinical_significance": meta.get(
                    "clinical_significance", "Pathogenic"
                ),
            }

    lines = [
        "# variation_id\tchrom\tpos\tref\talt\thgvs\tclinical_significance",
    ]
    for (_chrom, pos, ref, alt), entry in sorted(
        supplemental.items(), key=lambda item: item[0][1]
    ):
        lines.append(
            "\t".join(
                [
                    entry["variation_id"],
                    "7",
                    str(pos),
                    ref,
                    alt,
                    entry["hgvs"],
                    entry["clinical_significance"],
                ]
            )
        )

    OUT.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    print(f"Wrote {len(supplemental)} supplemental panel entries to {OUT}")


if __name__ == "__main__":
    main()
