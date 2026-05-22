#!/usr/bin/env python3
"""Build local gnomAD population AF cache for the CFTR region."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from niome_subnet.genomics.population_af import (  # noqa: E402
    DEFAULT_CFTR_REGION,
    GNOMAD_GENOMES_CHR7_URL,
    build_population_af_from_gnomad,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data" / "cftr_population_af.tsv",
        help="Output TSV path (default: data/cftr_population_af.tsv)",
    )
    parser.add_argument(
        "--region",
        default=DEFAULT_CFTR_REGION,
        help=f"genomic region for bcftools -r (default: {DEFAULT_CFTR_REGION})",
    )
    parser.add_argument(
        "--source-vcf",
        default=GNOMAD_GENOMES_CHR7_URL,
        help="gnomAD VCF URL or local path",
    )
    args = parser.parse_args()
    records = build_population_af_from_gnomad(
        args.output,
        region=args.region,
        source_vcf=args.source_vcf,
    )
    print(f"Wrote {len(records)} allele entries to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
