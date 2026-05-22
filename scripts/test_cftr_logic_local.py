#!/usr/bin/env python3
"""Run the CFTR miner workflow against a local task JSON file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from niome_subnet.cftr_miner_logic import process_cftr_task_for_miner


def main() -> None:
    parser = argparse.ArgumentParser(description="Test CFTR miner logic locally")
    parser.add_argument("task_json", help="Path to a validator task JSON file")
    parser.add_argument(
        "--base-dir",
        default=".",
        help="Project base directory containing data/, panel/, reads/, work/, and outputs/",
    )
    args = parser.parse_args()

    task_path = Path(args.task_json)
    task_data = json.loads(task_path.read_text(encoding="utf-8"))
    result = process_cftr_task_for_miner(task_data, base_dir=args.base_dir, logger=print)

    print("\nVCF content")
    print(result["vcf_content"])

    print("cftr_annotations JSON")
    print(json.dumps(result["cftr_annotations"], indent=2, sort_keys=True))

    print("\ncounts")
    print(json.dumps(result["counts"], indent=2, sort_keys=True))

    print("\noutput paths")
    print(json.dumps(result["paths"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
