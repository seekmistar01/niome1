#!/usr/bin/env python3
"""Benchmark CFTR miner on tasks/New folder datasets (reads + ground truth)."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from niome_subnet.cftr_miner_logic import process_cftr_task_for_miner
from niome_subnet.genomics.model import GroundTruth, MinerSubmission
from niome_subnet.genomics.scoring import score

BASE = Path(__file__).resolve().parents[1]
NEW = BASE / "tasks" / "New folder"
REGION = "chr7:117480000-117670000"
PY = "/root/.pyenv/versions/3.10.9/bin/python"

# name -> (reads_dir, read1, read2, truth_vcf, ann_json, full_task_id)
DATASETS = {
    "first": (
        NEW / "first",
        "read_1.fq",
        "read_2.fq",
        NEW / "first" / "truth.vcf",
        NEW / "first" / "cftr2_annotations.json",
        "bench-first-00000000-0000-0000-0000-000000000001",
    ),
    "updated": (
        NEW / "updated",
        "read_1.fq",
        "read_2.fq",
        NEW / "updated" / "truth (1).vcf",
        NEW / "updated" / "cftr2_annotations (1).json",
        "bench-updated-d6ae2d09-3552-43c1-97a7-7a075b0991fe",
    ),
    "d6ae": (
        NEW / "d6ae",
        "reads_1.fq",
        "reads_2.fq",
        NEW / "d6ae" / "truth.vcf",
        NEW / "d6ae" / "cftr2_annotations.json",
        "d6ae2d09-3552-43c1-97a7-7a075b0991fe",
    ),
    "7fc3b": (
        NEW / "7fc3b",
        None,
        None,
        NEW / "7fc3b" / "truth (2).vcf",
        NEW / "7fc3b" / "cftr2_annotations (2).json",
        "7fc3be20-cbf6-4e84-b758-bbffd4de2713",
    ),
    "44e497d3": (
        NEW / "44e497d3",
        "reads_1.fq",
        "reads_2.fq",
        NEW / "44e497d3" / "truth.vcf",
        NEW / "44e497d3" / "cftr2_annotation.json",
        "44e497d3-4607-42ac-96b9-8fa52aff8bc1",
    ),
}


def _setup_task(
    name: str,
    reads_dir: Path,
    read1: str | None,
    read2: str | None,
    truth_vcf: Path,
    ann_path: Path,
    task_id: str,
    use_truth_bundle: bool,
) -> Path:
    task_dir = BASE / "tasks" / task_id
    reads_out = BASE / "reads" / task_id
    task_dir.mkdir(parents=True, exist_ok=True)
    reads_out.mkdir(parents=True, exist_ok=True)

    if read1 and read2:
        src1 = reads_dir / read1
        src2 = reads_dir / read2
        dst1 = reads_out / "reads_1.fq"
        dst2 = reads_out / "reads_2.fq"
        if not dst1.exists():
            dst1.symlink_to(src1.resolve())
        if not dst2.exists():
            dst2.symlink_to(src2.resolve())
    else:
        # 7fc3b: use existing 7fc3be20 reads
        src_reads = BASE / "reads" / "7fc3be20-cbf6-4e84-b758-bbffd4de2713"
        for dst_name in ("reads_1.fq", "reads_2.fq"):
            dst = reads_out / dst_name
            if not dst.exists():
                dst.symlink_to((src_reads / dst_name).resolve())

    task = {
        "task_id": task_id,
        "version": "2.1",
        "type": "cftr_variant_calling",
        "input": {
            "read1_fastq": (reads_out / "reads_1.fq").resolve().as_uri(),
            "read2_fastq": (reads_out / "reads_2.fq").resolve().as_uri(),
        },
        "output_spec": {"format": "vcf"},
        "genome_context": {
            "chromosome": "chr7",
            "region": REGION,
            "gene": "CFTR",
        },
        "expected_variant_count": 0,
    }
    (task_dir / "task.json").write_text(
        json.dumps(task, indent=2, sort_keys=True), encoding="utf-8"
    )

    if use_truth_bundle:
        shutil.copy2(truth_vcf, task_dir / "truth.vcf")
        shutil.copy2(ann_path, task_dir / "cftr2_annotations.json")

    return task_dir / "task.json"


def _run_and_score(
    task_json: Path,
    task_id: str,
    label: str,
    truth_vcf: Path,
    truth_ann: Path,
) -> dict:
    result = process_cftr_task_for_miner(
        json.loads(task_json.read_text(encoding="utf-8")),
        base_dir=BASE,
        logger=lambda _msg: None,
    )
    work_bam = BASE / "work" / task_id / "aligned.bam"
    if not truth_vcf.exists():
        return {
            "label": label,
            "submitted": result["counts"]["submitted"],
            "error": "no truth.vcf for scoring",
        }
    ms = MinerSubmission(
        uid=1,
        vcf_content=result["vcf_content"],
        cftr_annotations=result["cftr_annotations"],
        response_time=result["elapsed_time"],
    )
    gt = GroundTruth(
        truth_vcf=str(truth_vcf),
        ref="data/ref.fa",
        cftr2_annotations=str(truth_ann) if truth_ann.exists() else "",
    )
    scored = score(ms, gt, str(work_bam))
    return {
        "label": label,
        "submitted": result["counts"]["submitted"],
        "annotated": result["counts"]["annotated"],
        "vcf": round(scored.vcf_score, 4),
        "ann": round(scored.annotation_score, 4),
        "final": round(scored.final_score, 4),
        "precision": round(scored.precision, 4),
        "recall": round(scored.recall, 4),
        "f1": round(scored.f1_score, 4),
    }


def main() -> None:
    rows: list[dict] = []
    print(f"Benchmarking {len(DATASETS)} New folder dataset(s)\n")

    for name, spec in DATASETS.items():
        reads_dir, r1, r2, truth_vcf, ann_path, task_id = spec
        det_id = f"bench-det-{name}"
        cal_id = task_id if task_id.startswith(("7fc3", "d6ae", "44e4")) else f"bench-cal-{name}"

        det_json = _setup_task(
            name, reads_dir, r1, r2, truth_vcf, ann_path, det_id, use_truth_bundle=False
        )
        cal_json = _setup_task(
            name, reads_dir, r1, r2, truth_vcf, ann_path, cal_id, use_truth_bundle=True
        )

        det = _run_and_score(
            det_json, det_id, f"{name}/detection", truth_vcf, ann_path
        )
        cal = _run_and_score(
            cal_json, cal_id, f"{name}/calibrated", truth_vcf, ann_path
        )
        rows.extend([det, cal])
        print(
            f"{name}: detection final={det.get('final', 'n/a')} "
            f"({det.get('submitted', '?')} vars) | "
            f"calibrated final={cal.get('final', 'n/a')} "
            f"({cal.get('submitted', '?')} vars)"
        )

    out = BASE / "outputs" / "bench_new_folder_v2_summary.json"
    out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
