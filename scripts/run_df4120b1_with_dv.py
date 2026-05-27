#!/usr/bin/env python3
"""Re-run task df4120b1 with NIOME_STRATEGY_DEEPVARIANT=1 using cached local FASTQs."""

import json
import os
import sys
from pathlib import Path

ROOT = Path("/root/55miner/niome1")
sys.path.insert(0, str(ROOT))

os.environ.setdefault("PYTHONPATH", str(ROOT))
os.environ["NIOME_USE_STRATEGY_AUTO"] = "1"
os.environ["NIOME_STRATEGY_DEEPVARIANT"] = "1"
os.environ["NIOME_STRATEGY_SENSITIVITY"] = "sensitive"
os.environ["NIOME_STRATEGY_SUBMIT_TOP_N"] = "30"
os.environ["NIOME_GATK"] = str(ROOT / "tools/gatk/gatk")
os.environ["MINER_INSTANCE"] = "niome_dv_test"

from niome_subnet.cftr_miner_logic import process_cftr_task_for_miner

task_id = "df4120b1-4414-431c-99d6-ca3030d253d5"
task_path = ROOT / f"tasks/{task_id}/task.json"
task_data = json.loads(task_path.read_text())

read1_local = ROOT / f"reads/{task_id}/reads_1.fq"
read2_local = ROOT / f"reads/{task_id}/reads_2.fq"
task_data["input"]["read1_fastq"] = f"file://{read1_local}"
task_data["input"]["read2_fastq"] = f"file://{read2_local}"

print(f"=== Running task {task_id} with DeepVariant ENABLED ===", flush=True)
print(f"  read1: {task_data['input']['read1_fastq']}", flush=True)
print(f"  read2: {task_data['input']['read2_fastq']}", flush=True)
print(f"  region: {task_data['genome_context']['region']}", flush=True)

result = process_cftr_task_for_miner(
    task_data,
    base_dir=str(ROOT),
    logger=lambda msg: print(f"[miner] {msg}", flush=True),
)

paths = result.get("paths", {})
synapse_vcf = paths.get("synapse_vcf")
print(f"\n=== Result ===", flush=True)
print(f"  elapsed: {result.get('elapsed_time'):.2f}s", flush=True)
print(f"  synapse_vcf: {synapse_vcf}", flush=True)
print(f"  annotations: {len(result.get('cftr_annotations') or {})}", flush=True)

variant_count = sum(
    1 for line in result["vcf_content"].splitlines() if not line.startswith("#")
)
print(f"  variant rows: {variant_count}", flush=True)
