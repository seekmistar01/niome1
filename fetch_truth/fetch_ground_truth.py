#!/usr/bin/env python3
"""Fetch NIOME task + ground truth and save under fetch_truth/<task_id>/.

Usage (from subnet-niome root, with venv active and `pip install -e .`):

  python fetch_truth/fetch_ground_truth.py \\
    --wallet.name seekmistar3 --wallet.hotkey seekmistar01 --netuid 55

  python fetch_truth/fetch_ground_truth.py --fetch truth \\
    --task-id 768bdcb5-4269-4c72-8705-d77d729b3db8
"""

from __future__ import annotations

import argparse
import asyncio
import json
import re
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import aiohttp

from niome_subnet.genomics.model import GroundTruth, Task
from niome_subnet.utils.constants import (
    BASE_DELAY_SECONDS,
    GROUND_TRUTH_URL,
    MAINNET_UID,
    MAX_TASK_RETRIES,
    TASK_REQUEST_TIMEOUT,
    TASK_URL,
    TESTNET_UID,
)

SCRIPT_DIR = Path(__file__).resolve().parent


def build_signature_headers(
    signature: str, hotkey: str, timestamp: str, netuid: str
) -> dict[str, str]:
    return {
        "X-Signature": signature,
        "X-Hotkey": hotkey,
        "X-Netuid": netuid,
        "X-Timestamp": timestamp,
    }


def sign_request(wallet, netuid: int) -> dict[str, str]:
    timestamp = str(time.time())
    canonical = json.dumps(
        {
            "payload": "{}",
            "hotkey": wallet.hotkey.ss58_address,
            "netuid": str(netuid),
            "timestamp": timestamp,
        },
        separators=(",", ":"),
        sort_keys=True,
    )
    signature = wallet.hotkey.sign(canonical).hex()
    return build_signature_headers(
        signature=signature,
        hotkey=wallet.hotkey.ss58_address,
        timestamp=timestamp,
        netuid=str(netuid),
    )


def safe_task_id(task_id: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", task_id.strip())
    return cleaned[:120] or "unknown-task"


def fetch_json_url(url: str) -> dict[str, Any]:
    with urllib.request.urlopen(url, timeout=TASK_REQUEST_TIMEOUT) as resp:
        return json.loads(resp.read())


async def signed_post(wallet, netuid: int, api_url: str) -> dict[str, Any]:
    last_error: Exception | None = None

    for attempt in range(1, MAX_TASK_RETRIES + 1):
        headers = sign_request(wallet, netuid)
        try:
            async with aiohttp.ClientSession() as client:
                async with client.post(
                    api_url,
                    headers=headers,
                    json={},
                    timeout=aiohttp.ClientTimeout(total=TASK_REQUEST_TIMEOUT),
                ) as response:
                    body = await response.text()
                    if response.status != 201:
                        detail = body[:500]
                        try:
                            detail = json.loads(body).get("error", detail)
                        except json.JSONDecodeError:
                            pass
                        raise RuntimeError(
                            f"API returned {response.status}: {detail}"
                        )
                    return json.loads(body)
        except Exception as exc:
            last_error = exc
            if attempt < MAX_TASK_RETRIES:
                delay = BASE_DELAY_SECONDS * (2 ** (attempt - 1))
                print(
                    f"POST {api_url} attempt {attempt} failed ({exc}); "
                    f"retrying in {delay}s...",
                    file=sys.stderr,
                )
                await asyncio.sleep(delay)
            else:
                break

    raise RuntimeError(f"All {MAX_TASK_RETRIES} attempts failed") from last_error


async def fetch_task_from_api(wallet, netuid: int) -> tuple[Task, dict[str, Any]]:
    envelope = await signed_post(wallet, netuid, TASK_URL)
    task_url = envelope.get("task_url", "")
    if not task_url:
        raise RuntimeError(f"Missing task_url in API response: {envelope}")
    task_data = await asyncio.to_thread(fetch_json_url, task_url)
    return Task(**task_data), envelope


async def fetch_ground_truth_from_api(
    wallet, netuid: int
) -> tuple[GroundTruth, dict[str, Any], str]:
    envelope = await signed_post(wallet, netuid, GROUND_TRUTH_URL)
    ground_truth_url = envelope.get("ground_truth_url", "")
    if not ground_truth_url:
        raise RuntimeError(f"Missing ground_truth_url in API response: {envelope}")
    raw = await asyncio.to_thread(fetch_json_url, ground_truth_url)
    return GroundTruth(**raw), {**envelope, **raw}, ground_truth_url


def resolve_task_id(
    *,
    explicit: str | None,
    task: Task | None,
    truth_envelope: dict[str, Any] | None,
) -> str:
    if explicit:
        return explicit
    if task is not None:
        return task.task_id
    if truth_envelope:
        for key in ("task_id", "taskId", "id"):
            value = truth_envelope.get(key)
            if value:
                return str(value)
    raise RuntimeError(
        "Could not determine task_id. Pass --task-id or fetch task with --fetch task/both."
    )


def download_file(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


def save_task_artifacts(task: Task, task_dir: Path, *, download_files: bool) -> dict[str, Any]:
    record: dict[str, Any] = task.model_dump()
    task_dir.mkdir(parents=True, exist_ok=True)
    (task_dir / "task.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8"
    )

    if not download_files:
        return record

    read1 = download_file(task.input.read1_fastq, task_dir / "reads_1.fq")
    read2 = download_file(task.input.read2_fastq, task_dir / "reads_2.fq")
    record["local_paths"] = {
        "read1_fastq": str(read1),
        "read2_fastq": str(read2),
    }
    (task_dir / "task.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8"
    )
    return record


def save_ground_truth_artifacts(
    ground_truth: GroundTruth,
    task_dir: Path,
    *,
    download_files: bool,
    api_envelope: dict[str, Any],
    ground_truth_url: str,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "ground_truth_url": ground_truth_url,
        "api_envelope": api_envelope,
        **ground_truth.model_dump(),
    }
    task_dir.mkdir(parents=True, exist_ok=True)

    if not download_files:
        (task_dir / "ground_truth.json").write_text(
            json.dumps(record, indent=2), encoding="utf-8"
        )
        return record

    truth_path = download_file(ground_truth.truth_vcf, task_dir / "truth.vcf")
    ref_path = download_file(ground_truth.ref, task_dir / "ref.fa")
    ann_path = download_file(
        ground_truth.cftr2_annotations, task_dir / "cftr2_annotations.json"
    )
    record["local_paths"] = {
        "truth_vcf": str(truth_path),
        "ref": str(ref_path),
        "cftr2_annotations": str(ann_path),
    }
    (task_dir / "ground_truth.json").write_text(
        json.dumps(record, indent=2), encoding="utf-8"
    )
    return record


async def run(args: argparse.Namespace) -> int:
    import bittensor as bt

    wallet = bt.Wallet(name=args.wallet_name, hotkey=args.wallet_hotkey)
    netuid = args.netuid
    fetch_mode = args.fetch

    print(f"Wallet: {wallet.name}/{wallet.hotkey_str}")
    print(f"Hotkey: {wallet.hotkey.ss58_address}")
    print(f"Netuid: {netuid}")
    print(f"Fetch: {fetch_mode}")

    task: Task | None = None
    task_envelope: dict[str, Any] | None = None
    ground_truth: GroundTruth | None = None
    truth_envelope: dict[str, Any] | None = None
    ground_truth_url = ""

    if fetch_mode in ("task", "both"):
        print(f"Task API: {TASK_URL}")
        task, task_envelope = await fetch_task_from_api(wallet, netuid)
        print(f"Fetched task: {task.task_id}")

    if fetch_mode in ("truth", "both"):
        print(f"Ground truth API: {GROUND_TRUTH_URL}")
        ground_truth, truth_envelope, ground_truth_url = await fetch_ground_truth_from_api(
            wallet, netuid
        )
        print(f"Fetched ground truth metadata from {ground_truth_url}")

    task_id = resolve_task_id(
        explicit=args.task_id,
        task=task,
        truth_envelope=truth_envelope,
    )
    task_dir = Path(args.output_dir) / safe_task_id(task_id)
    print(f"Saving to: {task_dir.resolve()}")

    if task is not None:
        save_task_artifacts(
            task, task_dir, download_files=not args.metadata_only
        )
        print(f"Wrote task.json under {task_dir}")

    if ground_truth is not None and truth_envelope is not None:
        save_ground_truth_artifacts(
            ground_truth,
            task_dir,
            download_files=not args.metadata_only,
            api_envelope=truth_envelope,
            ground_truth_url=ground_truth_url,
        )
        print(f"Wrote ground_truth.json under {task_dir}")

    if task_envelope is not None:
        (task_dir / "api_task_envelope.json").write_text(
            json.dumps(task_envelope, indent=2), encoding="utf-8"
        )

    print("Done.")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch NIOME task/ground truth into fetch_truth/<task_id>/."
    )
    parser.add_argument(
        "--wallet.name",
        dest="wallet_name",
        default="default",
        help="Bittensor wallet name",
    )
    parser.add_argument(
        "--wallet.hotkey",
        dest="wallet_hotkey",
        default="default",
        help="Bittensor wallet hotkey",
    )
    parser.add_argument(
        "--netuid",
        type=int,
        default=MAINNET_UID,
        help=f"Subnet netuid (default: {MAINNET_UID} mainnet, {TESTNET_UID} testnet)",
    )
    parser.add_argument(
        "--fetch",
        choices=("task", "truth", "both"),
        default="both",
        help="What to fetch from the API (default: both)",
    )
    parser.add_argument(
        "--task-id",
        default=None,
        help="Override task folder name when API does not return task_id (truth-only)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR,
        help=f"Base directory for per-task folders (default: {SCRIPT_DIR})",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Only save JSON manifests, skip downloading reads/VCF/ref",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        raise SystemExit(asyncio.run(run(args)))  # bt imported inside run()
    except KeyboardInterrupt:
        raise SystemExit(130) from None
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
