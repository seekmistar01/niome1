#!/usr/bin/env python3
"""Build real (FASTQ, truth) training samples from Niome leaderboard.

Two sources of "truth":

1. OFFICIAL truth (when leaderboard publishes it):
   GET /api/tasks/ground_truth_urls -> {"reads_1": <url>, "reads_2": <url>,
   "truth": <url>, "annotations": <url>}. Downloads via
   /api/tasks/ground_truth_file?url=<url>.

2. CONSENSUS proxy truth (always available for scored tasks):
   GET /api/miner_scores/task/<task_id> returns all miner submissions with
   final_score / precision / recall. We take variants that ≥ AGREEMENT_FRAC of
   top miners agree on. If at least one miner scored >= TOP_SCORE_THRESHOLD,
   their VCF is high-fidelity proxy truth.

Outputs into <out_dir>/<task_id>/:
   reads_1.fq.gz, reads_2.fq.gz   (if available)
   truth.vcf                      (official OR consensus proxy)
   annotations.json               (if official)
   metadata.json                  (source, top_score, n_submissions, hap stats)

Usage:
    ./venv/bin/python scripts/fetch_niome_truth.py \\
        --out /root/55miner/niome1/realdata \\
        --task-id df4120b1-4414-431c-99d6-ca3030d253d5 \\
        --source consensus

    # auto-detect tasks for which we have local FASTQs
    ./venv/bin/python scripts/fetch_niome_truth.py \\
        --out /root/55miner/niome1/realdata --auto-local
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

NIOME_ROOT = Path("/root/55miner/niome1")
LB_BASE = "https://niome-leaderboard.genomes.io"
HEADERS = {
    "Origin": LB_BASE,
    "Referer": f"{LB_BASE}/",
    "Accept": "application/json",
    "User-Agent": "niome-miner-trainer/1.0",
}

# Score thresholds
TOP_SCORE_THRESHOLD = 0.85    # any miner above this counts as "top"
AGREEMENT_FRAC = 0.66         # variant kept if ≥ this fraction of top miners agree


def http_get_json(url: str) -> Optional[dict]:
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as exc:
        # 404 with body might still carry an error message
        try:
            return json.loads(exc.read())
        except Exception:
            return {"error": str(exc)}
    except Exception as exc:
        return {"error": str(exc)}


def http_get_bytes(url: str) -> Optional[bytes]:
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            return r.read()
    except Exception as exc:
        print(f"  download error for {url[:120]}: {exc}", file=sys.stderr)
        return None


def fetch_all_miner_scores(task_id: str) -> List[dict]:
    url = f"{LB_BASE}/api/miner_scores/task/{urllib.parse.quote(task_id)}"
    data = http_get_json(url)
    if not data or "error" in data:
        print(f"  no scores for {task_id}: {data}", file=sys.stderr)
        return []
    items = data.get("items", data) if isinstance(data, dict) else data
    return items if isinstance(items, list) else []


def parse_vcf_from_log(log: str) -> List[Tuple[str, int, str, str, str]]:
    """Extract (chrom, pos, ref, alt, gt) tuples from a miner's log field."""
    out = []
    in_body = False
    for line in log.splitlines():
        if line.startswith("#CHROM"):
            in_body = True
            continue
        if not in_body or not line.strip():
            continue
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        chrom, pos_s, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        try:
            pos = int(pos_s)
        except ValueError:
            continue
        # GT is the first field of the SAMPLE column, split by ':'
        gt = parts[9].split(":")[0]
        # Normalize phased pipe to slash and sort for phase-insensitive compare
        gt_norm = "/".join(sorted(gt.replace("|", "/").split("/")))
        out.append((chrom, pos, ref, alt, gt_norm))
    return out


def build_consensus_truth(scores: List[dict]) -> Tuple[List[tuple], dict]:
    """Return (truth_records, stats). truth_records = sorted list of
    (chrom, pos, ref, alt, gt). Uses top miners (final_score >= TOP_SCORE_THRESHOLD)
    and keeps variants supported by ≥ AGREEMENT_FRAC.
    """
    top = [s for s in scores if (s.get("final_score") or 0) >= TOP_SCORE_THRESHOLD]
    if not top:
        # Fall back: take top 10 by final_score
        top = sorted(scores, key=lambda s: -(s.get("final_score") or 0))[:10]

    # Count (chrom, pos, ref, alt) occurrences and GT votes
    presence: Counter = Counter()
    gt_votes: Dict[Tuple[str, int, str, str], Counter] = defaultdict(Counter)
    for s in top:
        log = s.get("log") or ""
        seen_in_this_submission = set()
        for chrom, pos, ref, alt, gt in parse_vcf_from_log(log):
            key = (chrom, pos, ref, alt)
            if key in seen_in_this_submission:
                continue
            seen_in_this_submission.add(key)
            presence[key] += 1
            gt_votes[key][gt] += 1

    n_top = len(top)
    threshold = max(2, int(n_top * AGREEMENT_FRAC))
    truth = []
    for key, count in presence.items():
        if count < threshold:
            continue
        gt_winner, _ = gt_votes[key].most_common(1)[0]
        chrom, pos, ref, alt = key
        truth.append((chrom, pos, ref, alt, gt_winner))
    truth.sort(key=lambda r: (r[0], r[1]))

    top_score = max((s.get("final_score") or 0) for s in top) if top else 0.0
    stats = {
        "n_submissions": len(scores),
        "n_top_submissions": n_top,
        "agreement_threshold": threshold,
        "agreement_frac": AGREEMENT_FRAC,
        "top_score_threshold": TOP_SCORE_THRESHOLD,
        "top_score_observed": round(top_score, 4),
        "consensus_variant_count": len(truth),
    }
    return truth, stats


def write_truth_vcf(path: Path, records: List[tuple], task_id: str, source: str) -> None:
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##source=niome_{source}_truth\n")
        f.write(f"##task_id={task_id}\n")
        f.write("##contig=<ID=chr7,length=159345973>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for chrom, pos, ref, alt, gt in records:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def try_fetch_official_truth() -> Optional[dict]:
    """Try /api/tasks/ground_truth_urls. Returns dict with reads/truth/annotations
    URLs, or None if not yet published."""
    data = http_get_json(f"{LB_BASE}/api/tasks/ground_truth_urls")
    if not data or "error" in data:
        return None
    # Expected keys: reads_1, reads_2, truth, annotations
    required = {"reads_1", "reads_2", "truth"}
    if required.issubset(data.keys()):
        return data
    return None


_TASK_ID_RX = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
    re.IGNORECASE,
)


def extract_task_id_from_urls(urls: dict) -> Optional[str]:
    """Pull a uuid-style task_id from any of the truth URLs (S3 path)."""
    for key in ("truth", "annotations", "reads_1", "reads_2"):
        val = urls.get(key) or ""
        if isinstance(val, dict):
            val = val.get("url") or val.get("path") or ""
        m = _TASK_ID_RX.search(val)
        if m:
            return m.group(0).lower()
    return None


def proxy_via_ground_truth_file(remote_url: str) -> Optional[bytes]:
    """Use the leaderboard's CORS proxy to download a file that the JS would."""
    proxied = (
        f"{LB_BASE}/api/tasks/ground_truth_file?"
        f"url={urllib.parse.quote(remote_url, safe='')}"
    )
    return http_get_bytes(proxied)


def stage_official_files(task_dir: Path, urls: dict) -> dict:
    out = {}
    for key, fname in [
        ("reads_1", "reads_1.fq.gz"),
        ("reads_2", "reads_2.fq.gz"),
        ("truth", "truth.vcf"),
        ("annotations", "annotations.json"),
    ]:
        url = urls.get(key)
        if not url:
            continue
        body = proxy_via_ground_truth_file(url)
        if body is None:
            continue
        target = task_dir / fname
        target.write_bytes(body)
        out[key] = str(target)
    return out


def process_one_task(
    task_id: str,
    out_root: Path,
    *,
    use_official: bool = True,
    use_consensus: bool = True,
    local_reads_root: Optional[Path] = None,
) -> dict:
    task_dir = out_root / task_id
    task_dir.mkdir(parents=True, exist_ok=True)

    result = {"task_id": task_id, "sources": [], "files": {}}

    # 1. Try official truth (only if it's currently published)
    if use_official:
        urls = try_fetch_official_truth()
        if urls:
            files = stage_official_files(task_dir, urls)
            if files.get("truth"):
                result["sources"].append("official")
                result["files"].update(files)

    # 2. Always also build consensus proxy truth from miner_scores
    if use_consensus and "truth" not in result["files"]:
        scores = fetch_all_miner_scores(task_id)
        if scores:
            truth, stats = build_consensus_truth(scores)
            if truth:
                truth_path = task_dir / "truth.vcf"
                write_truth_vcf(truth_path, truth, task_id, "consensus")
                result["sources"].append("consensus")
                result["files"]["truth"] = str(truth_path)
                result["consensus_stats"] = stats

    # 3. Stage local FASTQs if available and not already present
    if local_reads_root and "reads_1" not in result["files"]:
        local_dir = local_reads_root / task_id
        for src_name, dst_name in [("reads_1.fq", "reads_1.fq.gz"),
                                    ("reads_2.fq", "reads_2.fq.gz")]:
            src = local_dir / src_name
            dst = task_dir / dst_name
            if not src.exists():
                continue
            with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            result["files"][src_name.replace(".fq", "")] = str(dst)
            result["sources"].append("local_fastq")

    # 4. Write metadata
    (task_dir / "metadata.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--task-id", action="append", default=[],
                        help="specific task_id(s) to fetch; repeatable")
    parser.add_argument("--auto-local", action="store_true",
                        help="fetch all tasks for which we have local FASTQs under "
                             "niome1/reads/")
    parser.add_argument("--official-latest", action="store_true",
                        help="poll the leaderboard's currently-available official "
                             "truth (ground_truth_urls) and download all 4 files. "
                             "No-op if nothing new is published. Designed for cron.")
    parser.add_argument("--source", choices=["both", "official", "consensus"],
                        default="both")
    parser.add_argument("--local-reads",
                        default=str(NIOME_ROOT / "reads"),
                        help="root dir containing local FASTQ per task")
    args = parser.parse_args()

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    local_reads = Path(args.local_reads) if args.local_reads else None

    task_ids = list(args.task_id)
    if args.auto_local and local_reads:
        for d in sorted(local_reads.iterdir()):
            if d.is_dir() and (d / "reads_1.fq").exists():
                task_ids.append(d.name)

    # --official-latest mode: poll endpoint, derive task_id from URL, dedupe
    if args.official_latest:
        urls = try_fetch_official_truth()
        if not urls:
            print("[official-latest] not available yet", flush=True)
            return
        tid = extract_task_id_from_urls(urls)
        if not tid:
            tid = "unknown_" + time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
            print(f"[official-latest] could not extract task_id from URLs, "
                  f"using {tid}", flush=True)
        target_dir = out_root / tid
        if (target_dir / "truth.vcf").exists():
            print(f"[official-latest] {tid} already downloaded, skipping",
                  flush=True)
            return
        target_dir.mkdir(parents=True, exist_ok=True)
        print(f"[official-latest] downloading {tid}", flush=True)
        files = stage_official_files(target_dir, urls)
        meta = {
            "task_id": tid,
            "source": "official",
            "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "files": files,
            "url_sources": urls,
        }
        (target_dir / "metadata.json").write_text(json.dumps(meta, indent=2))
        print(f"[official-latest] saved {len(files)} files to {target_dir}",
              flush=True)
        return

    if not task_ids:
        sys.exit("ERROR: no task_ids specified (use --task-id, --auto-local, "
                 "or --official-latest)")

    use_official = args.source in ("both", "official")
    use_consensus = args.source in ("both", "consensus")

    print(f"Processing {len(task_ids)} task(s) -> {out_root}", flush=True)
    results = []
    for tid in task_ids:
        print(f"  {tid}", flush=True)
        r = process_one_task(
            tid, out_root,
            use_official=use_official,
            use_consensus=use_consensus,
            local_reads_root=local_reads,
        )
        results.append(r)
        srcs = ",".join(r["sources"]) or "none"
        stats = r.get("consensus_stats", {})
        msg = f"    sources={srcs}"
        if stats:
            msg += (
                f"  top_score={stats['top_score_observed']}"
                f"  n_top={stats['n_top_submissions']}"
                f"  variants={stats['consensus_variant_count']}"
            )
        print(msg, flush=True)

    (out_root / "_index.json").write_text(json.dumps({
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "results": results,
    }, indent=2))
    print(f"\nDone. Index written to {out_root}/_index.json", flush=True)


if __name__ == "__main__":
    main()
