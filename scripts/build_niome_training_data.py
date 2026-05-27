#!/usr/bin/env python3
"""Generate Niome-like synthetic training data for DeepVariant fine-tuning.

Each sample is a per-task directory:
    <out>/sample_XXXXX/
        reads_1.fq.gz, reads_2.fq.gz   # paired-end Illumina-like reads
        truth.vcf                       # ground-truth variants with assigned GT
        sample.json                     # metadata (seed, hap_frac, depth, variant_count)

Calibrated to real Niome measurements taken from the user's archived tasks:
    - 75bp paired-end
    - mean coverage ~16x across chr7:117480000-117670000
    - 15-85% allele imbalance between haplotypes
    - variants sourced from local cftr_clinvar_panel.tsv

Requires: bcftools, samtools, wgsim (samtools-bundled). All present at
/usr/bin/. Output dir is created if missing.

Example:
    ./venv/bin/python scripts/build_niome_training_data.py \\
        --out /root/55miner/traindatasets_first \\
        --count 100 --workers 12 --seed 42
"""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import os
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import List, Tuple

NIOME_ROOT = Path("/root/55miner/niome1")
REGION_CHROM = "chr7"
REGION_START = 117480000
REGION_END = 117670000
REGION = f"{REGION_CHROM}:{REGION_START}-{REGION_END}"
REF_FA = NIOME_ROOT / "data/chr7.fa"
CLINVAR_PANEL = NIOME_ROOT / "panel/cftr_clinvar_panel.tsv"
SUPPLEMENTAL_PANEL = NIOME_ROOT / "panel/cftr_supplemental_panel.tsv"

# Calibrated from real Niome measurements
READ_LENGTH = 75
INSERT_SIZE_MEAN = 300
INSERT_SIZE_SD = 30
TARGET_TOTAL_PAIRS = 20000        # ~16x coverage on 190kb region
ERROR_RATE = 0.005                # Illumina ~0.5%
CONTAM_FRACTION = 0.05            # ~5% off-region reads
# Calibrated 2026-05-27 from 6 real tasks' consensus truth (top-miner agreement
# on tasks with ≥0.85 score). Real Niome truth has 19-31 variants per task with
# ~33% hom-alt and ~39% indels.
VARIANTS_PER_SAMPLE_MIN = 20
VARIANTS_PER_SAMPLE_MAX = 32

# GT assignment distribution
HOM_FRAC = 0.33
HET_FRAC = 0.67

# Indel oversampling weight: ClinVar pool is ~15% indels, but real Niome truth
# is ~39% indels. Weight each indel candidate by this factor at selection time.
INDEL_WEIGHT = 3.0

# Off-region for contamination — same chromosome, non-CFTR window
CONTAM_REGION = "chr7:60000000-60010000"


def run(cmd, cwd=None, check=True, capture=True):
    """Run a subprocess with friendly error messages."""
    if isinstance(cmd, str):
        return subprocess.run(
            cmd, shell=True, cwd=cwd, check=check,
            capture_output=capture, text=True,
        )
    return subprocess.run(
        cmd, cwd=cwd, check=check, capture_output=capture, text=True,
    )


def load_panel_variants() -> List[Tuple[int, str, str]]:
    """Return [(pos, ref, alt)] for in-region panel variants ≤ 50bp."""
    rows = []
    for panel in (CLINVAR_PANEL, SUPPLEMENTAL_PANEL):
        if not panel.exists():
            continue
        with panel.open() as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                if fields[0].lower() in {"variation_id", "variationid", "id"}:
                    continue
                try:
                    pos = int(fields[2])
                except ValueError:
                    continue
                if pos < REGION_START or pos > REGION_END:
                    continue
                ref, alt = fields[3], fields[4]
                if not ref or not alt or len(ref) > 50 or len(alt) > 50:
                    continue
                if any(c not in "ACGT" for c in ref + alt):
                    continue
                rows.append((pos, ref, alt))
    # Deduplicate
    return sorted(set(rows))


def write_vcf(path: Path, records: List[Tuple[int, str, str, str]]) -> None:
    """Write a minimal valid VCF with the supplied (pos, ref, alt, gt) records."""
    with path.open("w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={REGION_CHROM},length=159345973>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for pos, ref, alt, gt in sorted(records, key=lambda r: r[0]):
            f.write(
                f"{REGION_CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
            )


def pick_sample_variants(
    panel: List[Tuple[int, str, str]],
    rng: random.Random,
) -> List[Tuple[int, str, str, str]]:
    """Sample variants and assign per-haplotype genotype labels.

    Returns truth records as (pos, ref, alt, gt) and ensures non-overlapping
    positions (variants spaced ≥ 50 bp apart so haplotype construction is safe).
    """
    n_target = rng.randint(VARIANTS_PER_SAMPLE_MIN, VARIANTS_PER_SAMPLE_MAX)
    # Weighted sampling without replacement (Efraimidis-Spirakis A-Res):
    # key_i = u_i ** (1 / w_i), take top-N by key. Larger w → larger key in
    # expectation → more likely to be picked.
    weights = [
        INDEL_WEIGHT if (len(ref) != 1 or len(alt) != 1) else 1.0
        for pos, ref, alt in panel
    ]
    keyed = [
        (rng.random() ** (1.0 / max(w, 1e-9)), idx)
        for idx, w in enumerate(weights)
    ]
    keyed.sort(reverse=True)  # largest keys first
    shuffled_indices = [idx for _, idx in keyed]

    selected = []
    used_positions = []  # sorted
    for idx in shuffled_indices:
        pos, ref, alt = panel[idx]
        # Spacing check: keep at least 50bp gap from previously chosen variant
        if any(abs(pos - up) < 50 for up in used_positions):
            continue
        gt = "1/1" if rng.random() < HOM_FRAC else "0/1"
        selected.append((pos, ref, alt, gt))
        used_positions.append(pos)
        if len(selected) >= n_target:
            break

    return sorted(selected, key=lambda r: r[0])


def variants_to_haplotypes(
    truth: List[Tuple[int, str, str, str]],
    rng: random.Random,
) -> Tuple[List[Tuple[int, str, str]], List[Tuple[int, str, str]]]:
    """Split truth into hap1 / hap2 variant lists.

    1/1 → present on both. 0/1 → randomly assigned to one of two haps.
    """
    hap1, hap2 = [], []
    for pos, ref, alt, gt in truth:
        if gt == "1/1":
            hap1.append((pos, ref, alt))
            hap2.append((pos, ref, alt))
        else:  # 0/1
            if rng.random() < 0.5:
                hap1.append((pos, ref, alt))
            else:
                hap2.append((pos, ref, alt))
    return hap1, hap2


def make_haplotype_fasta(
    hap_vcf: Path,
    out_fasta: Path,
    workdir: Path,
) -> None:
    """Apply variants to reference region → haplotype FASTA."""
    # bcftools needs bgzipped+tabixed VCF
    run(f"bgzip -f -c {hap_vcf} > {hap_vcf}.gz", capture=True)
    run(["tabix", "-f", "-p", "vcf", f"{hap_vcf}.gz"])
    # Pull region from reference, then apply consensus
    run(
        f"samtools faidx {REF_FA} {REGION} | "
        f"bcftools consensus {hap_vcf}.gz > {out_fasta}",
        capture=True,
    )
    # Force contig name back to chr7 (samtools faidx names it chr7:start-end)
    # wgsim/aligners will then map to chr7 positions correctly when we rewrite
    rename_contig(out_fasta)


def rename_contig(fasta: Path) -> None:
    """Replace the region-named header with chr7 so reads align to genome coords."""
    with fasta.open() as f:
        text = f.read()
    new_text = text.replace(f">{REGION}", f">{REGION_CHROM}")
    with fasta.open("w") as f:
        f.write(new_text)


def shift_wgsim_reads(
    raw_fq: Path,
    out_fq: Path,
    region_start: int,
) -> None:
    """wgsim writes positions relative to the haplotype FASTA (which starts at 1).
    We need positions relative to the full chromosome. We rename the contig header
    in the FASTA (above) so reads carry the chr7 name; the position offset is
    embedded in the FASTQ read name from wgsim, e.g. chr7_4321_4895_0:0:0_... We
    rewrite the embedded position to be region_start + pos - 1 so downstream
    alignment is to the right chromosome coords.

    Actually, this is unnecessary if the haplotype FASTA is a subset starting at
    region_start. bwa mem will simply align reads against the full chr7.fa and
    discover the correct coordinates regardless of the read name. So this
    function is a no-op pass-through; we keep it for future flexibility.
    """
    if raw_fq == out_fq:
        return
    shutil.move(str(raw_fq), str(out_fq))


def simulate_hap_reads(
    hap_fasta: Path,
    n_pairs: int,
    out_r1: Path,
    out_r2: Path,
    seed: int,
) -> None:
    """Run wgsim on a haplotype FASTA to generate paired-end reads."""
    cmd = [
        "wgsim",
        "-N", str(int(n_pairs)),
        "-1", str(READ_LENGTH),
        "-2", str(READ_LENGTH),
        "-d", str(INSERT_SIZE_MEAN),
        "-s", str(INSERT_SIZE_SD),
        "-e", str(ERROR_RATE),
        "-r", "0",          # don't add extra mutations; haplotype already has them
        "-R", "0",          # indel rate within reads
        "-X", "0",          # extension rate
        "-S", str(seed),
        str(hap_fasta),
        str(out_r1),
        str(out_r2),
    ]
    run(cmd, capture=True)


def simulate_contamination(
    n_pairs: int,
    out_r1: Path,
    out_r2: Path,
    seed: int,
    workdir: Path,
) -> None:
    """Generate ~5% reads from a non-CFTR region (light contamination)."""
    contam_fa = workdir / "contam.fa"
    run(
        f"samtools faidx {REF_FA} {CONTAM_REGION} > {contam_fa}",
        capture=True,
    )
    rename_contig(contam_fa)
    simulate_hap_reads(contam_fa, n_pairs, out_r1, out_r2, seed)


def append_fastq(src: Path, dst_fp) -> None:
    """Append plain FASTQ src to an already-open gzipped destination file."""
    with src.open("rb") as f:
        shutil.copyfileobj(f, dst_fp)


def build_one_sample(args):
    """Worker function: generate one full sample directory."""
    idx, out_root, seed_base, panel = args
    seed = seed_base + idx
    rng = random.Random(seed)

    sample_id = f"sample_{idx:05d}"
    sample_dir = out_root / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    workdir = Path(tempfile.mkdtemp(prefix=f"niome_sim_{idx:05d}_"))
    try:
        # 1. Pick truth variants and split into haplotypes
        truth = pick_sample_variants(panel, rng)
        hap1_vars, hap2_vars = variants_to_haplotypes(truth, rng)

        # 2. Write per-hap VCF (homozygous 1/1 so bcftools consensus applies)
        hap1_vcf = workdir / "hap1.vcf"
        hap2_vcf = workdir / "hap2.vcf"
        write_vcf(hap1_vcf, [(p, r, a, "1/1") for p, r, a in hap1_vars])
        write_vcf(hap2_vcf, [(p, r, a, "1/1") for p, r, a in hap2_vars])

        # 3. Build haplotype FASTAs
        hap1_fa = workdir / "hap1.fa"
        hap2_fa = workdir / "hap2.fa"
        make_haplotype_fasta(hap1_vcf, hap1_fa, workdir)
        make_haplotype_fasta(hap2_vcf, hap2_fa, workdir)

        # 4. Pick haplotype fraction and total coverage
        # Beta(2,2) puts mass near 0.5 but allows tails; clamp to [0.15, 0.85]
        f = max(0.15, min(0.85, rng.betavariate(2, 2)))
        # Total pairs: small jitter around target
        total_pairs = int(TARGET_TOTAL_PAIRS * rng.uniform(0.9, 1.1))
        contam_pairs = int(total_pairs * CONTAM_FRACTION)
        primary_pairs = total_pairs - contam_pairs
        hap1_pairs = int(primary_pairs * f)
        hap2_pairs = primary_pairs - hap1_pairs

        # 5. Simulate reads from each hap
        h1_r1 = workdir / "h1_R1.fq"
        h1_r2 = workdir / "h1_R2.fq"
        h2_r1 = workdir / "h2_R1.fq"
        h2_r2 = workdir / "h2_R2.fq"
        c_r1 = workdir / "c_R1.fq"
        c_r2 = workdir / "c_R2.fq"
        simulate_hap_reads(hap1_fa, hap1_pairs, h1_r1, h1_r2, seed + 1000)
        simulate_hap_reads(hap2_fa, hap2_pairs, h2_r1, h2_r2, seed + 2000)
        simulate_contamination(contam_pairs, c_r1, c_r2, seed + 3000, workdir)

        # 6. Merge into final paired FASTQs (gzipped)
        out_r1 = sample_dir / "reads_1.fq.gz"
        out_r2 = sample_dir / "reads_2.fq.gz"
        with gzip.open(out_r1, "wb") as g1, gzip.open(out_r2, "wb") as g2:
            for src1, src2 in [(h1_r1, h1_r2), (h2_r1, h2_r2), (c_r1, c_r2)]:
                append_fastq(src1, g1)
                append_fastq(src2, g2)

        # 7. Truth VCF and metadata
        write_vcf(sample_dir / "truth.vcf", truth)
        meta = {
            "sample_id": sample_id,
            "seed": seed,
            "hap1_fraction": round(f, 3),
            "total_pairs": total_pairs,
            "hap1_pairs": hap1_pairs,
            "hap2_pairs": hap2_pairs,
            "contam_pairs": contam_pairs,
            "variant_count": len(truth),
            "hom_count": sum(1 for r in truth if r[3] == "1/1"),
            "het_count": sum(1 for r in truth if r[3] == "0/1"),
            "region": REGION,
            "read_length": READ_LENGTH,
            "insert_mean": INSERT_SIZE_MEAN,
            "error_rate": ERROR_RATE,
        }
        (sample_dir / "sample.json").write_text(json.dumps(meta, indent=2))
        return idx, True, None
    except subprocess.CalledProcessError as exc:
        return idx, False, f"{exc}\nSTDERR: {exc.stderr[:500]}"
    except Exception as exc:
        return idx, False, str(exc)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--count", type=int, default=100, help="number of samples")
    parser.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 2))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--start-idx", type=int, default=1,
                        help="first sample index (for appending to a batch)")
    args = parser.parse_args()

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    panel = load_panel_variants()
    if not panel:
        raise SystemExit("ERROR: no panel variants loaded — check panel/ files")
    print(f"Loaded {len(panel)} eligible panel variants from CFTR region", flush=True)

    print(
        f"Generating {args.count} samples to {out_root} "
        f"with {args.workers} workers (seed={args.seed})",
        flush=True,
    )

    manifest = {
        "count_requested": args.count,
        "seed": args.seed,
        "start_idx": args.start_idx,
        "settings": {
            "region": REGION,
            "read_length": READ_LENGTH,
            "insert_mean": INSERT_SIZE_MEAN,
            "insert_sd": INSERT_SIZE_SD,
            "target_total_pairs": TARGET_TOTAL_PAIRS,
            "error_rate": ERROR_RATE,
            "contam_fraction": CONTAM_FRACTION,
            "variants_per_sample_min": VARIANTS_PER_SAMPLE_MIN,
            "variants_per_sample_max": VARIANTS_PER_SAMPLE_MAX,
            "hom_frac": HOM_FRAC,
        },
        "panel_size": len(panel),
        "generation_start": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2))

    jobs = [
        (args.start_idx + i, out_root, args.seed, panel)
        for i in range(args.count)
    ]

    start = time.perf_counter()
    ok = 0
    failures = []
    with mp.Pool(processes=args.workers) as pool:
        for done, (idx, success, err) in enumerate(
            pool.imap_unordered(build_one_sample, jobs), 1
        ):
            if success:
                ok += 1
            else:
                failures.append((idx, err))
            if done % 10 == 0 or done == args.count:
                elapsed = time.perf_counter() - start
                rate = done / elapsed
                eta = (args.count - done) / rate if rate else float("inf")
                print(
                    f"  [{done}/{args.count}] ok={ok} fail={len(failures)} "
                    f"rate={rate:.1f}/s eta={eta:.0f}s",
                    flush=True,
                )

    manifest["generation_end"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["count_succeeded"] = ok
    manifest["count_failed"] = len(failures)
    manifest["failures"] = failures[:20]  # cap
    manifest["wall_seconds"] = round(time.perf_counter() - start, 1)
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(
        f"\nDone: {ok}/{args.count} samples in "
        f"{manifest['wall_seconds']}s. Output: {out_root}",
        flush=True,
    )
    if failures:
        print(f"FAILURES ({len(failures)}):", flush=True)
        for idx, err in failures[:5]:
            print(f"  sample_{idx:05d}: {err[:200]}", flush=True)


if __name__ == "__main__":
    main()
