"""Shared VCF compression and bcftools norm helpers (validator + miner)."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional, Union


def ensure_reference_index(reference_fasta: Union[str, Path], *, require_bwa: bool = False) -> Path:
    """Create samtools .fai for bcftools norm; optionally bwa index for alignment."""
    ref_path = Path(reference_fasta).resolve()
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_path}")

    fai_path = Path(f"{ref_path}.fai")
    if not fai_path.exists():
        _run(
            ["samtools", "faidx", str(ref_path)],
            f"index reference with samtools faidx ({ref_path})",
        )

    if require_bwa and not Path(f"{ref_path}.bwt").exists():
        _run(
            ["bwa", "index", str(ref_path)],
            f"index reference with bwa ({ref_path})",
        )
    return ref_path


def preprocess_vcf(vcf_path: Union[str, Path]) -> Path:
    """Bgzip and tabix a plain VCF; pass through if already compressed."""
    source = Path(vcf_path).resolve()
    if str(source).endswith(".gz"):
        if not Path(f"{source}.tbi").exists():
            _run(
                ["tabix", "-f", "-p", "vcf", str(source)],
                f"index VCF {source}",
            )
        return source

    output = Path(f"{source}.gz")
    if output.exists():
        output.unlink()
    if Path(f"{output}.tbi").exists():
        Path(f"{output}.tbi").unlink()

    with source.open("rb") as handle_in, output.open("wb") as handle_out:
        subprocess.run(
            ["bgzip", "-c"],
            stdin=handle_in,
            stdout=handle_out,
            check=True,
        )
    _run(["tabix", "-f", "-p", "vcf", str(output)], f"index VCF {output}")
    return output


def normalize_vcf(
    vcf_in: Union[str, Path],
    reference_fasta: Union[str, Path],
    output_vcf: Union[str, Path],
) -> Path:
    """Left-align and normalize alleles against reference (validator uses -c x)."""
    ref_path = ensure_reference_index(reference_fasta)
    vcf_path = Path(vcf_in).resolve()
    out_path = Path(output_vcf).resolve()

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    if str(vcf_path).endswith(".vcf") and not str(vcf_path).endswith(".vcf.gz"):
        vcf_path = preprocess_vcf(vcf_path)
    elif not Path(f"{vcf_path}.tbi").exists():
        _run(["tabix", "-f", "-p", "vcf", str(vcf_path)], f"index VCF {vcf_path}")

    if out_path.exists():
        out_path.unlink()
    tbi = Path(f"{out_path}.tbi")
    if tbi.exists():
        tbi.unlink()

    _run(
        [
            "bcftools",
            "norm",
            "-f",
            str(ref_path),
            "-c",
            "x",
            "-m",
            "-both",
            "-Oz",
            "-o",
            str(out_path),
            str(vcf_path),
        ],
        f"bcftools norm {vcf_path}",
    )
    _run(["bcftools", "index", "-f", str(out_path)], f"index normalized VCF {out_path}")
    return out_path


def _run(command: list[str], description: str) -> None:
    try:
        subprocess.run(command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        stdout = (exc.stdout or "").strip()
        details = stderr or stdout or f"exit code {exc.returncode}"
        raise RuntimeError(f"Failed to {description}: {details}") from exc


def verify_vcf_for_validator_scoring(
    vcf_text: str,
    reference_fasta: Union[str, Path],
    *,
    work_dir: Optional[Union[str, Path]] = None,
) -> None:
    """Raise RuntimeError if validator scoring bcftools norm -c x would fail on this VCF."""
    if not (vcf_text or "").strip():
        return
    ref_path = ensure_reference_index(reference_fasta)
    base = Path(work_dir or Path.cwd() / "work" / "validator_vcf_check")
    base.mkdir(parents=True, exist_ok=True)
    draft = base / "miner.check.vcf"
    norm_out = base / "miner.check.norm.vcf.gz"
    draft.write_text(vcf_text, encoding="utf-8")
    if norm_out.exists():
        norm_out.unlink()
    tbi = Path(f"{norm_out}.tbi")
    if tbi.exists():
        tbi.unlink()
    try:
        gz_path = preprocess_vcf(draft)
        normalize_vcf(gz_path, ref_path, norm_out)
    except Exception as exc:
        raise RuntimeError(
            f"Miner VCF failed bcftools norm (validator would score 0): {exc}"
        ) from exc
