"""Canonical subnet reference layout (validator expects data/ref.fa)."""

from __future__ import annotations

import os
from pathlib import Path


def ensure_canonical_reference(base_dir: Path) -> tuple[Path, str]:
    """
    Resolve the reference FASTA used for alignment, norm, and VCF headers.

    Returns (absolute_fasta_path, header_path_for_vcf) where header is always
    the relative path validators use: data/ref.fa
    """
    base_dir = base_dir.resolve()
    data_dir = base_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    canonical = data_dir / "ref.fa"
    header_path = "data/ref.fa"

    env_ref = os.environ.get("NIOME_CFTR_REF", "").strip()
    if env_ref:
        source = Path(env_ref)
        if not source.is_absolute():
            source = base_dir / source
        source = source.resolve()
    else:
        chr7 = data_dir / "chr7.fa"
        if chr7.exists() or chr7.is_symlink():
            source = chr7.resolve()
        elif canonical.exists() or canonical.is_symlink():
            source = canonical.resolve()
        else:
            source = chr7.resolve()

    if not source.exists():
        raise FileNotFoundError(
            f"Reference FASTA not found: {source}. "
            "Set NIOME_CFTR_REF or place data/chr7.fa / data/ref.fa."
        )

    if canonical.exists() or canonical.is_symlink():
        if canonical.resolve() != source:
            canonical.unlink()
            canonical.symlink_to(source)
    else:
        canonical.symlink_to(source)

    chr7 = data_dir / "chr7.fa"
    if not chr7.exists() and not chr7.is_symlink():
        chr7.symlink_to(source)

    return source, header_path
