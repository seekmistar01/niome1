#!/usr/bin/env python3
"""Prepare CFTR reference, ClinVar, and starter drug-response resources."""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path
from typing import List
from urllib.request import Request, urlopen


BASE_DIR = Path(__file__).resolve().parents[1]
REFERENCE_FASTA = BASE_DIR / "data" / "chr7.fa"
RAW_DIR = BASE_DIR / "db" / "raw"
DRUG_LABEL_DIR = BASE_DIR / "db" / "drug_labels"
PANEL_DIR = BASE_DIR / "panel"

CLINVAR_VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_TBI_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
VARIANT_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
)

DRUG_LABEL_URLS = {
    "ivacaftor_label.pdf": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2020/203188s029lbl.pdf",
    "tezacaftor_ivacaftor_label.pdf": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/210491s014lbl.pdf",
    "elexacaftor_tezacaftor_ivacaftor_label.pdf": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/212273s019lbl.pdf",
    "lumacaftor_ivacaftor_label.pdf": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2020/206038s012lbl.pdf",
}

STARTER_DRUG_ROWS = [
    ("7181", "non_responsive", "responsive", "responsive", "non_responsive"),
    ("53312", "non_responsive", "non_responsive", "non_responsive", "non_responsive"),
    ("53382", "responsive", "responsive", "responsive", "non_responsive"),
    ("2453648", "non_responsive", "non_responsive", "non_responsive", "non_responsive"),
    ("53502", "non_responsive", "non_responsive", "non_responsive", "non_responsive"),
    ("908328", "non_responsive", "non_responsive", "non_responsive", "non_responsive"),
    ("943683", "non_responsive", "non_responsive", "non_responsive", "non_responsive"),
    ("53856", "non_responsive", "non_responsive", "non_responsive", "non_responsive"),
]


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    DRUG_LABEL_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)

    ensure_reference_indexes(REFERENCE_FASTA)
    download_file(CLINVAR_VCF_URL, RAW_DIR / "clinvar_GRCh38.vcf.gz")
    download_file(CLINVAR_TBI_URL, RAW_DIR / "clinvar_GRCh38.vcf.gz.tbi")
    download_file(VARIANT_SUMMARY_URL, RAW_DIR / "variant_summary.txt.gz")
    build_clinvar_panel()
    download_drug_labels()
    write_starter_drug_panel(PANEL_DIR / "cftr_drug_response_by_id.csv")


def ensure_reference_indexes(reference_fasta: Path) -> None:
    if not reference_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")
    ensure_executable("samtools")
    ensure_executable("bwa")
    if not Path(f"{reference_fasta}.fai").exists():
        run(["samtools", "faidx", str(reference_fasta)])
    if not Path(f"{reference_fasta}.bwt").exists():
        run(["bwa", "index", str(reference_fasta)])


def download_file(url: str, output_path: Path) -> None:
    if output_path.exists() and output_path.stat().st_size > 0:
        print(f"Skipping existing {output_path}")
        return
    print(f"Downloading {url} -> {output_path}")
    request = Request(url, headers={"User-Agent": "niome-cftr-resource-prep/1.0"})
    with urlopen(request, timeout=300) as response, output_path.open("wb") as output_file:
        shutil.copyfileobj(response, output_file)


def build_clinvar_panel() -> None:
    ensure_executable("bcftools")
    ensure_executable("tabix")
    raw_vcf = RAW_DIR / "clinvar_GRCh38.vcf.gz"
    normalized_vcf = PANEL_DIR / "clinvar_cftr.norm.vcf.gz"
    panel_tsv = PANEL_DIR / "cftr_clinvar_panel.tsv"

    view_process = subprocess.Popen(
        ["bcftools", "view", "-r", "7:117480000-117670000", str(raw_vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    norm_process = subprocess.Popen(
        [
            "bcftools",
            "norm",
            "-m",
            "-both",
            "-Oz",
            "-o",
            str(normalized_vcf),
        ],
        stdin=view_process.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if view_process.stdout is not None:
        view_process.stdout.close()
    norm_stdout, norm_stderr = norm_process.communicate()
    view_stderr = view_process.stderr.read() if view_process.stderr is not None else b""
    view_returncode = view_process.wait()
    if view_returncode != 0:
        raise RuntimeError(view_stderr.decode(errors="replace"))
    if norm_process.returncode != 0:
        detail = (norm_stderr or norm_stdout).decode(errors="replace")
        raise RuntimeError(detail)

    run(["tabix", "-f", "-p", "vcf", str(normalized_vcf)])
    query_format = "%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CLNHGVS\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\t%INFO/CLNDN\\n"
    with panel_tsv.open("w", encoding="utf-8") as output_file:
        subprocess.run(
            ["bcftools", "query", "-f", query_format, str(normalized_vcf)],
            check=True,
            text=True,
            stdout=output_file,
        )


def download_drug_labels() -> None:
    for filename, url in DRUG_LABEL_URLS.items():
        download_file(url, DRUG_LABEL_DIR / filename)


def write_starter_drug_panel(output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(
            [
                "variation_id",
                "ivacaftor",
                "tezacaftor_ivacaftor",
                "elexacaftor_tezacaftor_ivacaftor",
                "lumacaftor_ivacaftor",
            ]
        )
        writer.writerows(STARTER_DRUG_ROWS)


def ensure_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def run(command: List[str]) -> None:
    print("+ " + " ".join(command))
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
