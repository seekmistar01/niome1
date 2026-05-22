#!/usr/bin/env python3
"""Prepare CFTR drug-response panel inputs for miner annotations."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Dict, List, Set, Tuple


DRUG_KEYS = [
    "ivacaftor",
    "tezacaftor_ivacaftor",
    "elexacaftor_tezacaftor_ivacaftor",
    "lumacaftor_ivacaftor",
]


LABELS = {
    "ivacaftor": {
        "pdf": "db/drug_labels/kalydeco_ivacaftor.pdf",
        "txt": "db/drug_labels/kalydeco_ivacaftor.txt",
        "url": "https://pi.vrtx.com/files/uspi_ivacaftor.pdf",
    },
    "tezacaftor_ivacaftor": {
        "pdf": "db/drug_labels/symdeko_tezacaftor_ivacaftor.pdf",
        "txt": "db/drug_labels/symdeko_tezacaftor_ivacaftor.txt",
        "url": "https://pi.vrtx.com/files/uspi_tezacaftor_ivacaftor.pdf",
    },
    "elexacaftor_tezacaftor_ivacaftor": {
        "pdf": "db/drug_labels/trikafta_elexacaftor_tezacaftor_ivacaftor.pdf",
        "txt": "db/drug_labels/trikafta_elexacaftor_tezacaftor_ivacaftor.txt",
        "url": "https://pi.vrtx.com/files/uspi_elexacaftor_tezacaftor_ivacaftor.pdf",
    },
    "lumacaftor_ivacaftor": {
        "pdf": "db/drug_labels/orkambi_lumacaftor_ivacaftor.pdf",
        "txt": "db/drug_labels/orkambi_lumacaftor_ivacaftor.txt",
        "url": "https://pi.vrtx.com/files/uspi_lumacaftor_ivacaftor.pdf",
    },
}


# Starter rows match the sample annotation format and can be expanded over time.
STARTER_BY_ID = {
    "7181": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "53312": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "53382": {
        "ivacaftor": "responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "2453648": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "53502": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "908328": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "943683": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "53856": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
}


AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Stop": "*",
}


def download(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"[OK] exists: {out_path}")
        return
    print(f"[DOWNLOAD] {url} -> {out_path}")
    urllib.request.urlretrieve(url, out_path)


def run(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    print("[CMD]", " ".join(map(str, cmd)))
    return subprocess.run(cmd, text=True, check=check)


def normalize_text(s: str) -> str:
    return (
        s.replace("\u2192", ">")
        .replace("\u2212", "-")
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u00a0", " ")
        .replace(" ", "")
        .upper()
    )


def convert_pdf_to_text(pdf: Path, txt: Path) -> None:
    if txt.exists() and txt.stat().st_size > 0:
        print(f"[OK] text exists: {txt}")
        return
    if shutil.which("pdftotext") is None:
        print("[WARN] pdftotext missing. Install: apt-get install -y poppler-utils")
        return
    run(["pdftotext", "-layout", str(pdf), str(txt)], check=False)


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def load_cftr_clinvar_ids(panel_path: Path) -> Set[str]:
    ids: Set[str] = set()
    with panel_path.open("r", encoding="utf-8", errors="replace") as panel_file:
        for line in panel_file:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            vid = parts[0].strip()
            if vid and vid != ".":
                ids.add(vid)
    return ids


def read_variant_summary(variant_summary_path: Path, wanted_ids: Set[str]) -> Dict[str, Dict[str, str]]:
    result: Dict[str, Dict[str, str]] = {}
    with gzip.open(variant_summary_path, "rt", encoding="utf-8", errors="replace") as summary_file:
        header = summary_file.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}

        def get(row: List[str], name: str) -> str:
            i = idx.get(name)
            if i is None or i >= len(row):
                return ""
            return row[i]

        for line in summary_file:
            row = line.rstrip("\n").split("\t")
            vid = get(row, "VariationID")
            if vid not in wanted_ids:
                continue
            result[vid] = {
                "variation_id": vid,
                "name": get(row, "Name"),
                "gene": get(row, "GeneSymbol"),
                "clinical_significance": get(row, "ClinicalSignificance"),
                "phenotype": get(row, "PhenotypeList"),
                "protein_change": get(row, "ProteinChange"),
            }
    return result


def short_protein_token(token: str) -> str:
    token = token.strip().replace("p.", "").replace("(", "").replace(")", "")

    match = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", token)
    if match:
        a1, pos, a2 = match.groups()
        if a1 in AA3_TO_1 and a2 in AA3_TO_1:
            return f"{AA3_TO_1[a1]}{pos}{AA3_TO_1[a2]}"

    match = re.match(r"^([A-Z][a-z]{2})(\d+)del$", token)
    if match:
        a1, pos = match.groups()
        if a1 in AA3_TO_1:
            return f"{AA3_TO_1[a1]}{pos}del"

    return token


def aliases_from_summary(summary: Dict[str, str]) -> Set[str]:
    aliases: Set[str] = set()
    text = " ".join([summary.get("name", ""), summary.get("protein_change", "")])

    for match in re.findall(r"c\.[A-Za-z0-9_\-\+\*>]+", text):
        aliases.add(match)
        aliases.add(match.replace("c.", ""))

    for match in re.findall(r"(?:p\.)?[A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|del|Ter|Stop|\*)", text):
        aliases.add(match)
        aliases.add(short_protein_token(match))

    for match in re.findall(r"\b[A-Z]\d{1,4}(?:[A-Z]|\*|del)\b", text):
        aliases.add(match)

    return {alias for alias in aliases if alias and len(alias) >= 3}


def infer_auto_response(
    aliases: Set[str],
    label_texts: Dict[str, str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    response = {key: "non_responsive" for key in DRUG_KEYS}
    evidence = {key: "" for key in DRUG_KEYS}
    normalized_aliases = {normalize_text(alias) for alias in aliases}

    for drug_key in DRUG_KEYS:
        label_norm = normalize_text(label_texts.get(drug_key, ""))
        hits = [alias for alias in normalized_aliases if alias and alias in label_norm]
        if hits:
            response[drug_key] = "responsive"
            evidence[drug_key] = ",".join(sorted(set(hits))[:20])

    # ORKAMBI is indicated for homozygous F508del; per-variant annotation marks
    # F508del as responsive and leaves genotype-level handling to scoring logic.
    if "F508DEL" in normalized_aliases:
        response["lumacaftor_ivacaftor"] = "responsive"
        evidence["lumacaftor_ivacaftor"] = "F508DEL"

    return response, evidence


def write_response_csv(path: Path, rows: Dict[str, Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=["variation_id"] + DRUG_KEYS)
        writer.writeheader()
        for vid in sorted(rows, key=lambda item: int(item) if item.isdigit() else item):
            row = {"variation_id": vid}
            row.update({key: rows[vid].get(key, "non_responsive") for key in DRUG_KEYS})
            writer.writerow(row)


def write_candidate_review(
    path: Path,
    summaries: Dict[str, Dict[str, str]],
    auto_rows: Dict[str, Dict[str, str]],
    evidence_rows: Dict[str, Dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output_file:
        fieldnames = [
            "variation_id",
            "name",
            "protein_change",
            "aliases",
            "ivacaftor",
            "ivacaftor_evidence",
            "tezacaftor_ivacaftor",
            "tezacaftor_ivacaftor_evidence",
            "elexacaftor_tezacaftor_ivacaftor",
            "elexacaftor_tezacaftor_ivacaftor_evidence",
            "lumacaftor_ivacaftor",
            "lumacaftor_ivacaftor_evidence",
        ]
        writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for vid in sorted(summaries, key=lambda item: int(item) if item.isdigit() else item):
            summary = summaries[vid]
            aliases = sorted(aliases_from_summary(summary))
            resp = auto_rows.get(vid, {key: "non_responsive" for key in DRUG_KEYS})
            ev = evidence_rows.get(vid, {key: "" for key in DRUG_KEYS})
            writer.writerow(
                {
                    "variation_id": vid,
                    "name": summary.get("name", ""),
                    "protein_change": summary.get("protein_change", ""),
                    "aliases": ",".join(aliases),
                    "ivacaftor": resp["ivacaftor"],
                    "ivacaftor_evidence": ev["ivacaftor"],
                    "tezacaftor_ivacaftor": resp["tezacaftor_ivacaftor"],
                    "tezacaftor_ivacaftor_evidence": ev["tezacaftor_ivacaftor"],
                    "elexacaftor_tezacaftor_ivacaftor": resp[
                        "elexacaftor_tezacaftor_ivacaftor"
                    ],
                    "elexacaftor_tezacaftor_ivacaftor_evidence": ev[
                        "elexacaftor_tezacaftor_ivacaftor"
                    ],
                    "lumacaftor_ivacaftor": resp["lumacaftor_ivacaftor"],
                    "lumacaftor_ivacaftor_evidence": ev["lumacaftor_ivacaftor"],
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clinvar-panel", default="panel/cftr_clinvar_panel.tsv")
    parser.add_argument("--variant-summary", default="db/raw/variant_summary.txt.gz")
    parser.add_argument("--out-manual", default="panel/cftr_drug_response_by_id.csv")
    parser.add_argument("--out-auto", default="panel/cftr_drug_response_by_id.auto.csv")
    parser.add_argument("--out-review", default="panel/cftr_drug_response_candidates.tsv")
    parser.add_argument("--accept-auto", action="store_true")
    parser.add_argument("--overwrite-manual", action="store_true")
    args = parser.parse_args()

    clinvar_panel = Path(args.clinvar_panel)
    variant_summary = Path(args.variant_summary)
    if not clinvar_panel.exists():
        raise FileNotFoundError(f"Missing ClinVar panel: {clinvar_panel}")
    if not variant_summary.exists():
        raise FileNotFoundError(f"Missing variant summary: {variant_summary}")

    print("[1/5] Download drug labels")
    for meta in LABELS.values():
        download(meta["url"], Path(meta["pdf"]))

    print("[2/5] Convert PDFs to text")
    for meta in LABELS.values():
        convert_pdf_to_text(Path(meta["pdf"]), Path(meta["txt"]))

    print("[3/5] Load CFTR ClinVar IDs and variant names")
    wanted_ids = load_cftr_clinvar_ids(clinvar_panel)
    summaries = read_variant_summary(variant_summary, wanted_ids)
    print(f"[INFO] CFTR ClinVar IDs: {len(wanted_ids)}")
    print(f"[INFO] Variant summary rows matched: {len(summaries)}")

    print("[4/5] Build automatic candidate drug-response table")
    label_texts = {key: read_text(Path(meta["txt"])) for key, meta in LABELS.items()}
    auto_rows: Dict[str, Dict[str, str]] = {}
    evidence_rows: Dict[str, Dict[str, str]] = {}
    for vid, summary in summaries.items():
        aliases = aliases_from_summary(summary)
        response, evidence = infer_auto_response(aliases, label_texts)
        auto_rows[vid] = response
        evidence_rows[vid] = evidence

    final_manual = dict(auto_rows)
    for vid, row in STARTER_BY_ID.items():
        final_manual[vid] = row

    write_response_csv(Path(args.out_auto), final_manual)
    write_candidate_review(Path(args.out_review), summaries, auto_rows, evidence_rows)

    print("[5/5] Write manual runtime panel")
    manual_path = Path(args.out_manual)
    if manual_path.exists() and not args.overwrite_manual and not args.accept_auto:
        print(f"[OK] manual panel exists, not overwriting: {manual_path}")
        print("Use --overwrite-manual or --accept-auto if you want to replace it.")
    else:
        write_response_csv(manual_path, final_manual)
        print(f"[OK] wrote runtime drug panel: {manual_path}")

    print("")
    print("[DONE]")
    print(f"Runtime panel:   {manual_path}")
    print(f"Auto panel:      {args.out_auto}")
    print(f"Review file:     {args.out_review}")
    print("")
    print("Recommended next step:")
    print(f"  less {args.out_review}")
    print("")
    print("Use the runtime panel in miner:")
    print("  export NIOME_CFTR_DRUG_PANEL=$(pwd)/panel/cftr_drug_response_by_id.csv")


if __name__ == "__main__":
    main()
