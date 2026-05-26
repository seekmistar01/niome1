#!/usr/bin/env bash
# Download GRCh38 chromosome 7 from NCBI RefSeq and prepare data/chr7.fa for miners.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="${ROOT}/data"
RAW_GZ="${DATA}/NC_000007.14.fna.gz"
OUT_FA="${DATA}/chr7.fa"
VENV="${ROOT}/venv"
GATK="${ROOT}/tools/gatk/gatk"

# NCBI RefSeq GRCh38.p14 chromosome 7 (primary assembly)
NCBI_ACC="NC_000007.14"
EFETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${NCBI_ACC}&rettype=fasta&retmode=text"

mkdir -p "${DATA}"

echo "==> Downloading ${NCBI_ACC} from NCBI..."
curl -fsSL "${EFETCH_URL}" | gzip -c > "${RAW_GZ}"

echo "==> Writing ${OUT_FA} with contig name chr7 (subnet tasks use chr7:...)..."
gzip -dc "${RAW_GZ}" | awk '
  /^>/ { print ">chr7"; next }
  { print }
' > "${OUT_FA}.tmp"
mv "${OUT_FA}.tmp" "${OUT_FA}"

echo "==> Indexing (samtools + bwa)..."
command -v samtools >/dev/null && command -v bwa >/dev/null
rm -f "${OUT_FA}.fai" "${OUT_FA}.bwt" "${OUT_FA}.amb" "${OUT_FA}.ann" "${OUT_FA}.pac" "${OUT_FA}.sa" "${DATA}/chr7.dict"
samtools faidx "${OUT_FA}"
bwa index "${OUT_FA}"

if [[ -x "${GATK}" && -x "${VENV}/bin/python" ]]; then
  echo "==> GATK sequence dictionary..."
  export NIOME_GATK="${GATK}" PYTHONPATH="${ROOT}"
  ROOT="${ROOT}" "${VENV}/bin/python" - <<'PY'
import os
from pathlib import Path
from niome_subnet.cftr_miner_logic import _config_from_env, _ensure_reference_indexes

base = Path(os.environ["ROOT"])
_ensure_reference_indexes(_config_from_env(base), print)
PY
fi

ln -sfn chr7.fa "${DATA}/ref.fa"

echo "==> Done: ${OUT_FA}"
head -1 "${OUT_FA}"
wc -c "${OUT_FA}"
ls -lh "${OUT_FA}" "${OUT_FA}.fai" "${DATA}/chr7.dict" 2>/dev/null || true
