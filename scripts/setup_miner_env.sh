#!/usr/bin/env bash
# Install system tools, Python venv, GATK, and CFTR reference for NIOME miners.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="${ROOT}/venv"
GATK_DIR="${ROOT}/tools/gatk"
GATK_ZIP="${GATK_DIR}/gatk-4.6.2.0.zip"
GATK_URL="https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip"
# CHR7_SOURCE=ncbi (default, RefSeq NC_000007.14) or ucsc (hg38 chr7.fa.gz)
CHR7_SOURCE="${CHR7_SOURCE:-ncbi}"
CHR7_GZ="${ROOT}/data/chr7.fa.gz"
CHR7_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz"

export DEBIAN_FRONTEND=noninteractive

echo "==> Installing system packages..."
apt-get update -qq
apt-get install -y --no-install-recommends \
  screen \
  bwa \
  samtools \
  bcftools \
  tabix \
  default-jre-headless \
  curl \
  wget \
  unzip \
  ca-certificates \
  python3 \
  python3-venv \
  python3-pip

mkdir -p "${ROOT}/data" "${ROOT}/tools" "${ROOT}/logs" "${ROOT}/db/raw"

echo "==> Creating Python venv at ${VENV}..."
if [[ ! -x "${VENV}/bin/python" ]]; then
  python3 -m venv "${VENV}"
fi
"${VENV}/bin/pip" install --upgrade pip wheel
"${VENV}/bin/pip" install -r "${ROOT}/requirements.txt"

echo "==> Installing btcli (separate venv; conflicts with miner scalecodec)..."
BTCLI_VENV="${ROOT}/venv-btcli"
if [[ ! -x "${BTCLI_VENV}/bin/btcli" ]]; then
  python3 -m venv "${BTCLI_VENV}"
  "${BTCLI_VENV}/bin/pip" install --upgrade pip
  "${BTCLI_VENV}/bin/pip" install bittensor-cli
fi
cat > "${VENV}/bin/btcli" <<'WRAP'
#!/usr/bin/env bash
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec "${ROOT}/scripts/btcli.sh" "$@"
WRAP
chmod +x "${VENV}/bin/btcli"

echo "==> Installing GATK 4.6.2.0..."
if [[ ! -x "${GATK_DIR}/gatk" ]]; then
  mkdir -p "${GATK_DIR}"
  if [[ ! -f "${GATK_ZIP}" ]]; then
    wget -q --show-progress -O "${GATK_ZIP}" "${GATK_URL}"
  fi
  unzip -qo "${GATK_ZIP}" -d "${GATK_DIR}"
  gatk_bin=""
  for candidate in \
    "${GATK_DIR}/gatk-4.6.2.0/gatk" \
    "${GATK_DIR}/gatk-4.6.2.0/gatk.py" \
    "${GATK_DIR}/gatk"; do
    if [[ -f "$candidate" ]]; then
      gatk_bin="$candidate"
      break
    fi
  done
  if [[ -z "$gatk_bin" ]]; then
    echo "GATK binary not found after unzip" >&2
    find "${GATK_DIR}" -maxdepth 3 -type f 2>/dev/null | head -20 >&2
    exit 1
  fi
  chmod +x "$gatk_bin"
  ln -sf "$(realpath "$gatk_bin")" "${GATK_DIR}/gatk"
fi

echo "==> Fetching chr7 reference (source: ${CHR7_SOURCE})..."
if [[ ! -s "${ROOT}/data/chr7.fa" ]]; then
  if [[ "${CHR7_SOURCE}" == "ncbi" ]]; then
    bash "${ROOT}/scripts/fetch_chr7_ref_ncbi.sh"
  else
    wget -q --show-progress -O "${CHR7_GZ}" "${CHR7_URL}"
    gunzip -c "${CHR7_GZ}" > "${ROOT}/data/chr7.fa"
    rm -f "${CHR7_GZ}"
    echo "==> Indexing reference (samtools + bwa)..."
    samtools faidx "${ROOT}/data/chr7.fa"
    echo "    (bwa index on full chr7 may take several minutes...)"
    bwa index "${ROOT}/data/chr7.fa"
    export NIOME_GATK="${GATK_DIR}/gatk" PYTHONPATH="${ROOT}"
    if [[ ! -f "${ROOT}/data/chr7.dict" ]]; then
      ROOT="${ROOT}" NIOME_GATK="${NIOME_GATK}" PYTHONPATH="${ROOT}" \
        "${VENV}/bin/python" - <<'PY'
import os
from pathlib import Path
from niome_subnet.cftr_miner_logic import _config_from_env, _ensure_reference_indexes

base = Path(os.environ["ROOT"])
_ensure_reference_indexes(_config_from_env(base), print)
PY
    fi
    ln -sfn chr7.fa "${ROOT}/data/ref.fa"
  fi
else
  echo "    data/chr7.fa already present; skip download (re-fetch: bash scripts/fetch_chr7_ref_ncbi.sh)"
fi

echo "==> Verifying tools..."
for cmd in screen bwa samtools bcftools tabix java; do
  command -v "$cmd" >/dev/null || { echo "Missing: $cmd"; exit 1; }
done
"${GATK_DIR}/gatk" --version | head -1
"${VENV}/bin/python" -c "import bittensor; print('bittensor', bittensor.__version__)"

echo ""
echo "Setup complete."
echo "  ROOT:     ${ROOT}"
echo "  Python:   ${VENV}/bin/python"
echo "  GATK:     ${GATK_DIR}/gatk"
echo "  Reference: ${ROOT}/data/chr7.fa"
echo ""
echo "Next: create/register wallets, then:"
echo "  bash scripts/btcli.sh wallet list"
echo "  export NIOME_WALLET_MAIN=your_coldkey NIOME_WALLET_ALT=your_coldkey2  # optional"
echo "  bash scripts/start_miners_screen.sh"
