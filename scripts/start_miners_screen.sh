#!/usr/bin/env bash
# Start four NIOME miners in detached screen sessions.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=miner_env.sh
source "${ROOT}/scripts/miner_env.sh"
VENV="${ROOT}/venv"
PY="${NIOME_PYTHON:-${VENV}/bin/python}"
GATK="${NIOME_GATK}"
CALLER="${NIOME_CALLER}"

# Registered on SN55 (finney):
#   seekmistar001 / seekmistar01  UID 145  -> 50007
#   seekmistar001 / seekmistar02  UID 120  -> 50006
#   seekmistar001 / seekmistar03  UID 116  -> 50005
#   seekmistar3   / seekmistar01  UID 7    -> 50004
WALLET_MAIN="${NIOME_WALLET_MAIN:-seekmistar001}"
WALLET_ALT="${NIOME_WALLET_ALT:-seekmistar3}"

die() { echo "ERROR: $*" >&2; exit 1; }

[[ -x "$PY" ]] || die "Python not found at $PY — run: bash scripts/setup_miner_env.sh"
[[ -x "$GATK" ]] || die "GATK not found at $GATK — run: bash scripts/setup_miner_env.sh"
[[ -f "${ROOT}/data/chr7.fa" ]] || die "Missing ${ROOT}/data/chr7.fa — run: bash scripts/setup_miner_env.sh"
for suffix in .fai .amb .ann .bwt .pac .sa; do
  [[ -f "${ROOT}/data/chr7.fa${suffix}" ]] || die "Missing BWA/FAI index ${ROOT}/data/chr7.fa${suffix} — run: bash scripts/fetch_chr7_ref_ncbi.sh"
done
command -v screen >/dev/null || die "screen not installed — run: bash scripts/setup_miner_env.sh"
for cmd in bwa samtools bcftools tabix java; do
  command -v "$cmd" >/dev/null || die "$cmd not on PATH — run: bash scripts/setup_miner_env.sh"
done

cd "$ROOT"
mkdir -p logs

"$PY" -c "import bittensor"

pkill -f "neurons/miner.py --netuid 55" 2>/dev/null || true
sleep 2

start_one() {
  local name="$1"
  local instance="$2"
  shift 2
  screen -dmS "$name" bash -lc \
    "cd '$ROOT' && source scripts/miner_env.sh && export MINER_INSTANCE='$instance' && exec '$PY' neurons/miner.py $* >> logs/${name}.log 2>&1"
}

start_one niome_m1 niome_m1 \
  --netuid 55 --subtensor.network finney \
  --wallet.name "$WALLET_MAIN" --wallet.hotkey seekmistar01 \
  --axon.port 50007 --logging.info

start_one niome_m2 niome_m2 \
  --netuid 55 --subtensor.network finney \
  --wallet.name "$WALLET_MAIN" --wallet.hotkey seekmistar02 \
  --axon.port 50006 --logging.info

start_one niome_m3 niome_m3 \
  --netuid 55 --subtensor.network finney \
  --wallet.name "$WALLET_MAIN" --wallet.hotkey seekmistar03 \
  --axon.port 50005 --logging.info

start_one niome_m4 niome_m4 \
  --netuid 55 --subtensor.network finney \
  --wallet.name "$WALLET_ALT" --wallet.hotkey seekmistar01 \
  --axon.port 50004 --logging.info

sleep 8
echo "=== screen ==="
screen -ls
echo "=== processes ==="
ps aux | grep "neurons/miner.py --netuid 55" | grep -v grep || true
echo "=== strategy env ==="
echo "  NIOME_USE_STRATEGY_AUTO=${NIOME_USE_STRATEGY_AUTO}"
echo "  NIOME_STRATEGY_SENSITIVITY=${NIOME_STRATEGY_SENSITIVITY}"
echo "  NIOME_STRATEGY_SUBMIT_TOP_N=${NIOME_STRATEGY_SUBMIT_TOP_N}"
echo "  NIOME_STRATEGY_USE_RANK_GT=${NIOME_STRATEGY_USE_RANK_GT}"
echo "=== niome_m1 log (last 15 lines) ==="
tail -15 logs/niome_m1.log 2>/dev/null || echo "(no log yet)"
