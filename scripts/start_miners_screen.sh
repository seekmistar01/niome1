#!/usr/bin/env bash
# Start all four Niome miners in screen with the pyenv Python that has bittensor.
set -euo pipefail

ROOT="/root/55miner/niome_1/subnet-niome"
PY="/root/.pyenv/versions/3.10.9/bin/python"
GATK="${NIOME_GATK:-/root/manualtest_55/gatk-4.6.2.0/gatk}"
CALLER="${NIOME_CALLER:-gatk}"

cd "$ROOT"
mkdir -p logs

pkill -f "neurons/miner.py --netuid 55" 2>/dev/null || true
sleep 2

"$PY" -c "import bittensor"  # fail fast if env is wrong

start_one() {
  local name="$1"
  shift
  screen -dmS "$name" bash -lc \
    "cd $ROOT && export PYTHONPATH=\$(pwd) NIOME_GATK=$GATK NIOME_CALLER=$CALLER && exec $PY neurons/miner.py $* >> logs/${name}.log 2>&1"
}

start_one niome_m1 \
  --netuid 55 --subtensor.network finney \
  --wallet.name seekmistar001 --wallet.hotkey seekmistar01 \
  --axon.port 20051 --logging.info

start_one niome_m2 \
  --netuid 55 --subtensor.network finney \
  --wallet.name seekmistar001 --wallet.hotkey seekmistar02 \
  --axon.port 20056 --logging.info

start_one niome_m3 \
  --netuid 55 --subtensor.network finney \
  --wallet.name seekmistar001 --wallet.hotkey seekmistar03 \
  --axon.port 20054 --logging.info

start_one niome_m4 \
  --netuid 55 --subtensor.network finney \
  --wallet.name seekmistar3 --wallet.hotkey seekmistar01 \
  --axon.port 20050 --logging.info

sleep 8
echo "=== screen ==="
screen -ls
echo "=== processes ==="
ps aux | grep "neurons/miner.py --netuid 55" | grep -v grep || true
echo "=== niome_m2 log (last 5 lines) ==="
tail -5 logs/niome_m2.log
