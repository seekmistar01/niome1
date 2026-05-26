#!/usr/bin/env bash
# Wallet CLI (separate venv — cannot share scalecodec with bittensor-cli).
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BTCLI_VENV="${ROOT}/venv-btcli"
if [[ ! -x "${BTCLI_VENV}/bin/btcli" ]]; then
  echo "Installing bittensor-cli into ${BTCLI_VENV}..." >&2
  python3 -m venv "${BTCLI_VENV}"
  "${BTCLI_VENV}/bin/pip" install -q --upgrade pip bittensor-cli
fi
exec "${BTCLI_VENV}/bin/btcli" "$@"
