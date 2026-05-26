#!/usr/bin/env bash
# Restart all SN55 miners with current code and scripts/miner_env.sh settings.
set -euo pipefail
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/start_miners_screen.sh"
