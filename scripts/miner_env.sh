#!/usr/bin/env bash
# Shared environment for live NIOME miners (source from start/restart scripts).
# Override any variable before sourcing if needed.

_MINER_ENV_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export NIOME_ROOT="${NIOME_ROOT:-$_MINER_ENV_ROOT}"
export PYTHONPATH="${NIOME_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
export NIOME_GATK="${NIOME_GATK:-${NIOME_ROOT}/tools/gatk/gatk}"
export NIOME_CALLER="${NIOME_CALLER:-gatk}"

# Sensitive multi-caller strategy (bcftools + FreeBayes + weakscan + rank-aligned submit).
export NIOME_USE_STRATEGY_AUTO="${NIOME_USE_STRATEGY_AUTO:-1}"
export NIOME_STRATEGY_SENSITIVITY="${NIOME_STRATEGY_SENSITIVITY:-sensitive}"
export NIOME_STRATEGY_SUBMIT_TOP_N="${NIOME_STRATEGY_SUBMIT_TOP_N:-30}"
export NIOME_STRATEGY_SUBMIT_ALLOWLIST="${NIOME_STRATEGY_SUBMIT_ALLOWLIST:-1}"
export NIOME_STRATEGY_USE_RANK_GT="${NIOME_STRATEGY_USE_RANK_GT:-1}"
export NIOME_STRATEGY_HYBRID_HOM_GT="${NIOME_STRATEGY_HYBRID_HOM_GT:-1}"
export NIOME_STRATEGY_DEEPVARIANT="${NIOME_STRATEGY_DEEPVARIANT:-0}"

# Evidence floors (strategy path uses light accept for strategy_auto sources in code).
export NIOME_EVIDENCE_MIN_ALT_DEPTH="${NIOME_EVIDENCE_MIN_ALT_DEPTH:-1}"

unset _MINER_ENV_ROOT
