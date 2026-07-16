#!/usr/bin/env bash
# Step 4 — predicted-fold degradation (run from repo root).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

# OpenMM / pdbfixer live in the Poetry venv; Boltz CLI in Anaconda base (adjust if needed).
export PATH="${HOME}/anaconda3/bin:${PATH}"
PY="${PY:-${HOME}/Library/Caches/pypoetry/virtualenvs/peptideforge-IqirF6bS-py3.11/bin/python}"
export PYTHONPATH=benchmarks:core/src

if ! command -v boltz >/dev/null; then
  echo "FAIL: boltz not on PATH. Install: pip install boltz" >&2
  exit 1
fi
if ! "$PY" -c "import openmm, pdbfixer" 2>/dev/null; then
  echo "FAIL: openmm+pdbfixer missing in $PY" >&2
  exit 1
fi
if ! python -c "import torch" 2>/dev/null; then
  echo "FAIL: torch broken in Anaconda (boltz needs it). Try: pip install --force-reinstall 'torch>=2.2'" >&2
  exit 1
fi

echo "== predict_folds (5 unique PDBs; CPU Boltz is slow) =="
"$PY" benchmarks/skempi/predict_folds.py "$@"

echo "== run_fold_degradation =="
"$PY" benchmarks/skempi/run_fold_degradation.py

echo "== stratify =="
"$PY" benchmarks/skempi/stratify.py

echo "== authorization_build =="
"$PY" -m peptideforge.authorization_build

echo "Done. See benchmarks/skempi/data/fold_degradation_last_run.json"
