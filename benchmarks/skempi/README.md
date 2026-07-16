# SKEMPI within-target ΔΔG benchmark

Product-relevant co-primary gate (see `core/src/peptideforge/eval/ACCEPTANCE.md`).
Cross-target affinity is **not** replaced by this axis.

## Data (manual download — fail loud if missing)

1. Download SKEMPI v2.0 from https://life.bsc.es/pid/skempi2/
2. Convert / subset to a TSV with columns required by
   `peptideforge_benchmarks.skempi.load_skempi_ddg`
   (`record_id`, `pdb_id`, `mutant`, `ddg_kcal_mol`, …).
3. Place WT crystal structures as `{pdb_id}.pdb` under `data/structures/`.

Fixture-only path for plumbing: `benchmarks/fixtures/skempi_ddg_v1.tsv`
(does **not** constitute a gate PASS).

## Run

```bash
PYTHONPATH=benchmarks:core/src python -m benchmarks.skempi.run_skempi_ddg \
  --skempi-tsv benchmarks/skempi/data/skempi_v2.tsv \
  --structure-dir benchmarks/skempi/data/structures \
  --out benchmarks/skempi/data/skempi_ddg_last_run.json \
  --gb-model gbn2 --solute-dielectric 1.0
```

Or: `python benchmarks/skempi/run_skempi_ddg.py ...`

## Step 4 — Fold→score degradation

**Run from the repo root** (not `~`). Paths like `benchmarks/skempi/...` are relative to
`quantum-drug-discovery/`.

```bash
cd /path/to/quantum-drug-discovery/quantum-drug-discovery

# One env: Poetry python (OpenMM) + Anaconda boltz on PATH
export PATH="$HOME/anaconda3/bin:$PATH"
PY="$HOME/Library/Caches/pypoetry/virtualenvs/peptideforge-IqirF6bS-py3.11/bin/python"
export PYTHONPATH=benchmarks:core/src

# Smoke (1 hold-out row → 1 WT Boltz fold):
"$PY" benchmarks/skempi/predict_folds.py --max-complexes 1

# Full hold-out (5 unique PDBs × Boltz CPU — expect hours):
"$PY" benchmarks/skempi/predict_folds.py
"$PY" benchmarks/skempi/run_fold_degradation.py
"$PY" benchmarks/skempi/stratify.py
"$PY" -m peptideforge.authorization_build
```

Or: `bash benchmarks/skempi/run_step4.sh` (optional `--max-complexes 1` forwarded to predict).

If `boltz` fails on `import torch` with `libtorch_cpu.dylib` missing, reinstall:
`pip install --force-reinstall 'torch>=2.2'` in the same env that owns `boltz`.

Matched hold-out: `data/skempi_powered_holdout_v1.json` (same N=100 as experimental PASS).
