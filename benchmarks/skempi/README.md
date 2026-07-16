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
