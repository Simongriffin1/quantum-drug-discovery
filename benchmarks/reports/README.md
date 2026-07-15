# Benchmark reports (P12)

Generated credibility artifacts. Regenerate with:

```bash
make benchmark-report
```

- `benchmark_report.md` — human-readable report (oracle Spearmen, ECE, loop spend, …)
- `benchmark_report.json` — machine-readable `BenchmarkReport`
- `artifacts/` — JSON sources for synthetic_* sections (surrogate / acquisition / loop)

**Oracle affinity numbers come from** `../fixtures/structures/oracle_validity_last_run.json`
(and optional MLflow run IDs). Fake numbers are not allowed.
