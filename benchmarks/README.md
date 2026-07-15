# Benchmarks

Public benchmark **loaders** and **fixtures** for PeptideForge oracle-validity
and loop evaluation.

## Layout

```
benchmarks/
  fixtures/                     # local only — no live network in CI
    LICENSE.md
    pdbbind_peptide_affinity_v1.tsv
    skempi_ddg_v1.tsv
    clusters_v1.tsv
  peptideforge_benchmarks/      # Python package
  tests/
```

## Usage

```python
from peptideforge_benchmarks import (
    load_pdbbind_peptide_affinity,
    load_skempi_ddg,
    homology_aware_split,
)

records = load_pdbbind_peptide_affinity()
clusters = {r.record_id: r.cluster_id for r in records if r.cluster_id}
split = homology_aware_split(tuple(clusters), clusters, test_fraction=0.25, seed=0)
```

## License

See `fixtures/LICENSE.md`. Fixtures are format-faithful examples for CI plumbing.
They do **not** mean the oracle-validity stage gate has passed — that is P3+ with
a real physics oracle and pre-registered thresholds in
`core/src/peptideforge/eval/ACCEPTANCE.md`.
