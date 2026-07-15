# Structure fixtures for OpenMM oracle-validity (P3)

## Contents

| File | Role |
|---|---|
| `structure_manifest_v1.tsv` | Named subset `pdbbind_peptide_affinity_v1_structures_openmm` |
| `PP*__*_interface.pdb` | Trimmed public MHC–peptide interfaces (OpenMM/PDBFixer prepared) |
| `alanine_dipeptide_implicit.pdb` | OpenMM reference system for energy-conservation golden test |
| `oracle_validity_last_run.json` | Latest real Spearman / RMSE report (regenerate locally) |
| `*_full.pdb` | Optional source dumps (gitignored; rebuild interfaces with script) |

## Rebuild interfaces (no network)

```bash
cd core
poetry run python ../benchmarks/scripts/prepare_interfaces_from_full.py --cutoff 12.0
```

Requires local `*_full.pdb` beside this README.

## Run oracle-validity

```bash
cd core
poetry run python -m peptideforge.oracles.run_oracle_validity --mlflow
```

Reports a **real** Spearman vs experimental pK. Gate status is honest — see
`core/src/peptideforge/eval/ACCEPTANCE.md` (ρ ≥ 0.40). Latest run is written to
`oracle_validity_last_run.json` and logged to SQLite MLflow.
