# PeptideForge Benchmark Report

_Generated at 2026-07-15T21:35:51.455907+00:00_ · git `unknown`

## Caveats

- Oracle-validity FAIL means binding campaigns are not authorized (CURSOR_PROJECT_CONTEXT §8).
- Surrogate / loop / acquisition sections marked synthetic_* validate algorithmic plumbing only — not physics fidelity.
- All numbers below are traced to JSON artifacts and/or MLflow run IDs; none are invented.

## Gate summary

| Section | Status |
|---|---|
| Oracle validity (affinity) | **FAIL** |
| Stability / ddG (SKEMPI) | **NOT_RUN** |
| Surrogate calibration (synthetic plumbing) | **PASS** |
| Acquisition (Branin–Currin hypervolume) | **PASS** |
| Loop efficiency (simulations-to-target) | **PASS** |

## Oracle validity (affinity)

**Status:** FAIL

Subset `pdbbind_peptide_affinity_v1_structures_openmm` N=5: Spearman=-0.3000 (threshold ≥ 0.4). Gate FAILED — reported honestly.

| Metric | Value | Source |
|---|---|---|
| `spearman` | -0.3 | `file:benchmarks/fixtures/structures/oracle_validity_last_run.json; mlflow_run_id=1a1889d0aae54001b9026066f6bc2024`, data_version=`pdbbind_peptide_affinity_v1_structures_openmm` |
| `rmse` | 23.344 | `file:benchmarks/fixtures/structures/oracle_validity_last_run.json; mlflow_run_id=1a1889d0aae54001b9026066f6bc2024`, data_version=`pdbbind_peptide_affinity_v1_structures_openmm` |
| `n` | 5 count | `file:benchmarks/fixtures/structures/oracle_validity_last_run.json; mlflow_run_id=1a1889d0aae54001b9026066f6bc2024`, data_version=`pdbbind_peptide_affinity_v1_structures_openmm` |
| `spearman_threshold` | 0.4 | `ACCEPTANCE.md` |
| `passed_threshold` | 0 | `file:benchmarks/fixtures/structures/oracle_validity_last_run.json; mlflow_run_id=1a1889d0aae54001b9026066f6bc2024` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/fixtures/structures/oracle_validity_last_run.json'
mlflow: {'run_id': '1a1889d0aae54001b9026066f6bc2024', 'experiment_id': '1', 'tracking_uri': 'sqlite:////Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/fixtures/structures/oracle_validity_mlflow.db'}
protocol: {'method': 'OpenMM MM-GBSA (amber14 + OBC2), ΔG ≈ E_c − E_r − E_l after minimize', 'invert_sign': True, 'minimize_max_iterations': 0, 'platform': 'CPU', 'seed': 0, 'manifest': '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/fixtures/structures/structure_manifest_v1.tsv', 'structures': 'trimmed public MHC–peptide interface PDBs (see REMARK headers)'}
record_ids: ['PP001', 'PP003', 'PP005', 'PP009', 'PP011']
subset_name: 'pdbbind_peptide_affinity_v1_structures_openmm'
```

</details>


## Stability / ddG (SKEMPI)

**Status:** NOT_RUN

No logged SKEMPI ΔΔG oracle-validity run. Threshold pre-registered in ACCEPTANCE.md (Spearman ≥ 0.30). Refusing to invent a correlation.

| Metric | Value | Source |
|---|---|---|
| `ddg_spearman_threshold` | 0.3 | `ACCEPTANCE.md` |

<details><summary>Details (JSON-backed)</summary>

```
fixture: 'benchmarks/fixtures/skempi_ddg_v1.tsv'
```

</details>


## Surrogate calibration (synthetic plumbing)

**Status:** PASS

ECE=0.0375 (threshold < 0.1); red-team passed=True. This validates UQ plumbing on synthetic_* labels — NOT physics affinity.

| Metric | Value | Source |
|---|---|---|
| `ece` | 0.0375 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/surrogate_acceptance_last_run.json`, data_version=`synthetic_surrogate_v1` |
| `empirical_coverage` | 0.9375 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/surrogate_acceptance_last_run.json`, data_version=`synthetic_surrogate_v1` |
| `ece_threshold` | 0.1 | `ACCEPTANCE.md` |
| `red_team_passed` | 1 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/surrogate_acceptance_last_run.json` |
| `model_rho` | 0.979412 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/surrogate_acceptance_last_run.json`, data_version=`synthetic_surrogate_v1` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/surrogate_acceptance_last_run.json'
data_version: 'synthetic_surrogate_v1'
notes: 'synthetic_learnable_oracle_not_physics'
```

</details>


## Acquisition (Branin–Currin hypervolume)

**Status:** PASS

qNEHVI HV=2067.0869 vs random HV=1868.4346 (same budget). Validates acquisition logic, not biology.

| Metric | Value | Source |
|---|---|---|
| `hv_qnehvi` | 2067.09 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/acquisition_branin_currin_last_run.json` — synthetic_branin_currin_perfect_surrogate |
| `hv_random` | 1868.43 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/acquisition_branin_currin_last_run.json` — synthetic_branin_currin_perfect_surrogate |
| `acquisition_passed` | 1 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/acquisition_branin_currin_last_run.json` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/acquisition_branin_currin_last_run.json'
batch_size: 2
n_init: 6
n_pool: 60
n_rounds: 5
reference: [-193.8906704566085, -14.311968516995535]
```

</details>


## Loop efficiency (simulations-to-target)

**Status:** PASS

qNEHVI calls=12 best=-4.553916231220688; random calls=18 best=-4.553916231220688; target=-4.5. Mode=synthetic_physics_label_shared_pool.

| Metric | Value | Source |
|---|---|---|
| `oracle_calls_qnehvi` | 12 count | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json`, data_version=`synthetic_physics_label_shared_pool` |
| `oracle_calls_random` | 18 count | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json`, data_version=`synthetic_physics_label_shared_pool` |
| `best_qnehvi` | -4.55392 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json`, data_version=`synthetic_physics_label_shared_pool` |
| `best_random` | -4.55392 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json`, data_version=`synthetic_physics_label_shared_pool` |
| `target_value` | -4.5 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json` |
| `spend_gate_passed` | 1 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/reports/artifacts/loop_spend_gate_last_run.json'
notes: 'simulation_mode_spend_gate_shared_pool'
```

</details>

---

Regenerate with:

```bash
make benchmark-report
# or: cd core && poetry run python -m peptideforge.reports.generate_benchmark_report
```
