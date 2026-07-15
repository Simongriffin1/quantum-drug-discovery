# PeptideForge Benchmark Report

_Generated at 2026-07-15T22:47:01.365681+00:00_ · git `894bf38c7522b568f520f8d8f4bd32586ebe73c7`

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

Subset `peptide_affinity_v2_experimental_openmm` N=31 partition=test: Spearman=0.1153 CI95%=[-0.270,0.475] (threshold ≥ 0.4, require CI_low>0, N≥30, red-team). Gate FAILED — reported honestly.

| Metric | Value | Source |
|---|---|---|
| `spearman` | 0.115323 | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `rmse` | 3301.39 | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `n` | 31 count | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `spearman_threshold` | 0.4 | `ACCEPTANCE.md` |
| `passed_threshold` | 0 | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07` |
| `spearman_ci_low` | -0.269822 | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `spearman_ci_high` | 0.474878 | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `pearson` | -0.0103843 | `file:benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json; mlflow_run_id=c89d92290ab84589ae60895542edbc07`, data_version=`peptide_affinity_v2_experimental_openmm` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json'
measurable: True
mlflow: {'run_id': 'c89d92290ab84589ae60895542edbc07', 'experiment_id': '1', 'tracking_uri': 'sqlite:////Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/peptide_affinity/data/oracle_validity_v2_mlflow.db'}
protocol: {'method': 'OpenMM MM-GBSA on prepared experimental structures', 'minimize_max_iterations': 0, 'platform': 'CPU', 'seed': 0, 'n_bootstrap': 1000, 'catalog': 'benchmarks/peptide_affinity/data/peptide_affinity_catalog_v2.tsv', 'prep_manifest': 'benchmarks/peptide_affinity/data/prepared/prep_manifest_scoreable.tsv', 'splits': 'benchmarks/peptide_affinity/data/splits_v2.json'}
record_ids: ['PA_2LSK_B', 'PA_2KOH_B', 'PA_4U2W_B', 'PA_5LY3_B', 'PA_6TYT_C', 'PA_4A2A_C', 'PA_4ODN_B', 'PA_1BXL_B', 'PA_2LP8_B', 'PA_3ZQI_C', 'PA_2N1G_B', 'PA_2MPS_B', 'PA_1G1E_A', 'PA_4ESG_C', 'PA_2N3K_B', 'PA_4J2C_B', 'PA_3V30_B', 'PA_3KTR_B', 'PA_5I22_B', 'PA_4WJQ_B', 'PA_4LG6_B', 'PA_4TZQ_B', 'PA_1F47_A', 'PA_4KMD_B', 'PA_2N9X_B', 'PA_2NM1_B', 'PA_6GZL_B', 'PA_6HOL_C', 'PA_2K00_B', 'PA_4CY1_C', 'PA_6F55_B']
red_team: {'passed': False, 'label_shuffle_passed': True, 'label_shuffle_rho': -0.08830645161290322, 'trivial_baseline_passed': False, 'model_rho': 0.1153225806451613, 'baseline_rho': 0.06318261512929894, 'leakage_passed': True, 'leakage_findings': ["trivial baseline(s) win or tie within Δρ: ['net_charge']; model_rho=0.1153 baselines={'peptide_length': 0.06318261512929894, 'net_charge': 0.4473853322373469}"], 'notes': "trivial baseline(s) win or tie within Δρ: ['net_charge']; model_rho=0.1153 baselines={'peptide_length': 0.06318261512929894, 'net_charge': 0.4473853322373469}"}
subset_name: 'peptide_affinity_v2_experimental_openmm'
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
