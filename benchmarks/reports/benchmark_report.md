# PeptideForge Benchmark Report

_Generated at 2026-07-16T02:31:46.810065+00:00_ · git `ea0d919ed7a918f26e9d0dbc9c804ca2be61beab`

## Caveats

- Oracle-validity FAIL means binding campaigns are not authorized (CURSOR_PROJECT_CONTEXT §8).
- Surrogate / loop / acquisition sections marked synthetic_* validate algorithmic plumbing only — not physics fidelity.
- All numbers below are traced to JSON artifacts and/or MLflow run IDs; none are invented.

## Gate summary

| Section | Status |
|---|---|
| Oracle validity (affinity) | **FAIL** |
| Stability / ddG (SKEMPI within-target) | **FAIL** |
| Surrogate calibration (synthetic plumbing) | **PASS** |
| Acquisition (Branin–Currin hypervolume) | **PASS** |
| Loop efficiency (simulations-to-target) | **PASS** |

## Oracle validity (affinity)

**Status:** FAIL

Subset `peptide_affinity_v2_experimental_openmm` N=40 partition=test: Spearman=0.1870 CI95%=[-0.160,0.514] (threshold ≥ 0.4, require CI_low>0, N≥30, red-team). Gate FAILED — reported honestly.

| Metric | Value | Source |
|---|---|---|
| `spearman` | 0.186987 | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `rmse` | 9462.18 | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `n` | 40 count | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `spearman_threshold` | 0.4 | `ACCEPTANCE.md` |
| `passed_threshold` | 0 | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001` |
| `spearman_ci_low` | -0.159511 | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `spearman_ci_high` | 0.513896 | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001`, data_version=`peptide_affinity_v2_experimental_openmm` |
| `pearson` | -0.0437044 | `file:benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json; mlflow_run_id=aeefc1bdfd1d4e289b923912c8639001`, data_version=`peptide_affinity_v2_experimental_openmm` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: 'benchmarks/peptide_affinity/data/oracle_validity_v3_oneshot_test.json'
measurable: True
mlflow: {'run_id': 'aeefc1bdfd1d4e289b923912c8639001', 'experiment_id': '1', 'tracking_uri': 'sqlite:///benchmarks/peptide_affinity/data/oracle_validity_v3.db'}
protocol: {'name': 'gbsa_gbn2_eps1_salt0_min0', 'gb_model': 'gbn2', 'solute_dielectric': 1.0, 'salt_conc_M': 0.0, 'minimize_max_iterations': 0, 'chosen_from': 'elec_sweep_traindev_v3 + entropy_waters_ab (no retune after test)', 'method': 'OpenMM MM-GBSA on prepared experimental structures', 'platform': 'CPU', 'seed': 0, 'n_bootstrap': 1000, 'catalog': 'benchmarks/peptide_affinity/data/peptide_affinity_catalog_v3.tsv', 'prep_manifest': 'benchmarks/peptide_affinity/data/prepared/prep_manifest_scoreable.tsv', 'splits': 'benchmarks/peptide_affinity/data/splits_v3.json', 'oneshot_test': True}
record_ids: ['PA_3PKN_B', 'PA_2LOZ_B', 'PA_5H5S_B', 'PA_4FI9_B', 'PA_2M41_A', 'PA_4RRV_B', 'PA_6HOL_C', 'PA_2K00_B', 'PA_4CY1_C', 'PA_5CQX_C', 'PA_3BU3_B', 'PA_4MNY_C', 'PA_4MNX_B', 'PA_4MNW_B', 'PA_5WQD_H', 'PA_2KOH_B', 'PA_3NFL_E', 'PA_2KYM_B', 'PA_5V6Y_E', 'PA_5IY4_B', 'PA_4A2A_C', 'PA_4ODN_B', 'PA_1G1E_A', 'PA_1ABT_B', 'PA_4ESG_C', 'PA_2N3K_B', 'PA_3KTR_B', 'PA_5I22_B', 'PA_4WJQ_B', 'PA_4LG6_B', 'PA_3QN7_B', 'PA_6E49_D', 'PA_4KMD_B', 'PA_2N9X_B', 'PA_2NM1_B', 'PA_6GZL_B', 'PA_6O21_B', 'PA_5CQY_C', 'PA_3BU5_B', 'PA_1S5Q_A']
red_team: {'passed': False, 'label_shuffle_passed': True, 'label_shuffle_rho': -0.17113102388354878, 'trivial_baseline_passed': False, 'model_rho': 0.18698691370609252, 'baseline_rho': 0.20460044983108497, 'leakage_passed': False, 'leakage_findings': ["identical sequence leakage: ['PA_2KYM_B']", "sequence identity >30% leakage: ['PA_5H5S_B~PA_3V30_B', 'PA_4CY1_C~PA_3ZQI_C', 'PA_4MNY_C~PA_4U2W_B', 'PA_4MNW_B~PA_4U2W_B', 'PA_3NFL_E~PA_1YBO_C', 'PA_2KYM_B~PA_1YBO_C', 'PA_4LG6_B~PA_6MEW_B', 'PA_6O21_B~PA_4U2W_B']", "trivial baseline(s) win or tie within Δρ: ['peptide_length']; model_rho=0.1870 baselines={'peptide_length': 0.20460044983108497, 'net_charge': 0.13641841244700895}"], 'notes': "trivial baseline(s) win or tie within Δρ: ['peptide_length']; model_rho=0.1870 baselines={'peptide_length': 0.20460044983108497, 'net_charge': 0.13641841244700895}"}
subset_name: 'peptide_affinity_v2_experimental_openmm'
```

</details>


## Stability / ddG (SKEMPI within-target)

**Status:** FAIL

Within-target held-out SKEMPI ΔΔG: N=16, Spearman ρ=0.48823529411764705 (95% CI [-0.03869047619047619, 0.8295964125560537]). Pre-registered gate ρ≥0.30 with CI_low>0: FAILED. Red-team=pass. Cross-target affinity gate is reported separately and is not replaced.

| Metric | Value | Source |
|---|---|---|
| `skempi_ddg_spearman` | 0.488235 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/skempi/data/skempi_ddg_last_run.json` — N=16; CI=[-0.03869047619047619, 0.8295964125560537] |
| `skempi_ddg_spearman_ci_low` | -0.0386905 | `file:/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/skempi/data/skempi_ddg_last_run.json` |
| `skempi_ddg_gate_threshold` | 0.3 | `ACCEPTANCE.md` |

<details><summary>Details (JSON-backed)</summary>

```
artifact: '/Users/simongriffin/quantum-drug-discovery/quantum-drug-discovery/benchmarks/skempi/data/skempi_ddg_last_run.json'
gate_pass: False
protocol: {'gb_model': 'gbn2', 'solute_dielectric': 1.0, 'salt_conc_M': 0.0}
red_team: {'passed': True, 'label_shuffle_passed': True, 'trivial_baseline_passed': True, 'leakage_passed': True, 'model_rho': 0.48823529411764705, 'baseline_rho': 0.0}
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
