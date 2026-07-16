# PeptideForge evaluation acceptance thresholds (P2)

These thresholds are **pre-registered before** running real physics oracles (P3+)
or surrogates (P7+). Do not relax them after inspecting results.

Source of truth for stage gates: `CURSOR_PROJECT_CONTEXT.md` §8 / §9.

---

## Oracle-validity gate (affinity)

| Metric | Dataset | Threshold | Notes |
|---|---|---|---|
| Spearman ρ (predicted vs experimental pKd) | Held-out cold-start test of `peptide_affinity_v2_experimental_openmm` (experimental prepared structures; PepBench PpI_ba labels ∩ RCSB) | **ρ ≥ 0.40** with bootstrap 95% **CI lower bound > 0** | **N ≥ 40** required before PASS/FAIL is measurable (raised from 30 after Step 2) |
| Pearson r | Same held-out test | Reported with bootstrap 95% CI | Diagnostic |
| RMSE | Same subset | Reported; no hard fail in M2 | Diagnostic |
| Red-team | label_shuffle + trivial baselines (length, charge, contacts, SASA) + leakage (@30% ID) | **All must pass**; oracle must beat **net_charge** | If a trivial baseline wins → halt; task/split flawed |

### Predicted-fold within-target authorization — pre-registered 2026-07-16 (Step 4)

**Committed before any predicted-fold scoring verdict.** Do not loosen after seeing numbers.

Live (predicted-fold) within-target campaigns are **authorized iff** predicted-fold
Spearman ρ ≥ **0.30** with bootstrap 95% **CI_low > 0** on the **same** SKEMPI
homology-aware hold-out used for the experimental PASS (`skempi_ddg_powered_last_run`,
N=100 reference ρ=0.381 on crystals) — either:

1. **Unconditionally** on all predicted folds in that hold-out, **or**
2. **Restricted** to folds above a stated confidence threshold (ipTM / DockQ / interface-pLDDT)
   with **surviving N ≥ 30** and the gate still holding on that stratum.

Otherwise predicted-fold campaigns remain **BLOCKED**, even though experimental-structure
scoring is validated. Cross-target absolute affinity remains **not authorized**.

Primary measurement mode: **mutate_in_place** on a Boltz-2 predicted WT (denovo re-fold is
comparison only). Experimental ρ=0.381 is the reference and must never be overwritten or
re-cited as the predicted number.

| Field | Value |
|---|---|
| Pre-registration commit intent | Step 4A.5 before Phase 4A.3/4A.4 verdict |
| Experimental reference | ρ=0.381, CI [0.189, 0.556], N=100, crystals |
| Predicted gate | ρ ≥ 0.30, CI_low > 0; if stratified, N_surviving ≥ 30 |



Legacy diagnostic subset (N=5 MHC interfaces): `pdbbind_peptide_affinity_v1_structures_openmm`
— statistically uninterpretable alone; kept for continuity.

Failure of both co-primary gates **blocks** unauthorized binding claims. A within-target
PASS with cross-target FAIL may authorize **scoped** within-target campaigns only, with
the cross-target limitation documented. Point estimates without a CI **must not** be
used to declare PASS/FAIL.

## Stability / ddG gate

| Metric | Dataset | Threshold |
|---|---|---|
| Spearman ρ (predicted vs experimental ΔΔG) | Named SKEMPI mutation subset | **ρ ≥ 0.30** |

## Surrogate calibration gate (P7)

| Metric | Threshold |
|---|---|
| Expected Calibration Error (ECE) | **ECE < 0.10** at target coverage (default 0.90) |
| Red-team controls | All must **pass** (see below) |

Implementation: `DeepEnsembleSurrogate` (bootstrap ridge ensemble + split conformal).
Evaluate with `run_surrogate_acceptance` on homology-aware holdouts. Passing the
gate on a `synthetic_*` learnable oracle validates **plumbing + UQ math only** —
it does **not** satisfy oracle-validity (P3) or authorize binding campaigns.

## Acquisition validation (P8)

| Metric | Benchmark | Threshold |
|---|---|---|
| Final hypervolume | Synthetic Branin–Currin (maximize −Branin, −Currin) | **qNEHVI HV > random HV** (same budget) |
| Batch constraints | Cost + `max_*` / `min_*` / `exclude_ids` | Honored; no overspend |

Implementation: stdlib Monte Carlo qNEHVI (BoTorch optional; fails loud if
`use_botorch=True` without install). Validates acquisition **logic**, not biology.

## Loop efficiency gate (P9)

| Metric | Threshold |
|---|---|
| Simulations-to-target vs random | Loop reaches target affinity/hv in **fewer** oracle calls than random (seeded, shared pool) |

Validated in **simulation mode** (`synthetic_physics_label`) on a shared discrete
pool so acquisition — not generator luck — is compared. Public fixture path uses
PDBbind-peptide **sequences** as seeds with the same synthetic oracle (does **not**
claim the P3 affinity gate). On that small pool, ties in simulations-to-target are
accepted as non-inferiority; the strict beat-random gate is the synthetic shared pool.

## Agent attribution gate (P10)

| Rule | Threshold |
|---|---|
| Reported numbers | **Every** numeric claim must appear in a `ToolResult` attribution ledger |
| Stage gates | Oracle-validity / spend pauses must fire; human (or explicit simulation skip) required |
| Hallucination test | Summary with invented floats (e.g. Spearman 0.99) **must fail** ledger validation |

## Platform provenance gate (P11)

| Artifact | Requirement |
|---|---|
| Campaign / Pareto / structure / calibration / trace | Each response includes `provenance` with `data_version` + `tool_versions` (+ `git_sha` when available) |
| Non-simulation | Rejected until oracle-validity gate passes |
| UI panels | Pareto, Mol* (or CA fallback), calibration, agent trace all load from API mocks in tests |

## Benchmark report gate (P12)

| Requirement | Detail |
|---|---|
| Regenerable | `make benchmark-report` / `python -m peptideforge.reports.generate_benchmark_report` |
| Traceability | Every metric cites `file:…` and/or `mlflow_run_id=…` (or `ACCEPTANCE.md` for thresholds) |
| Honesty | Oracle FAIL / ddG NOT_RUN must appear as such — no padded Spearmans |
| Synthetic sections | Surrogate ECE, Branin–Currin HV, loop spend: labeled as plumbing, not physics |



## Quantum gate (P4)

| Metric | Threshold |
|---|---|
| H2 ground-state (VQE vs PySCF classical at R=0.7414 Å, STO-3G FCI) | Within **0.01 hartree** |
| Classical baseline | **Always** returned with every VQE result |
| He atom HF/STO-3G | Within **0.01 hartree** of literature (−2.80774745 Eh) |
| Missing deps | Raise (`PySCFUnavailableError` / `PennyLaneUnavailableError`) — never fabricate |

---

## Red-team controls (always required before believing a result)

1. **label_shuffle_control** — Spearman of predictions vs *shuffled* labels must drop
   below `|ρ| < 0.20` (or lose ≥ 0.40 absolute Spearman vs unshuffled). Flags
   label-memorizing / hardcoded experimental leaks.
2. **trivial_baseline_check** — Model Spearman must beat a trivial baseline
   (mean predictor / sequence-length heuristic) by **Δρ ≥ 0.05**.
3. **leakage_audit** — No train/test cluster overlap; no identical `record_id` or
   peptide sequence shared across train and test. Any leakage → **FAIL**.

Deliberately leaky toy sets used in unit tests must be flagged by these controls.

---

## Fixture note

`benchmarks/fixtures/*_v1.tsv` are **format fixtures** for loader/harness CI.
They do **not** constitute a passed oracle-validity gate.

P3 runs OpenMM MM-GBSA on `benchmarks/fixtures/structures/structure_manifest_v1.tsv`
(`pdbbind_peptide_affinity_v1_structures_openmm`), logs Spearman to MLflow / 
`oracle_validity_last_run.json`, and compares against the thresholds above.
Gate failure is reported honestly — never padded or fabricated.
