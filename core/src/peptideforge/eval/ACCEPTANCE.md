# PeptideForge evaluation acceptance thresholds (P2)

These thresholds are **pre-registered before** running real physics oracles (P3+)
or surrogates (P7+). Do not relax them after inspecting results.

Source of truth for stage gates: `CURSOR_PROJECT_CONTEXT.md` §8 / §9.

---

## Oracle-validity gate (affinity)

| Metric | Dataset | Threshold | Notes |
|---|---|---|---|
| Spearman ρ (predicted vs experimental pK) | `pdbbind_peptide_affinity_v1_structures_openmm` (trimmed public MHC–peptide interfaces in `benchmarks/fixtures/structures/`) | **ρ ≥ 0.40** | Pre-registered for M2; report exact subset + N |
| RMSE (same scale as predictions) | Same subset | Reported; no hard fail in M2 | Useful diagnostic |

Failure of the Spearman threshold **blocks** downstream campaigns — every surrogate
and loop result is meaningless without a real oracle.

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
