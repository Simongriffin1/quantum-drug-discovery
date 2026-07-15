# PeptideForge

Physics-grounded, agentic, closed-loop peptide design platform.

**Oracle = physics** (OpenMM MM-GBSA/MD, PySCF). No proprietary training data.
Validated against public benchmarks (PDBbind, SKEMPI). Humans own stage gates.

Read [`CURSOR_PROJECT_CONTEXT.md`](./CURSOR_PROJECT_CONTEXT.md) before contributing.

## Layout

| Path | Role |
|---|---|
| `core/` | `peptideforge` Python package (contracts first) |
| `backend/` | FastAPI |
| `frontend/` | Next.js + TypeScript + Tailwind |
| `infra/` | docker-compose (Postgres+pgvector, MinIO, MLflow) |
| `benchmarks/` | Public benchmark loaders (P2) |
| `docs/` | Documentation |

## Quick start

```bash
# Infra
make docker-up

# Python 3.11 + Poetry required
pipx install poetry==1.8.3   # if needed
make install
make schemas
make test
make typecheck
```

**Status:** M1–M5 scaffolding complete (**P0–P12**). Expanded oracle-validity
(Step 2) on experimental PepBench∩RCSB structures is **measurable and FAIL**:
held-out N=31, Spearman **0.12** (95% CI −0.27–0.47), net_charge trivial
baseline wins — binding campaigns remain unauthorized. Artifacts:
`benchmarks/peptide_affinity/data/oracle_validity_v2_last_run.json` and
`benchmarks/reports/benchmark_report.md`. Protocol sweep (Step 3) is in
progress on train/dev only.
Acceptance for benchmark report (P12): `make benchmark-report` regenerates
Markdown/JSON from oracle JSON (+ optional MLflow) and seeded synthetic_*
artifacts; every metric cites a file path or run ID.

Acceptance for platform (P11): POST `/campaigns` runs simulation DBTL + agent;
GET pareto / structures / calibration / trace return real artifacts with
provenance; frontend workspace mounts all four panels.

Acceptance for agent (P10): `PeptideForgeAgent` runs synthetic-mode end-to-end;
every reported number traces to a `ToolResult`; stage-gate pauses fire; invented
floats are rejected by the attribution ledger.

Acceptance for closed loop (P9): `ClosedLoopOrchestrator` persists `LoopState`;
simulation mode runs without heavy deps; spend gate passes on shared pool;
reports attribute numbers to tools.

Acceptance for acquisition (P8): `QNEHVIAcquisition` / `RandomAcquisition` honor
`AcquisitionFunction`; batches respect budget + constraints; qNEHVI dominates
random by hypervolume on synthetic Branin–Currin.

Acceptance for surrogate UQ (P7): `DeepEnsembleSurrogate` honors `Surrogate`;
empirical ECE vs pre-registered threshold; red-team controls; JSON calibration /
acceptance report artifact.

Acceptance for generation/folding (P6): `MutationGenerator` / `PeptideGenerator`
return valid diverse `Candidates`; `FixtureStructurePredictor` returns cached
interface PDB + confidence or raises; `Boltz2StructurePredictor` raises when CLI
missing (no fabricated coordinates).

Acceptance for developability (P5): five algorithmic predictors honor
`DevelopabilityPredictor`; `PeptideDevelopabilityEvaluator` returns per-property
(value, uncertainty) vector without scalar collapse; NetMHC not bundled (license).

Acceptance for scaffold (P1): compose stack runs; mypy clean; contract tests
validate/reject shapes.

Acceptance for eval harness (P2): loaders parse fixtures (no network in CI);
harness computes Spearman/RMSE; red-team flags deliberately leaked toy sets.

Acceptance for physics oracle (P3): contract + cost cap; energy-conservation golden;
real Spearman logged to MLflow.

## Non-negotiables

No fabricated results. No hallucinated APIs. Fail loud. Cost caps. Spec/eval first.
See §8 in `CURSOR_PROJECT_CONTEXT.md`.
