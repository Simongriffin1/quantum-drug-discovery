# PeptideForge — Cursor Project Context

> Reference this file in every prompt. Non-negotiable rules (§8) apply to all work.

## What we are building

We are building **PeptideForge**: a physics-grounded, agentic, closed-loop peptide design
platform. Physics simulation (OpenMM MM-GBSA/MD, PySCF quantum chemistry) is the **ORACLE** —
there is **NO proprietary training data**; ground truth is computed from first principles and
validated against **PUBLIC benchmarks** (PDBbind, SKEMPI). Modality: therapeutic peptides
(5–50 residues).

**Thesis:** Traditional AI-drug-discovery tools win by owning proprietary training sets. We
don't have one — so we don't build that kind of tool. Binding/stability/energy labels come
from docking, MD, MM-GBSA, and quantum chemistry (first principles). That makes physics an
**infinite, self-labeling data generator**. Cheap surrogate models train on physics output;
expensive simulation only labels; the surrogate does mass screening. Credibility comes from
showing the physics oracle **reproduces known experimental affinities on public benchmarks**.

## The loop

```
GENERATE (ESM sampling)
  → FOLD (Boltz-2 peptide–target complex)
  → PHYSICS ORACLE (tiered: docking → MM-GBSA → MD → FEP → quantum chem)
  → DEVELOPABILITY (aggregation, solubility, MHC-II immunogenicity, synthesizability)
  → SURROGATE + calibrated UQ (deep ensemble + conformal)
  → multi-objective Bayesian acquisition (BoTorch qNEHVI)
  → back to FOLD
```

An LLM agent orchestrates and explains; **humans own stage gates**. The versioned
(sequence, structure, physics-label) dataset is the core asset.

## Stack

| Layer | Tool |
|---|---|
| Language / core | Python 3.11, strict typing (mypy) |
| Embeddings / generation | ESM-2 (ESM-3 only if license verified) |
| Structure / complex fold | Boltz-2 (or Chai-1); avoid AF3 weights commercially |
| MD / MM-GBSA | OpenMM |
| Quantum chemistry | PySCF; PennyLane/Qiskit VQE (always vs PySCF baseline) |
| Surrogate / UQ / acquisition | PyTorch, GPyTorch, BoTorch |
| Orchestration | Ray, Prefect |
| Agent | LLM API + tool-calling harness |
| Data | PostgreSQL + pgvector, DVC/LakeFS, MinIO/S3 |
| Backend | FastAPI, Pydantic |
| Frontend | Next.js, TypeScript (strict), Tailwind, Mol\*, Plotly |
| Experiment tracking | MLflow |
| Infra | Docker, docker-compose, cloud GPU |

**License rule of thumb:** prefer MIT/Apache/BSD (OpenMM, PySCF, Vina, Boltz, ProteinMPNN,
ESM-2, BoTorch). Verify all licenses at build time. Keep every third-party model **pluggable**.
Avoid/verify: AlphaFold3 weights, Rosetta, FoldX, ESM-3 community license, NetMHC commercial
terms. Data licenses: check PDBbind redistribution; prefer permissive sets (e.g. SKEMPI/PDB)
where possible.

## Monorepo layout

```
/core        — peptideforge Python package (contracts, oracles, surrogates, loop, eval, reports)
/backend     — FastAPI
/frontend    — Next.js + TS + Tailwind
/infra       — Docker, docker-compose (Postgres+pgvector, MinIO, MLflow)
/benchmarks  — public benchmark loaders
/docs        — documentation
```

## Style

- Strict typing: mypy (Python), TypeScript `strict`
- pytest for all core logic
- Docstrings explain biological/physical rationale, not just what the code does
- Spec-first, eval-first: contracts and acceptance criteria BEFORE implementation
- Seed everything; pin versions; log git SHA + data version + config to MLflow

---

## NON-NEGOTIABLE RULES (§8 — Development procedure)

Put these rules in every prompt. Do not proceed past stage gates without human sign-off.

### Hard rules

- **No fabricated or hardcoded results.** Synthetic data only when named `synthetic_*` and
  used solely for shape/plumbing tests. The agent reports only real tool outputs.
- **No hallucinated APIs.** Verify imports; fail loudly if a function/tool is unavailable —
  never silently fall back or invent numbers.
- **Physics golden tests are mandatory.** Validate each physics component against a known
  reference before trusting it.
- **Homology-aware splits** for any surrogate benchmark; **red-team** every claimed result
  (label-shuffle, trivial-baseline, leakage audit, ablation) before believing it.
- **Fail loud** (explicit exceptions), **seed everything**, **pin versions**, **log git SHA
  + data version + config** to MLflow.
- **Cost caps** on every job that uses GPU/quantum; default to the cheapest oracle tier.
- **Spec-first, eval-first:** write contracts and the evaluation/acceptance criteria BEFORE
  the implementation.

### Human stage gates (do not proceed without sign-off)

1. **Oracle-validity gate:** physics oracle reproduces known affinities on a public
   benchmark (Spearman above a pre-set threshold) — otherwise every downstream number is
   meaningless.
2. **Calibration gate:** surrogate intervals are calibrated (ECE below threshold) before
   they drive acquisition.
3. **Loop gate:** the full loop beats random search on a held-out public benchmark before
   any expensive campaign.
4. **Quantum gate:** VQE reproduces the H2 reference energy; every quantum result ships
   with its PySCF classical baseline.
5. **Spend gate:** confirm the loop reduces simulations-to-target vs. random before scaling
   compute.

### Golden reference tests (§9)

- **MM-GBSA / oracle validity:** on a public protein–peptide affinity subset (PDBbind /
  SKEMPI), predicted binding energies correlate with experiment (report Spearman + exact
  subset).
- **Stability (ddG):** on SKEMPI mutation data, predicted ddG correlates with measured ddG.
- **Quantum chemistry:** PySCF/VQE reproduces H2 ground-state energy at a known bond length
  within tolerance; small-molecule DFT matches a literature reference.
- **MD sanity:** energy conservation / equilibration checks on a known system.

---

## Build order (do not skip gates)

| Phase | Prompt | Milestone |
|---|---|---|
| P0 | This context file | — |
| P1 | Scaffold + contracts | M1 |
| P2 | Public benchmark loaders + oracle-validity harness | M1 |
| P3 | Physics oracle: MM-GBSA via OpenMM | M2 (oracle-validity gate) |
| P4 | Quantum chemistry tier | M2 (quantum gate) |
| P5 | Developability predictors | M3 |
| P6 | Generation + folding | M3 |
| P7 | Surrogate + calibrated UQ | M3 (calibration gate) |
| P8 | Acquisition (qNEHVI) | M3 |
| P9 | Closed-loop orchestrator | M3 (loop gate / spend gate) |
| P10 | LLM agent orchestrator | M4 |
| P11 | Platform, 3D viz, provenance | M4 |
| P12 | Benchmark report | M5 |

## Milestones

- **M1 (foundation):** P0–P2 — scaffold, contracts, public benchmarks, eval harness.
- **M2 (real oracle):** P3–P4 — MM-GBSA + quantum tier; **oracle-validity gate passed**.
- **M3 (smart loop, cheap):** P5–P9 in simulation/benchmark mode; loop **beats random**.
- **M4 (agentic + usable):** P10–P11 — agent + demoable web UI.
- **M5 (credible):** P12 — reproducible benchmark report / preprint basis.
