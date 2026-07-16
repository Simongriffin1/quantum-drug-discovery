# Oracle diagnosis (Step 3A–3B) — train/dev only

**Firewall:** held-out affinity TEST was not used for diagnosis or protocol tuning.

## Attrition audit (3A.2)

Funnel (v2): PepBench 1,433 → matched catalog 68 → prep OK 45 → OpenMM-scoreable **40**
(~2.8% retention). Dominant loss = **no RCSB exact peptide-chain match**.

Survivor mean net_charge (**0.33**) ≈ PepBench (**0.46**) — no large charge selection
skew. **Call: `oracle_genuinely_weak_vs_charge_baseline`** (selection bias does not
explain the test-set net_charge win by itself).

Artifact: `benchmarks/peptide_affinity/data/attrition_audit.json`.

## Energy decomposition (3A.1)

Method: **OpenMM force-group + charge/LJ scaling** (not AmberTools MMPBSA.py).
Partition: train+val, **N=9** (v2). Correlations of **−Δ** components vs experimental pKd:

| Term | Spearman ρ | 95% CI |
|---|---|---|
| −total | 0.30 | [−0.53, 0.93] |
| −vdW | 0.32 | [−0.63, 0.89] |
| −elec (Coulomb) | **−0.03** | [−0.71, 0.88] |
| −polar_solv (GB) | **0.00** | [−1.00, 0.66] |
| −(elec+polar) | 0.10 | [−0.69, 0.95] |
| −nonpolar | 0.37 | [−0.43, 0.90] |
| net_charge | −0.05 | [−0.81, 0.68] |

Artifact: `benchmarks/peptide_affinity/data/decomposition_traindev.json`.

## Data recovery (3A.3)

Catalog v3: PepBench ∩ RCSB relaxed matching → **113** catalog rows; prep recoveries
(trim 3.0 nm) + OpenMM scoreable filter → **69** scoreable complexes.
Cold-start 30% ID split (`splits_v3.json`): **train=23, val=6, test=40** (53 clusters).
Targets met: test N≥40, train/dev N=29 (≥25). Propedia/PepBDB full dumps were **not**
auto-downloaded (README marks manual) — refuse to fabricate; PepBench expansion used instead.

## Electrostatics sweep (3B.2) — v2 N=9 result

OpenMM amber14 GB XML ignores `soluteDielectric` / `implicitSolventSaltConc`; ε_in is
implemented via **charge scaling** (documented in `mm_gbsa.py`). On v2 train/dev N=9:

- Best: `gbsa_gbn2_eps1_salt0_min0` ρ≈**0.367**
- Raising ε_in ∈ {2,4,8} and salt 0.15 M **hurt** ρ

v3 re-sweep on N≈29 is logged to `elec_sweep_traindev_v3.json` (test firewalled).

## One-paragraph diagnosis

On train/dev, **van der Waals / nonpolar terms carry the weak affinity signal**;
**Coulomb and GB polar solvation do not track pKd** (ρ≈0). That localizes a broken
electrostatic/solvation balance relative to the held-out net_charge baseline (ρ≈0.45),
but damping ε_in via charge scaling did **not** recover the signal on N=9 — so either
the charge–pKd association is not encoded in the end-point Coulomb+GB terms for these
poses, or protonation/structure prep still corrupts them. Phase 3C (entropy, interface
waters) and the pre-registered within-target SKEMPI gate proceed with the best
electrostatics config found on enlarged train/dev; the cross-target gate remains the
stated hard limit.

## Electrostatics sweep (3B.2) — v3 N=29

Artifact: `elec_sweep_traindev_v3.json`. **Best: `gbsa_gbn2_eps1_salt0_min0` ρ≈−0.016**
(CI crosses 0). Raising ε_in / salt **monotonically hurt** ρ (down to ≈−0.28).
Electrostatics tuning did **not** move train/dev toward the gate. Chosen protocol for
one-shot test = this best-of-a-bad-grid config (GBn2, ε_in=1, salt=0, no minimize).

## Entropy + interface waters (3C)

Artifact: `entropy_waters_ab_traindev.json`. Truncated-NMA entropy ≈ noise (ρ 0.018 vs
baseline −0.016). Interface waters: many OpenMM HOH template failures; scored subset
n=12 with **worse** ρ. Flags off for the chosen protocol.

## SKEMPI within-target (3D.2) — pre-registered, parallel

Artifact: `benchmarks/skempi/data/skempi_ddg_last_run.json`. Held-out N=16, Spearman
**ρ=0.488**, 95% CI **[−0.039, 0.830]** — point estimate above 0.30 but **CI_low ≤ 0 →
FAIL** under pre-registered gate. Red-team passed. Mutation = side-chain strip + PDBFixer
(not rename-only).


## One-shot cross-target re-test (3D.3) — firewalled test, once

Artifact: `oracle_validity_v3_oneshot_test.json` (MLflow `aeefc1bdfd1d4e289b923912c8639001`).
Protocol locked from train/dev: `gbsa_gbn2_eps1_salt0_min0`.

| Gate | Result |
|---|---|
| Cross-target affinity | **FAIL** — N=40, ρ=0.187, CI [−0.160, 0.514] (CI_low ≤ 0; ρ < 0.40) |
| Within-target SKEMPI | **FAIL** — N=16, ρ=0.488, CI [−0.039, 0.830] (ρ≥0.30 but CI_low ≤ 0) |

Trivial baselines on test: peptide_length ρ≈0.205 beats oracle; net_charge ρ≈0.136
(no longer dominates as on the v2 N=31 set). Red-team also flags sequence-identity
leakage across the v3 cold-start split — split hygiene must be tightened before any
future re-test; thresholds were **not** loosened.

**Verdict:** Both co-primary gates fail after diagnosis-led electrostatics tuning,
entropy/water A/B, and data recovery. This is a **scientific finding**: single-structure
end-point MM-GBSA is insufficient for this peptide affinity target class as currently
posed. Escalation options (document only): short-MD ensembles, explicit-solvent MM-PBSA,
or a component re-scorer — **not** on ~40 points.

## Structure provenance (clarification)

All affinity oracle scores used **experimental RCSB crystal structures**
(`fold_method=experimental_prepared`, PepBench labels ∩ deposited PDBs).
**Not** Boltz / AlphaFold predicted complexes. Pose error is therefore **not**
confounded with oracle error on this fork — the cross-target failure (even if
optimistic under leakage) is an end-point physics / prep problem on experimental
poses.

## Leakage stop-the-line (splits_v3 → v4)

`splits_v3` clustered on **receptor only**; red-team found peptide-sequence
leakage (identical + >30% ID). **All affinity ρ from the v3 one-shot are
INVALIDATED** (optimistic direction). `splits_v4` uses joint receptor∪peptide
30% identity clustering; leakage_audit passes. N_test=37 (<40) — integrity over
sample size; do not interpret a new affinity ρ until/unless re-scored once on v4
with the locked protocol (no retuning).
