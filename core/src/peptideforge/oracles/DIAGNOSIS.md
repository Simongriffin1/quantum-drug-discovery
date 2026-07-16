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
