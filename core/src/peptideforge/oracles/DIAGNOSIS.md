# Oracle diagnosis (Step 3A) — train/dev only

**Firewall:** held-out affinity TEST was not used for this diagnosis.

## Attrition audit (3A.2)

Funnel: PepBench 1,433 → matched catalog 68 → prep OK 45 → OpenMM-scoreable **40**
(~2.8% retention). The dominant loss is **no RCSB exact peptide-chain match**
(1,433→68), not prep/OpenMM.

Survivor mean net_charge (**0.33**) is close to PepBench (**0.46**) — no large
charge selection skew. **Call: `oracle_genuinely_weak_vs_charge_baseline`**
(selection bias does not explain the test-set net_charge win by itself).

Artifact: `benchmarks/peptide_affinity/data/attrition_audit.json`.

## Energy decomposition (3A.1)

Method: **OpenMM force-group + charge/LJ scaling** (not AmberTools MMPBSA.py).
Partition: train+val, **N=9**. Correlations of **−Δ** components vs experimental pKd:

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

## One-paragraph diagnosis

On train/dev, **van der Waals / nonpolar terms carry the weak affinity signal**;
**Coulomb and GB polar solvation do not track pKd** (ρ≈0). That is the smoking
gun for an electrostatics/solvation balance bug: the physics term that should
encode the charge-driven signal seen on the held-out test (net_charge ρ≈0.45)
is currently noise in the total ΔG. Therefore Phase 3B must lead with
**protonation rigor + ε_in / salt / GB-model electrostatics**, not a large
force-field×sampling grid. CIs are wide at N=9 — **data recovery (3A.3) is
mandatory** before trusting any tuned “winner.”
