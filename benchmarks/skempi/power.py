"""Phase 5A.1 — PDB-level power analysis (cheap, before any folding spend).

Step 4 numbers were mutation-level and therefore anti-conservative when many
mutations share the same PDB/fold (pseudo-replication). This script estimates
how many *unique PDBs* are needed for:

1) A pooled Spearman ρ (predicted vs experimental ΔΔG) to have a clustered
   (PDB-block) bootstrap CI with CI_low > 0, at an effect size near the
   observed experimental value (≈0.38).

2) A paired per-PDB within-target comparison (experimental vs predicted fold)
   to detect a plausible degradation in per-PDB ρ with adequate power.

This is intentionally assumption-driven; it is used to set a curation target
N_pdb* before spending compute.
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from typing import Iterable

from peptideforge.eval.metrics import spearman_rho


@dataclass(frozen=True)
class PowerAssumptions:
    """Assumptions for Monte Carlo power estimation.

    - n_mut_per_pdb: average usable mutations per PDB (after filters)
    - rho: true pooled Spearman correlation between predicted and experimental
    - icc: intra-PDB correlation of residual noise (0 = iid mutations; 1 = identical)
    - n_boot: bootstrap replicates for clustered CI
    """

    n_mut_per_pdb: int = 10
    rho: float = 0.38
    icc: float = 0.50
    n_boot: int = 600
    seed: int = 0


def _rank(x: list[float]) -> list[float]:
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    for k, i in enumerate(order):
        r[i] = float(k)
    return r


def simulate_clustered_dataset(
    *,
    n_pdb: int,
    n_mut_per_pdb: int,
    rho: float,
    icc: float,
    rng: random.Random,
) -> tuple[list[float], list[float], list[int]]:
    """Generate a synthetic clustered (PDB-blocked) dataset.

    We model an experimental latent score x and a predicted latent score y with
    correlation rho at the mutation level, then add a shared per-PDB noise
    component to induce within-PDB dependence controlled by icc.
    """
    if not (0.0 <= icc <= 1.0):
        raise ValueError("icc must be in [0,1]")
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1,1]")
    if n_pdb < 1 or n_mut_per_pdb < 3:
        raise ValueError("need n_pdb>=1 and n_mut_per_pdb>=3")

    exp: list[float] = []
    pred: list[float] = []
    cluster: list[int] = []

    # Shared noise strength chosen so that (shared / total) ≈ icc in expectation.
    # Total noise variance = shared^2 + iid^2. Choose shared = sqrt(icc).
    shared_scale = math.sqrt(icc)
    iid_scale = math.sqrt(max(0.0, 1.0 - icc))

    for pdb_i in range(n_pdb):
        # Per-PDB shared components
        e_shared = rng.gauss(0.0, 1.0)
        p_shared = rng.gauss(0.0, 1.0)
        for _ in range(n_mut_per_pdb):
            # Base correlated normals for mutation-level signal
            z1 = rng.gauss(0.0, 1.0)
            z2 = rng.gauss(0.0, 1.0)
            x0 = z1
            y0 = rho * z1 + math.sqrt(max(0.0, 1.0 - rho * rho)) * z2

            # Add clustered noise
            x = x0 + shared_scale * e_shared + iid_scale * rng.gauss(0.0, 1.0)
            y = y0 + shared_scale * p_shared + iid_scale * rng.gauss(0.0, 1.0)
            exp.append(x)
            pred.append(y)
            cluster.append(pdb_i)

    # Convert to ranks (Spearman acts on ranks)
    return _rank(exp), _rank(pred), cluster


def clustered_bootstrap_ci(
    pred: list[float],
    exp: list[float],
    cluster_ids: list[int],
    *,
    n_boot: int,
    seed: int,
) -> tuple[float, float, float]:
    """Clustered (PDB-block) bootstrap CI for Spearman ρ."""
    if not (len(pred) == len(exp) == len(cluster_ids)):
        raise ValueError("length mismatch")
    ids = sorted(set(cluster_ids))
    by_cluster: dict[int, list[int]] = {c: [] for c in ids}
    for i, c in enumerate(cluster_ids):
        by_cluster[c].append(i)
    if len(ids) < 3:
        raise ValueError("need ≥3 clusters")

    point = spearman_rho(pred, exp)
    rng = random.Random(seed)
    samples: list[float] = []
    for _ in range(n_boot):
        boot_clusters = [ids[rng.randrange(len(ids))] for _ in ids]
        idx: list[int] = []
        for c in boot_clusters:
            idx.extend(by_cluster[c])
        samples.append(spearman_rho([pred[i] for i in idx], [exp[i] for i in idx]))
    samples.sort()
    lo = samples[int(0.025 * len(samples))]
    hi = samples[int(0.975 * len(samples))]
    return point, lo, hi


def estimate_n_pdb_for_ci_low_positive(
    *,
    assumptions: PowerAssumptions,
    n_pdb_grid: Iterable[int],
    n_mc: int = 250,
    require_prob: float = 0.80,
) -> dict[str, float | int | dict[str, float]]:
    """Find smallest N_pdb where P(CI_low>0) ≥ require_prob under the assumptions."""
    rng = random.Random(assumptions.seed)
    for n_pdb in n_pdb_grid:
        ok = 0
        for _ in range(n_mc):
            pred, exp, clusters = simulate_clustered_dataset(
                n_pdb=n_pdb,
                n_mut_per_pdb=assumptions.n_mut_per_pdb,
                rho=assumptions.rho,
                icc=assumptions.icc,
                rng=rng,
            )
            _, lo, _ = clustered_bootstrap_ci(
                pred, exp, clusters, n_boot=assumptions.n_boot, seed=rng.randrange(10**9)
            )
            if lo > 0.0:
                ok += 1
        prob = ok / float(n_mc)
        if prob >= require_prob:
            return {
                "n_pdb_star": n_pdb,
                "p_ci_low_gt_0": prob,
                "assumptions": {
                    "rho": assumptions.rho,
                    "n_mut_per_pdb": assumptions.n_mut_per_pdb,
                    "icc": assumptions.icc,
                    "n_boot": assumptions.n_boot,
                    "n_mc": n_mc,
                    "require_prob": require_prob,
                },
            }
    return {
        "n_pdb_star": -1,
        "p_ci_low_gt_0": 0.0,
        "assumptions": {
            "rho": assumptions.rho,
            "n_mut_per_pdb": assumptions.n_mut_per_pdb,
            "icc": assumptions.icc,
            "n_boot": assumptions.n_boot,
            "n_mc": n_mc,
            "require_prob": require_prob,
        },
    }


def estimate_paired_per_pdb_power(
    *,
    n_pdb: int,
    rho_exp: float,
    rho_pred: float,
    sd_per_pdb: float = 0.35,
    n_mc: int = 2000,
    alpha: float = 0.05,
    seed: int = 0,
) -> dict[str, float | int]:
    """Paired per-PDB power via Fisher-z approx with per-PDB variability.

    We treat each PDB's within-target Spearman as a noisy draw around a mean
    (rho_exp / rho_pred) with sd_per_pdb. We then test mean(z_exp - z_pred) > 0.
    """
    rng = random.Random(seed)
    # Fisher transform guard
    def z(r: float) -> float:
        r = min(0.999, max(-0.999, r))
        return 0.5 * math.log((1 + r) / (1 - r))

    z_exp_mu = z(rho_exp)
    z_pred_mu = z(rho_pred)
    # Approximate z-space sd (use delta method around mu)
    z_sd = sd_per_pdb
    thr = abs(_normal_quantile(1 - alpha / 2))
    sig = 0
    for _ in range(n_mc):
        diffs = [
            (rng.gauss(z_exp_mu, z_sd) - rng.gauss(z_pred_mu, z_sd)) for _ in range(n_pdb)
        ]
        mean = sum(diffs) / n_pdb
        sd = math.sqrt(sum((d - mean) ** 2 for d in diffs) / max(1, n_pdb - 1))
        t = 0.0 if sd == 0 else mean / (sd / math.sqrt(n_pdb))
        if abs(t) > thr and mean > 0:
            sig += 1
    return {"n_pdb": n_pdb, "power": sig / float(n_mc), "alpha": alpha, "n_mc": n_mc}


def _normal_quantile(p: float) -> float:
    """Inverse CDF for standard normal (Acklam approximation)."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0,1)")
    # Coefficients from Peter John Acklam's approximation
    a = [
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    ]
    b = [
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    ]
    plow = 0.02425
    phigh = 1 - plow
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        )
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        )
    q = p - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        / (
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho", type=float, default=0.38)
    parser.add_argument("--mut-per-pdb", type=int, default=10)
    parser.add_argument("--icc", type=float, default=0.5)
    parser.add_argument("--require-prob", type=float, default=0.8)
    parser.add_argument("--n-mc", type=int, default=250)
    parser.add_argument("--n-boot", type=int, default=600)
    parser.add_argument("--grid-max", type=int, default=60)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--rho-pred-plausible", type=float, default=0.28)
    args = parser.parse_args()

    assumptions = PowerAssumptions(
        n_mut_per_pdb=args.mut_per_pdb,
        rho=args.rho,
        icc=args.icc,
        n_boot=args.n_boot,
        seed=args.seed,
    )
    grid = list(range(5, args.grid_max + 1, 1))
    pooled = estimate_n_pdb_for_ci_low_positive(
        assumptions=assumptions, n_pdb_grid=grid, n_mc=args.n_mc, require_prob=args.require_prob
    )

    # Paired per-PDB plausibility check: detect degradation from rho to rho_pred_plausible
    paired = estimate_paired_per_pdb_power(
        n_pdb=max(5, int(pooled["n_pdb_star"]) if int(pooled["n_pdb_star"]) > 0 else 30),
        rho_exp=args.rho,
        rho_pred=args.rho_pred_plausible,
        seed=args.seed,
    )

    print("=== Phase 5A.1 power analysis (assumption-driven) ===")
    print(
        f"Assumptions: rho≈{assumptions.rho:.2f}, mut/pdb≈{assumptions.n_mut_per_pdb}, "
        f"icc≈{assumptions.icc:.2f}, clustered_boot={assumptions.n_boot}, mc={args.n_mc}, "
        f"require P(CI_low>0)≥{args.require_prob:.2f}"
    )
    n_star = int(pooled["n_pdb_star"])
    if n_star < 0:
        print(f"N_pdb* not reached up to grid_max={args.grid_max}. Increase --grid-max or revisit assumptions.")
    else:
        print(f"N_pdb* (pooled clustered CI_low>0): {n_star} unique PDBs (estimated).")
        print(f"  achieved P(CI_low>0)≈{float(pooled['p_ci_low_gt_0']):.2f} at N_pdb={n_star}")
    print("")
    print("Paired per-PDB degradation (Fisher-z approx; diagnostic):")
    print(
        f"  n_pdb={paired['n_pdb']}, test mean(z_exp - z_pred)>0 with alpha={paired['alpha']}, "
        f"power≈{paired['power']:.2f} for rho_exp≈{args.rho:.2f} vs rho_pred≈{args.rho_pred_plausible:.2f}"
    )


if __name__ == "__main__":
    main()

