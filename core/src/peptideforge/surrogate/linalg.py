"""Linear algebra helpers for ridge regression (stdlib only)."""

from __future__ import annotations

from collections.abc import Sequence


def mat_vec(a: list[list[float]], x: Sequence[float]) -> list[float]:
    return [sum(aij * xj for aij, xj in zip(row, x, strict=True)) for row in a]


def transpose(a: list[list[float]]) -> list[list[float]]:
    if not a:
        return []
    return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]


def mat_mul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    bt = transpose(b)
    return [[sum(ai * bj for ai, bj in zip(row, col, strict=True)) for col in bt] for row in a]


def eye(n: int, scale: float = 1.0) -> list[list[float]]:
    return [[scale if i == j else 0.0 for j in range(n)] for i in range(n)]


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[aij + bij for aij, bij in zip(ra, rb, strict=True)] for ra, rb in zip(a, b, strict=True)]


def solve_linear(a: list[list[float]], b: Sequence[float]) -> list[float]:
    """Gaussian elimination with partial pivoting for Ax = b."""
    n = len(a)
    if n == 0:
        return []
    # Augmented matrix
    m = [row[:] + [bj] for row, bj in zip(a, b, strict=True)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise ValueError("singular matrix in ridge solve")
        m[col], m[pivot] = m[pivot], m[col]
        piv = m[col][col]
        for j in range(col, n + 1):
            m[col][j] /= piv
        for row in range(n):
            if row == col:
                continue
            factor = m[row][col]
            for j in range(col, n + 1):
                m[row][j] -= factor * m[col][j]
    return [m[i][n] for i in range(n)]


def fit_ridge(
    x: list[list[float]],
    y: Sequence[float],
    *,
    l2: float = 1.0,
) -> list[float]:
    """Closed-form ridge: (XᵀX + λI)⁻¹ Xᵀy. Assumes bias column already in X."""
    if len(x) != len(y):
        raise ValueError("x and y length mismatch")
    if not x:
        raise ValueError("empty design matrix")
    xt = transpose(x)
    xtx = mat_mul(xt, x)
    dim = len(xtx)
    xtx_reg = mat_add(xtx, eye(dim, l2))
    xty = mat_vec(xt, y)
    return solve_linear(xtx_reg, xty)


def predict_linear(weights: Sequence[float], features: Sequence[float]) -> float:
    return sum(w * f for w, f in zip(weights, features, strict=True))
