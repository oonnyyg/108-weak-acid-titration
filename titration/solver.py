"""Weak acid / strong base titration equilibrium solver (25 C, Kw = 1e-14)."""

from __future__ import annotations

import math
from dataclasses import dataclass

KW = 1.0e-14
H_MIN = 1.0e-14
H_MAX = 1.0
SOLVER_TOL = 1.0e-12
SOLVER_MAX_ITER = 200


@dataclass(frozen=True)
class PointResult:
    volume_added_ml: float
    total_volume_ml: float
    acid_total: float
    base_cation: float
    h: float
    oh: float
    ha: float
    a_minus: float
    ph: float
    charge_residual: float
    mass_residual: float
    mass_residual: float


def solve_point(acid_total: float, base_cation: float, ka: float,
                volume_added_ml: float, total_volume_ml: float) -> PointResult:
    """Solve charge balance h + Cb = Kw/h + Ka*Ca/(Ka+h) for h > 0.

    The residual is monotonic in h, so bisection on geometric midpoints
    converges reliably from the initial point through the equivalence
    region into excess base, with one consistent equilibrium definition.
    """
    lo, hi = H_MIN, H_MAX

    def residual(h: float) -> float:
        return h + base_cation - KW / h - ka * acid_total / (ka + h)

    if residual(lo) >= 0.0:
        raise ArithmeticError("lower bracket residual is non-negative")
    converged = False
    for _ in range(SOLVER_MAX_ITER):
        mid = math.sqrt(lo * hi)
        if residual(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi / lo - 1.0 < SOLVER_TOL:
            converged = True
            break
    if not converged:
        raise ArithmeticError("solver did not converge within iteration limit")

    h = math.sqrt(lo * hi)
    oh = KW / h
    a_minus = ka * acid_total / (ka + h)
    ha = acid_total - a_minus
    charge_residual = (h + base_cation) - (oh + a_minus)
    mass_residual = acid_total - (ha + a_minus)
    return PointResult(
        volume_added_ml=volume_added_ml,
        total_volume_ml=total_volume_ml,
        acid_total=acid_total,
        base_cation=base_cation,
        h=h,
        oh=oh,
        ha=ha,
        a_minus=a_minus,
        ph=-math.log10(h),
        charge_residual=charge_residual,
        mass_residual=mass_residual,
    )
