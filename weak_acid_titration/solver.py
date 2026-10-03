"""弱酸-强碱混合体系的平衡求解。

模型（25 摄氏度）：
    HA <=> H+ + A-，Ka = [H+][A-]/[HA]
    H2O <=> H+ + OH-，Kw = 1e-14
    强碱（如 NaOH）完全解离，[Na+] = Cb
    电荷平衡: [H+] + [Na+] = [OH-] + [A-]
    物料守恒: [HA] + [A-] = Ca（稀释后的酸总浓度）

对 h = [H+] > 0 构造单调函数并用对数域二分法求解，
从初始点到过量碱区使用完全一致的平衡定义。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

KW = 1.0e-14
TEMPERATURE_C = 25.0

_RTOL = 1e-12
_MAX_ITER = 200


class SolverError(RuntimeError):
    """平衡求解未收敛或输入非法。"""


@dataclass
class TitrationPointResult:
    """单个滴定点的完整平衡结果。"""

    base_volume_ml: float
    total_volume_ml: float
    h: float
    oh: float
    na: float
    a_minus: float
    ha: float
    ph: float
    ca_diluted: float
    cb_diluted: float
    charge_residual: float
    mass_residual: float
    iterations: int = field(default=0)


def _charge_imbalance(h: float, ca: float, cb: float, ka: float) -> float:
    """f(h) = [H+] + [Na+] - [OH-] - [A-]，对 h>0 严格单调递增。"""
    return h + cb - KW / h - ka * ca / (ka + h)


def solve_hydrogen(ca: float, cb: float, ka: float) -> tuple[float, int]:
    """求正氢离子浓度，返回 (h, 迭代次数）。"""
    if not (math.isfinite(ca) and math.isfinite(cb) and math.isfinite(ka)):
        raise SolverError("浓度或解离常数不是有限值")
    if ca < 0.0 or cb < 0.0 or ka <= 0.0:
        raise SolverError("浓度必须非负且 Ka 必须为正")

    lo = 1e-300
    hi = max(ca, cb, math.sqrt(KW), 1e-20) * 1e6
    f_lo = _charge_imbalance(lo, ca, cb, ka)
    f_hi = _charge_imbalance(hi, ca, cb, ka)
    if not (f_lo < 0.0 < f_hi):
        raise SolverError("无法在正区间内夹出氢离子浓度的根")

    log_lo, log_hi = math.log(lo), math.log(hi)
    iterations = 0
    for iterations in range(1, _MAX_ITER + 1):
        log_mid = 0.5 * (log_lo + log_hi)
        h_mid = math.exp(log_mid)
        f_mid = _charge_imbalance(h_mid, ca, cb, ka)
        if f_mid > 0.0:
            log_hi = log_mid
        else:
            log_lo = log_mid
        if (log_hi - log_lo) <= _RTOL * abs(log_mid):
            break
    else:
        raise SolverError(f"二分求解未在 {_MAX_ITER} 次迭代内收敛")

    h = math.exp(0.5 * (log_lo + log_hi))
    return h, iterations


def solve_point(
    acid_volume_ml: float,
    acid_concentration: float,
    ka: float,
    base_concentration: float,
    base_volume_ml: float,
) -> TitrationPointResult:
    """对单个累计加碱体积求解完整平衡点。"""
    total_volume_ml = acid_volume_ml + base_volume_ml
    ca = acid_concentration * acid_volume_ml / total_volume_ml
    cb = base_concentration * base_volume_ml / total_volume_ml

    h, iterations = solve_hydrogen(ca, cb, ka)
    oh = KW / h
    a_minus = ka * ca / (ka + h)
    ha = ca - a_minus
    na = cb
    ph = -math.log10(h)

    charge_residual = (h + na) - (oh + a_minus)
    mass_residual = (ha + a_minus) - ca

    return TitrationPointResult(
        base_volume_ml=base_volume_ml,
        total_volume_ml=total_volume_ml,
        h=h,
        oh=oh,
        na=na,
        a_minus=a_minus,
        ha=ha,
        ph=ph,
        ca_diluted=ca,
        cb_diluted=cb,
        charge_residual=charge_residual,
        mass_residual=mass_residual,
        iterations=iterations,
    )


def equivalence_volume_ml(
    acid_volume_ml: float, acid_concentration: float, base_concentration: float
) -> float:
    """理论当量体积：n(HA) = n(OH-) 时的累计加碱体积（mL）。"""
    return acid_concentration * acid_volume_ml / base_concentration
