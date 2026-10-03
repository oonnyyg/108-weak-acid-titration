"""CSV 曲线与 JSON 汇总输出。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .batch import ExperimentResult
from .solver import KW, TEMPERATURE_C

CSV_COLUMNS = [
    "experiment",
    "base_volume_ml",
    "total_volume_ml",
    "ph",
    "h_mol_per_l",
    "oh_mol_per_l",
    "na_mol_per_l",
    "a_minus_mol_per_l",
    "ha_mol_per_l",
    "ca_diluted_mol_per_l",
    "cb_diluted_mol_per_l",
    "charge_residual_mol_per_l",
    "mass_residual_mol_per_l",
    "solver_iterations",
]


def write_csv(results: list[ExperimentResult], path: str | Path) -> None:
    """写出全部成功实验的滴定曲线 CSV，保留输入顺序与重复点。"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_COLUMNS)
        for result in results:
            if not result.ok:
                continue
            for point in result.points:
                writer.writerow(
                    [
                        result.name,
                        f"{point.base_volume_ml:.6g}",
                        f"{point.total_volume_ml:.6g}",
                        f"{point.ph:.6f}",
                        f"{point.h:.6e}",
                        f"{point.oh:.6e}",
                        f"{point.na:.6e}",
                        f"{point.a_minus:.6e}",
                        f"{point.ha:.6e}",
                        f"{point.ca_diluted:.6e}",
                        f"{point.cb_diluted:.6e}",
                        f"{point.charge_residual:.3e}",
                        f"{point.mass_residual:.3e}",
                        point.iterations,
                    ]
                )


def build_summary(results: list[ExperimentResult]) -> dict:
    """构造 JSON 汇总结构。"""
    experiments = []
    for result in results:
        entry: dict = {"name": result.name}
        if result.ok:
            max_charge = max(
                (abs(p.charge_residual) for p in result.points), default=0.0
            )
            max_mass = max(
                (abs(p.mass_residual) for p in result.points), default=0.0
            )
            entry.update(
                {
                    "status": "ok",
                    "equivalence_volume_ml": result.equivalence_volume_ml,
                    "point_count": len(result.points),
                    "max_abs_charge_residual_mol_per_l": max_charge,
                    "max_abs_mass_residual_mol_per_l": max_mass,
                    "max_solver_iterations": max(
                        (p.iterations for p in result.points), default=0
                    ),
                }
            )
        else:
            entry.update({"status": "error", "error": result.error})
        experiments.append(entry)

    return {
        "model": {
            "description": "一元弱酸-强碱滴定，完整平衡求解（非分段近似）",
            "temperature_c": TEMPERATURE_C,
            "kw": KW,
            "assumptions": [
                "强碱完全解离为一价阳离子",
                "忽略活度系数与体积收缩",
                "浓度 mol/L，体积 mL",
            ],
            "solver": "对数域二分法，相对精度 1e-12，上限 200 次迭代",
        },
        "experiments": experiments,
    }


def write_json(results: list[ExperimentResult], path: str | Path) -> None:
    summary = build_summary(results)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
