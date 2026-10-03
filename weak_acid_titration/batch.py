"""实验配置读取、校验与批量计算。"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from .solver import (
    SolverError,
    TitrationPointResult,
    equivalence_volume_ml,
    solve_point,
)


class ConfigError(ValueError):
    """单个实验配置非法。"""


@dataclass
class ExperimentConfig:
    name: str
    acid_volume_ml: float
    acid_concentration: float
    ka: float
    base_concentration: float
    base_volumes_ml: list[float]


@dataclass
class ExperimentResult:
    name: str
    equivalence_volume_ml: float = 0.0
    points: list[TitrationPointResult] = field(default_factory=list)
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


def _require_finite(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{label} 必须是数字，得到 {value!r}")
    result = float(value)
    if not math.isfinite(result):
        raise ConfigError(f"{label} 必须是有限值，得到 {value!r}")
    return result


def _require_positive(value: object, label: str) -> float:
    result = _require_finite(value, label)
    if result <= 0.0:
        raise ConfigError(f"{label} 必须为正，得到 {result}")
    return result


def parse_experiment(raw: object, index: int) -> ExperimentConfig:
    """解析并校验单个实验配置。"""
    if not isinstance(raw, dict):
        raise ConfigError(f"第 {index} 个实验必须是 JSON 对象")

    name = raw.get("name", f"experiment_{index}")
    if not isinstance(name, str) or not name:
        raise ConfigError(f"第 {index} 个实验的 name 必须是非空字符串")

    acid_volume = _require_positive(raw.get("acid_volume_ml"), "acid_volume_ml")
    acid_conc = _require_positive(raw.get("acid_concentration"), "acid_concentration")
    ka = _require_positive(raw.get("ka"), "ka")
    base_conc = _require_positive(raw.get("base_concentration"), "base_concentration")

    volumes_raw = raw.get("base_volumes_ml")
    if not isinstance(volumes_raw, list) or not volumes_raw:
        raise ConfigError("base_volumes_ml 必须是非空数组")
    volumes: list[float] = []
    for i, item in enumerate(volumes_raw):
        volume = _require_finite(item, f"base_volumes_ml[{i}]")
        if volume < 0.0:
            raise ConfigError(f"base_volumes_ml[{i}] 不能为负，得到 {volume}")
        volumes.append(volume)

    return ExperimentConfig(
        name=name,
        acid_volume_ml=acid_volume,
        acid_concentration=acid_conc,
        ka=ka,
        base_concentration=base_conc,
        base_volumes_ml=volumes,
    )


def load_config(path: str | Path) -> list[ExperimentConfig]:
    """读取 JSON 配置文件：{"experiments": [ {...}, ... ]}。"""
    text = Path(path).read_text(encoding="utf-8")
    data = json.loads(text)
    if not isinstance(data, dict) or not isinstance(data.get("experiments"), list):
        raise ConfigError("配置文件顶层必须是包含 experiments 数组的对象")
    return [parse_experiment(item, i) for i, item in enumerate(data["experiments"])]


def run_experiment(config: ExperimentConfig) -> ExperimentResult:
    """计算单个实验的完整滴定曲线，保留输入体积顺序与重复点。"""
    result = ExperimentResult(
        name=config.name,
        equivalence_volume_ml=equivalence_volume_ml(
            config.acid_volume_ml, config.acid_concentration, config.base_concentration
        ),
    )
    try:
        for volume in config.base_volumes_ml:
            result.points.append(
                solve_point(
                    acid_volume_ml=config.acid_volume_ml,
                    acid_concentration=config.acid_concentration,
                    ka=config.ka,
                    base_concentration=config.base_concentration,
                    base_volume_ml=volume,
                )
            )
    except SolverError as exc:
        result.points.clear()
        result.error = f"求解失败：{exc}"
    return result


def run_experiments(configs: list[ExperimentConfig]) -> list[ExperimentResult]:
    """批量计算；单个实验失败不影响其他实验。"""
    return [run_experiment(config) for config in configs]


def run_file(path: str | Path) -> list[ExperimentResult]:
    """从 JSON 文件批量计算，单个实验配置错误不影响其他实验。"""
    text = Path(path).read_text(encoding="utf-8")
    data = json.loads(text)
    if not isinstance(data, dict) or not isinstance(data.get("experiments"), list):
        raise ConfigError("配置文件顶层必须是包含 experiments 数组的对象")

    results: list[ExperimentResult] = []
    for index, item in enumerate(data["experiments"]):
        try:
            config = parse_experiment(item, index)
        except ConfigError as exc:
            name = item.get("name") if isinstance(item, dict) else None
            results.append(
                ExperimentResult(
                    name=name if isinstance(name, str) else f"experiment_{index}",
                    error=f"配置错误：{exc}",
                )
            )
            continue
        results.append(run_experiment(config))
    return results
