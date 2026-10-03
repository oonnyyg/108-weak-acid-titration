"""Experiment configuration, validation, and batch titration computation."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .solver import PointResult, solve_point


class ConfigError(ValueError):
    """Raised when one experiment configuration is invalid."""


def _finite_number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{name} must be a number, got {value!r}")
    result = float(value)
    if not math.isfinite(result):
        raise ConfigError(f"{name} must be finite, got {value!r}")
    return result


def _positive(value, name: str) -> float:
    result = _finite_number(value, name)
    if result <= 0.0:
        raise ConfigError(f"{name} must be positive, got {result}")
    return result


@dataclass(frozen=True)
class Experiment:
    name: str
    acid_volume_ml: float
    acid_concentration: float
    ka: float
    base_concentration: float
    additions_ml: tuple[float, ...] = field(default_factory=tuple)

    @property
    def equivalence_volume_ml(self) -> float:
        return self.acid_volume_ml * self.acid_concentration / self.base_concentration


def parse_experiment(raw: dict, index: int) -> Experiment:
    if not isinstance(raw, dict):
        raise ConfigError(f"experiment #{index} must be an object")
    label = raw.get("name", f"experiment-{index}")
    try:
        acid_volume = _positive(raw["acid_volume_ml"], "acid_volume_ml")
        acid_conc = _positive(raw["acid_concentration"], "acid_concentration")
        ka = _positive(raw["ka"], "ka")
        base_conc = _positive(raw["base_concentration"], "base_concentration")
        additions_raw = raw["additions_ml"]
    except KeyError as exc:
        raise ConfigError(f"missing required field {exc.args[0]!r}") from exc
    if not isinstance(additions_raw, list):
        raise ConfigError("additions_ml must be a list")
    additions = []
    for value in additions_raw:
        volume = _finite_number(value, "additions_ml entry")
        if volume < 0.0:
            raise ConfigError(f"additions_ml entry must be >= 0, got {volume}")
        additions.append(volume)
    return Experiment(
        name=str(label),
        acid_volume_ml=acid_volume,
        acid_concentration=acid_conc,
        ka=ka,
        base_concentration=base_conc,
        additions_ml=tuple(additions),
    )


@dataclass(frozen=True)
class ExperimentResult:
    experiment: Experiment
    points: tuple[PointResult, ...]


@dataclass(frozen=True)
class ExperimentFailure:
    name: str
    reason: str


def run_experiment(experiment: Experiment) -> ExperimentResult:
    """Compute the titration curve, preserving input order and duplicates."""
    points = []
    for added in experiment.additions_ml:
        total_volume = experiment.acid_volume_ml + added
        acid_total = (experiment.acid_concentration
                      * experiment.acid_volume_ml / total_volume)
        base_cation = experiment.base_concentration * added / total_volume
        points.append(solve_point(
            acid_total=acid_total,
            base_cation=base_cation,
            ka=experiment.ka,
            volume_added_ml=added,
            total_volume_ml=total_volume,
        ))
    return ExperimentResult(experiment=experiment, points=tuple(points))


def run_batch(raw_experiments: list) -> tuple[list[ExperimentResult], list[ExperimentFailure]]:
    """Run all experiments; one bad config never blocks the others."""
    if not isinstance(raw_experiments, list):
        raise ConfigError("top-level JSON value must be a list of experiments")
    results: list[ExperimentResult] = []
    failures: list[ExperimentFailure] = []
    for index, raw in enumerate(raw_experiments):
        label = raw.get("name", f"experiment-{index}") if isinstance(raw, dict) else f"experiment-{index}"
        try:
            experiment = parse_experiment(raw, index)
            results.append(run_experiment(experiment))
        except (ConfigError, ArithmeticError) as exc:
            failures.append(ExperimentFailure(name=str(label), reason=str(exc)))
    return results, failures
