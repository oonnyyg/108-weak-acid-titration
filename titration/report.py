"""CSV and JSON report generation for titration results."""

from __future__ import annotations

import csv
import io
import math

from .experiment import ExperimentFailure, ExperimentResult
from .solver import KW, SOLVER_TOL

CSV_COLUMNS = (
    "experiment",
    "volume_added_ml",
    "total_volume_ml",
    "acid_total_mol_L",
    "base_cation_mol_L",
    "h_mol_L",
    "oh_mol_L",
    "ha_mol_L",
    "a_minus_mol_L",
    "ph",
    "charge_residual_mol_L",
    "mass_residual_mol_L",
)


def render_csv(results: list[ExperimentResult]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(CSV_COLUMNS)
    for result in results:
        for point in result.points:
            writer.writerow((
                result.experiment.name,
                f"{point.volume_added_ml:.6g}",
                f"{point.total_volume_ml:.6g}",
                f"{point.acid_total:.6e}",
                f"{point.base_cation:.6e}",
                f"{point.h:.6e}",
                f"{point.oh:.6e}",
                f"{point.ha:.6e}",
                f"{point.a_minus:.6e}",
                f"{point.ph:.6f}",
                f"{point.charge_residual:.3e}",
                f"{point.mass_residual:.3e}",
            ))
    return buffer.getvalue()
    return buffer.getvalue()


def build_summary(results: list[ExperimentResult],
                  failures: list[ExperimentFailure]) -> dict:
    experiments = []
    for result in results:
        exp = result.experiment
        points = []
        for point in result.points:
            points.append({
                "volume_added_ml": point.volume_added_ml,
                "total_volume_ml": point.total_volume_ml,
                "concentrations_mol_L": {
                    "H+": point.h,
                    "OH-": point.oh,
                    "HA": point.ha,
                    "A-": point.a_minus,
                    "base_cation": point.base_cation,
                },
                "ph": point.ph,
                "charge_residual_mol_L": point.charge_residual,
                "mass_residual_mol_L": point.mass_residual,
            })
        experiments.append({
            "name": exp.name,
            "acid_volume_ml": exp.acid_volume_ml,
            "acid_concentration_mol_L": exp.acid_concentration,
            "ka": exp.ka,
            "pka": _pka(exp.ka),
            "base_concentration_mol_L": exp.base_concentration,
            "equivalence_volume_ml": exp.equivalence_volume_ml,
            "points": points,
        })
    return {
        "model": {
            "temperature_C": 25.0,
            "kw": KW,
            "solver": "bisection on charge balance, geometric midpoint",
            "solver_relative_tolerance": SOLVER_TOL,
            "assumptions": [
                "strong base fully dissociates to a monovalent cation",
                "activities equal concentrations",
                "volumes are additive (no contraction)",
            ],
        },
        "experiments": experiments,
        "failures": [{"name": f.name, "reason": f.reason} for f in failures],
    }


def _pka(ka: float) -> float:
    return -math.log10(ka)
