"""弱酸滴定工具的单元测试（标准库 unittest）。"""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from weak_acid_titration.batch import (
    ConfigError,
    ExperimentConfig,
    parse_experiment,
    run_experiment,
    run_file,
)
from weak_acid_titration.report import build_summary, write_csv
from weak_acid_titration.solver import KW, solve_hydrogen, solve_point

# 醋酸体系：25 mL 0.1 M，Ka = 1.8e-5，用 0.1 M NaOH 滴定
ACID_V = 25.0
ACID_C = 0.1
KA = 1.8e-5
BASE_C = 0.1


def make_config(volumes: list[float]) -> ExperimentConfig:
    return ExperimentConfig(
        name="acetic",
        acid_volume_ml=ACID_V,
        acid_concentration=ACID_C,
        ka=KA,
        base_concentration=BASE_C,
        base_volumes_ml=volumes,
    )


class SolverTest(unittest.TestCase):
    def test_zero_base_matches_weak_acid_formula(self):
        """零加碱时 [H+] 与弱酸解析解一致。"""
        point = solve_point(ACID_V, ACID_C, KA, BASE_C, 0.0)
        expected = (-KA + math.sqrt(KA * KA + 4 * KA * ACID_C)) / 2
        self.assertAlmostEqual(point.h, expected, delta=expected * 1e-6)
        self.assertAlmostEqual(point.ph, -math.log10(expected), places=5)
        self.assertAlmostEqual(point.na, 0.0)

    def test_half_equivalence_buffer(self):
        """半当量点 pH 接近 pKa。"""
        point = solve_point(ACID_V, ACID_C, KA, BASE_C, 12.5)
        self.assertAlmostEqual(point.ph, -math.log10(KA), delta=1e-3)

    def test_equivalence_point_is_basic(self):
        """当量点弱碱性，[OH-] 与 A- 水解解析解一致。"""
        point = solve_point(ACID_V, ACID_C, KA, BASE_C, 25.0)
        ca = ACID_C * ACID_V / (ACID_V + 25.0)
        kb = KW / KA
        oh = (-kb + math.sqrt(kb * kb + 4 * kb * ca)) / 2
        self.assertAlmostEqual(point.oh, oh, delta=oh * 1e-3)
        self.assertGreater(point.ph, 7.0)
        self.assertLess(point.ph, 9.0)

    def test_excess_base_region(self):
        """过量碱区 [OH-] 接近过量强碱浓度。"""
        point = solve_point(ACID_V, ACID_C, KA, BASE_C, 50.0)
        excess = (BASE_C * 50.0 - ACID_C * ACID_V) / (ACID_V + 50.0)
        self.assertAlmostEqual(point.oh, excess, delta=excess * 1e-3)
        self.assertGreater(point.ph, 11.0)

    def test_conservation_residuals_small(self):
        """每个点的电荷与物料守恒残差可忽略。"""
        for volume in [0.0, 5.0, 12.5, 24.9, 25.0, 25.1, 40.0]:
            point = solve_point(ACID_V, ACID_C, KA, BASE_C, volume)
            self.assertLess(abs(point.charge_residual), 1e-12)
            self.assertLess(abs(point.mass_residual), 1e-12)

    def test_ph_monotonic_increasing(self):
        """pH 随累计加碱体积单调上升。"""
        volumes = [float(i) for i in range(0, 51)]
        phs = [solve_point(ACID_V, ACID_C, KA, BASE_C, v).ph for v in volumes]
        for earlier, later in zip(phs, phs[1:]):
            self.assertLess(earlier, later)

    def test_solver_rejects_invalid_input(self):
        with self.assertRaises(Exception):
            solve_hydrogen(-1.0, 0.0, KA)
        with self.assertRaises(Exception):
            solve_hydrogen(0.1, 0.0, 0.0)
        with self.assertRaises(Exception):
            solve_hydrogen(float("nan"), 0.0, KA)


class BatchTest(unittest.TestCase):
    def test_order_and_duplicates_preserved(self):
        """保留输入体积顺序与重复点。"""
        volumes = [10.0, 0.0, 10.0, 25.0]
        result = run_experiment(make_config(volumes))
        self.assertTrue(result.ok)
        self.assertEqual([p.base_volume_ml for p in result.points], volumes)
        self.assertEqual(
            result.points[0].ph, result.points[2].ph
        )

    def test_equivalence_volume(self):
        result = run_experiment(make_config([0.0]))
        self.assertAlmostEqual(result.equivalence_volume_ml, 25.0)

    def test_config_validation(self):
        base = {
            "acid_volume_ml": 25.0,
            "acid_concentration": 0.1,
            "ka": KA,
            "base_concentration": 0.1,
            "base_volumes_ml": [0.0],
        }
        bad_cases = [
            {"acid_volume_ml": 0.0},
            {"acid_volume_ml": -1.0},
            {"acid_concentration": -0.1},
            {"ka": 0.0},
            {"ka": float("inf")} if False else {"base_concentration": "x"},
            {"base_volumes_ml": [-1.0]},
            {"base_volumes_ml": []},
            {"base_volumes_ml": [float("nan")]},
        ]
        for override in bad_cases:
            raw = {**base, **override}
            with self.assertRaises(ConfigError, msg=str(override)):
                parse_experiment(raw, 0)

    def test_run_file_fault_tolerance(self):
        """单个实验配置错误不影响其他实验。"""
        payload = {
            "experiments": [
                {
                    "name": "good",
                    "acid_volume_ml": 25.0,
                    "acid_concentration": 0.1,
                    "ka": KA,
                    "base_concentration": 0.1,
                    "base_volumes_ml": [0.0, 25.0],
                },
                {
                    "name": "bad",
                    "acid_volume_ml": -25.0,
                    "acid_concentration": 0.1,
                    "ka": KA,
                    "base_concentration": 0.1,
                    "base_volumes_ml": [0.0],
                },
            ]
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            results = run_file(path)
        self.assertEqual(len(results), 2)
        self.assertTrue(results[0].ok)
        self.assertFalse(results[1].ok)
        self.assertIn("配置错误", results[1].error)


class ReportTest(unittest.TestCase):
    def test_csv_and_summary(self):
        results = [run_experiment(make_config([0.0, 12.5, 25.0, 50.0]))]
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp) / "curve.csv"
            write_csv(results, csv_path)
            lines = csv_path.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(len(lines), 1 + 4)
        self.assertIn("charge_residual_mol_per_l", lines[0])

        summary = build_summary(results)
        entry = summary["experiments"][0]
        self.assertEqual(entry["status"], "ok")
        self.assertAlmostEqual(entry["equivalence_volume_ml"], 25.0)
        self.assertLess(entry["max_abs_charge_residual_mol_per_l"], 1e-12)
        self.assertEqual(summary["model"]["kw"], KW)


if __name__ == "__main__":
    unittest.main()
