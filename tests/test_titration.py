import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from titration.experiment import ConfigError, parse_experiment, run_batch, run_experiment
from titration.solver import KW, solve_point

KA = 1.8e-5  # acetic acid


def make_config(additions, **overrides):
    config = {
        "name": "acetic",
        "acid_volume_ml": 50.0,
        "acid_concentration": 0.1,
        "ka": KA,
        "base_concentration": 0.1,
        "additions_ml": additions,
    }
    config.update(overrides)
    return config


class SolverTests(unittest.TestCase):
    def test_zero_addition_matches_weak_acid_approximation(self):
        ca = 0.1
        point = solve_point(ca, 0.0, KA, 0.0, 50.0)
        h_approx = (-KA + math.sqrt(KA**2 + 4 * KA * ca)) / 2
        self.assertAlmostEqual(point.h, h_approx, delta=1e-9)
        self.assertAlmostEqual(point.ph, -math.log10(h_approx), places=6)

    def test_half_equivalence_ph_near_pka(self):
        result = run_experiment(parse_experiment(make_config([25.0]), 0))
        point = result.points[0]
        self.assertAlmostEqual(point.ph, -math.log10(KA), delta=0.02)

    def test_equivalence_point_is_basic(self):
        result = run_experiment(parse_experiment(make_config([50.0]), 0))
        point = result.points[0]
        self.assertGreater(point.ph, 7.0)
        self.assertLess(point.ph, 9.5)

    def test_excess_base_approaches_strong_base_limit(self):
        result = run_experiment(parse_experiment(make_config([100.0]), 0))
        point = result.points[0]
        oh_excess = 0.1 * 50.0 / 150.0
        self.assertAlmostEqual(point.oh, oh_excess, delta=1e-4)
        self.assertGreater(point.ph, 11.0)

    def test_residuals_are_tiny_across_curve(self):
        additions = [0.0, 5.0, 25.0, 49.9, 50.0, 50.1, 75.0, 150.0]
        result = run_experiment(parse_experiment(make_config(additions), 0))
        for point in result.points:
            self.assertLess(abs(point.charge_residual), 1e-10)
            self.assertLess(abs(point.mass_residual), 1e-10)
            self.assertAlmostEqual(point.h * point.oh, KW, delta=1e-16)
            self.assertAlmostEqual(point.h * point.oh, KW, delta=1e-16)


class ValidationTests(unittest.TestCase):
    def test_rejects_nonpositive_values(self):
        for field in ("acid_volume_ml", "acid_concentration", "ka", "base_concentration"):
            with self.assertRaises(ConfigError):
                parse_experiment(make_config([0.0], **{field: 0.0}), 0)

    def test_rejects_nonfinite_and_negative_addition(self):
        with self.assertRaises(ConfigError):
            parse_experiment(make_config([float("nan")]), 0)
        with self.assertRaises(ConfigError):
            parse_experiment(make_config([float("inf")]), 0)
        with self.assertRaises(ConfigError):
            parse_experiment(make_config([-1.0]), 0)

    def test_zero_addition_allowed(self):
        exp = parse_experiment(make_config([0.0]), 0)
        self.assertEqual(exp.additions_ml, (0.0,))


class BatchTests(unittest.TestCase):
    def test_order_and_duplicates_preserved(self):
        configs = [make_config([10.0, 0.0, 10.0, 10.0])]
        results, failures = run_batch(configs)
        self.assertEqual(failures, [])
        volumes = [p.volume_added_ml for p in results[0].points]
        self.assertEqual(volumes, [10.0, 0.0, 10.0, 10.0])

    def test_bad_experiment_does_not_block_others(self):
        configs = [
            make_config([0.0], name="good"),
            make_config([0.0], ka=-1.0, name="bad"),
            make_config([50.0], name="also-good"),
        ]
        results, failures = run_batch(configs)
        self.assertEqual([r.experiment.name for r in results], ["good", "also-good"])
        self.assertEqual(len(failures), 1)
        self.assertEqual(failures[0].name, "bad")
        self.assertIn("ka", failures[0].reason)


class CliTests(unittest.TestCase):
    def test_cli_writes_csv_and_json(self):
        configs = [make_config([0.0, 25.0, 50.0, 75.0])]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            input_path = tmp_path / "in.json"
            csv_path = tmp_path / "curve.csv"
            json_path = tmp_path / "summary.json"
            input_path.write_text(json.dumps(configs), encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, "-m", "titration", str(input_path),
                 "--csv", str(csv_path), "--json", str(json_path)],
                capture_output=True, text=True,
                cwd=Path(__file__).resolve().parents[1],
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            rows = csv_path.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(len(rows), 5)  # header + 4 points
            summary = json.loads(json_path.read_text(encoding="utf-8"))
            exp = summary["experiments"][0]
            self.assertAlmostEqual(exp["equivalence_volume_ml"], 50.0)
            self.assertEqual(len(exp["points"]), 4)
            self.assertEqual(summary["failures"], [])


if __name__ == "__main__":
    unittest.main()
