"""Command line entry point for the weak acid titration tool."""

from __future__ import annotations

import argparse
import json
import sys

from .experiment import ConfigError, run_batch
from .report import build_summary, render_csv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="titration",
        description="Weak acid / strong base titration curve calculator (25 C).",
    )
    parser.add_argument("input", help="JSON file with a list of experiment configs")
    parser.add_argument("--csv", dest="csv_path", help="write the CSV curve to this path")
    parser.add_argument("--json", dest="json_path", help="write the JSON summary to this path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        with open(args.input, "r", encoding="utf-8") as handle:
            raw = json.load(handle)
    except OSError as exc:
        print(f"error: cannot read input file: {exc}", file=sys.stderr)
        return 2
    except json.JSONDecodeError as exc:
        print(f"error: invalid JSON: {exc}", file=sys.stderr)
        return 2

    try:
        results, failures = run_batch(raw)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    for failure in failures:
        print(f"warning: experiment {failure.name!r} skipped: {failure.reason}",
              file=sys.stderr)

    csv_text = render_csv(results)
    summary_text = json.dumps(build_summary(results, failures),
                              ensure_ascii=False, indent=2)

    if args.csv_path:
        with open(args.csv_path, "w", encoding="utf-8", newline="") as handle:
            handle.write(csv_text)
    else:
        sys.stdout.write(csv_text)
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as handle:
            handle.write(summary_text + "\n")
    else:
        sys.stdout.write(summary_text + "\n")
    return 1 if failures and not results else 0


if __name__ == "__main__":
    raise SystemExit(main())
