"""命令行入口：python -m weak_acid_titration 配置.json --csv 曲线.csv --json 汇总.json"""

from __future__ import annotations

import argparse
import sys

from .batch import ConfigError, run_file
from .report import write_csv, write_json


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="weak_acid_titration",
        description="一元弱酸-强碱滴定曲线计算（25 摄氏度，完整平衡求解）",
    )
    parser.add_argument("config", help="实验配置 JSON 文件路径")
    parser.add_argument("--csv", dest="csv_path", help="滴定曲线 CSV 输出路径")
    parser.add_argument("--json", dest="json_path", help="汇总 JSON 输出路径")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        results = run_file(args.config)
    except (ConfigError, OSError, ValueError) as exc:
        print(f"错误：无法读取配置：{exc}", file=sys.stderr)
        return 2

    for result in results:
        if result.ok:
            print(
                f"[完成] {result.name}: {len(result.points)} 个点，"
                f"当量体积 {result.equivalence_volume_ml:.4f} mL"
            )
        else:
            print(f"[失败] {result.name}: {result.error}", file=sys.stderr)

    if args.csv_path:
        write_csv(results, args.csv_path)
        print(f"CSV 曲线已写入 {args.csv_path}")
    if args.json_path:
        write_json(results, args.json_path)
        print(f"JSON 汇总已写入 {args.json_path}")

    return 1 if any(not r.ok for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
