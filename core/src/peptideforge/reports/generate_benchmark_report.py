"""CLI: generate PeptideForge credibility benchmark report (P12).

Pulls oracle-validity from logged JSON (+ optional MLflow), regenerates
synthetic_* surrogate / acquisition / loop sections with fixed seeds, and
writes Markdown whose every number cites a source path or run ID.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from peptideforge.reports.collect import build_benchmark_report, default_oracle_artifact
from peptideforge.reports.render import render_markdown


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Markdown output path (default: benchmarks/reports/benchmark_report.md)",
    )
    parser.add_argument(
        "--artifacts-dir",
        type=Path,
        default=None,
        help="Directory for surrogate/loop/acquisition JSON artifacts",
    )
    parser.add_argument(
        "--oracle-artifact",
        type=Path,
        default=None,
        help="Path to oracle_validity_last_run.json",
    )
    parser.add_argument(
        "--mlflow-uri",
        type=str,
        default=None,
        help="Optional MLflow tracking URI to annotate latest oracle run",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--also-json",
        action="store_true",
        help="Also write the BenchmarkReport model as JSON next to the markdown",
    )
    args = parser.parse_args(argv)

    report = build_benchmark_report(
        artifacts_dir=args.artifacts_dir,
        oracle_artifact=args.oracle_artifact or default_oracle_artifact(),
        mlflow_uri=args.mlflow_uri,
        seed=args.seed,
    )
    md = render_markdown(report)

    out = args.out
    if out is None:
        out = (
            Path(__file__).resolve().parents[4]
            / "benchmarks"
            / "reports"
            / "benchmark_report.md"
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md, encoding="utf-8")
    print(f"wrote {out}")

    if args.also_json:
        json_path = out.with_suffix(".json")
        json_path.write_text(report.model_dump_json(indent=2), encoding="utf-8")
        print(f"wrote {json_path}")

    # Exit non-zero only if caller wants CI fail — here we always succeed if report wrote.
    # Gate status is honest inside the document (oracle may be FAIL).
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
