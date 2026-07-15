"""Benchmark report generation — P12.

Regenerate with ``python -m peptideforge.reports.generate_benchmark_report``.
"""

from peptideforge.reports.collect import build_benchmark_report, default_oracle_artifact
from peptideforge.reports.models import BenchmarkReport, ReportSection, TraceableNumber
from peptideforge.reports.render import render_markdown

__all__ = [
    "BenchmarkReport",
    "ReportSection",
    "TraceableNumber",
    "build_benchmark_report",
    "default_oracle_artifact",
    "render_markdown",
]
