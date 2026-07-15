"""Render a BenchmarkReport to deterministic-ish Markdown (stable section order)."""

from __future__ import annotations

from peptideforge.reports.models import BenchmarkReport, TraceableNumber


def _fmt_num(n: TraceableNumber) -> str:
    if n.value is None:
        val = "—"
    elif abs(n.value - round(n.value)) < 1e-12 and abs(n.value) < 1e9:
        val = str(int(round(n.value)))
    else:
        val = f"{n.value:.6g}"
    unit = f" {n.unit}" if n.unit else ""
    notes = f" — {n.notes}" if n.notes else ""
    dv = f", data_version=`{n.data_version}`" if n.data_version else ""
    return f"| `{n.name}` | {val}{unit} | `{n.source}`{dv}{notes} |"


def render_markdown(report: BenchmarkReport) -> str:
    lines: list[str] = [
        f"# {report.title}",
        "",
        f"_Generated at {report.generated_at}_"
        + (f" · git `{report.git_sha}`" if report.git_sha else " · git `unknown`"),
        "",
        "## Caveats",
        "",
    ]
    for c in report.caveats:
        lines.append(f"- {c}")
    lines.extend(["", "## Gate summary", "", "| Section | Status |", "|---|---|"])
    for s in report.sections:
        lines.append(f"| {s.title} | **{s.status}** |")

    for section in report.sections:
        lines.extend(
            [
                "",
                f"## {section.title}",
                "",
                f"**Status:** {section.status}",
                "",
                section.summary,
                "",
            ]
        )
        if section.numbers:
            lines.extend(
                [
                    "| Metric | Value | Source |",
                    "|---|---|---|",
                ]
            )
            for n in section.numbers:
                lines.append(_fmt_num(n))
            lines.append("")
        if section.details:
            lines.append("<details><summary>Details (JSON-backed)</summary>")
            lines.append("")
            lines.append("```")
            for key in sorted(section.details.keys()):
                lines.append(f"{key}: {section.details[key]!r}")
            lines.append("```")
            lines.append("")
            lines.append("</details>")
            lines.append("")

    lines.extend(
        [
            "---",
            "",
            "Regenerate with:",
            "",
            "```bash",
            "make benchmark-report",
            "# or: cd core && poetry run python -m peptideforge.reports.generate_benchmark_report",
            "```",
            "",
        ]
    )
    return "\n".join(lines)
