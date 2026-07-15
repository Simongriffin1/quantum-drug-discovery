"""P12 benchmark report tests — traced numbers only, deterministic reload."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from peptideforge.reports.collect import (
    build_benchmark_report,
    default_oracle_artifact,
    load_json_artifact,
)
from peptideforge.reports.render import render_markdown


@pytest.fixture
def oracle_src() -> Path:
    path = default_oracle_artifact()
    if not path.is_file():
        pytest.skip(f"oracle validity artifact missing: {path}")
    return path


@pytest.mark.eval
def test_report_traces_oracle_spearman(oracle_src: Path, tmp_path: Path) -> None:
    art = tmp_path / "artifacts"
    report = build_benchmark_report(
        artifacts_dir=art,
        oracle_artifact=oracle_src,
        seed=0,
    )
    oracle = next(s for s in report.sections if s.title.startswith("Oracle"))
    assert oracle.status == "FAIL"  # known last run
    spearman = next(n for n in oracle.numbers if n.name == "spearman")
    assert spearman.value is not None
    assert spearman.value == pytest.approx(-0.3, abs=0.05)
    assert "file:" in spearman.source
    assert "mlflow_run_id=" in spearman.source or oracle.details.get("mlflow")


@pytest.mark.eval
def test_report_regenerates_deterministically_from_artifacts(
    oracle_src: Path, tmp_path: Path
) -> None:
    art = tmp_path / "artifacts"
    r1 = build_benchmark_report(artifacts_dir=art, oracle_artifact=oracle_src, seed=0)
    r2 = build_benchmark_report(artifacts_dir=art, oracle_artifact=oracle_src, seed=0)
    # Reload path: second build must reuse JSON artifacts — numeric values equal
    def metric_map(report):  # type: ignore[no-untyped-def]
        out: dict[str, float | None] = {}
        for s in report.sections:
            for n in s.numbers:
                out[f"{s.title}:{n.name}"] = n.value
        return out

    assert metric_map(r1) == metric_map(r2)
    md = render_markdown(r1)
    assert "PeptideForge Benchmark Report" in md
    assert "Oracle validity" in md
    assert "Spearman" in md or "spearman" in md
    assert "none are invented" in md.lower() or "Never invent" in md or "invented" in md


def test_missing_artifact_fails_loud(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="benchmark artifact missing"):
        load_json_artifact(tmp_path / "nope.json")


def test_cli_writes_markdown(oracle_src: Path, tmp_path: Path) -> None:
    from peptideforge.reports.generate_benchmark_report import main

    out = tmp_path / "report.md"
    art = tmp_path / "art"
    # seed artifacts once
    build_benchmark_report(artifacts_dir=art, oracle_artifact=oracle_src, seed=1)
    rc = main(
        [
            "--out",
            str(out),
            "--artifacts-dir",
            str(art),
            "--oracle-artifact",
            str(oracle_src),
            "--seed",
            "1",
            "--also-json",
        ]
    )
    assert rc == 0
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert "Gate summary" in text
    assert out.with_suffix(".json").is_file()


def test_copies_oracle_into_sandbox(oracle_src: Path, tmp_path: Path) -> None:
    """Ensure report still works if oracle artifact is relocated."""
    dest = tmp_path / "oracle_validity_last_run.json"
    shutil.copy(oracle_src, dest)
    report = build_benchmark_report(
        artifacts_dir=tmp_path / "a",
        oracle_artifact=dest,
        seed=0,
    )
    assert report.sections[0].details["artifact"] == str(dest)
