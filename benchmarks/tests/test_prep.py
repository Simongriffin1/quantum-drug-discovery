"""Offline prep smoke test on a tiny fixture complex."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
INTERFACE = REPO / "fixtures" / "structures" / "PP011_2CLR_interface.pdb"


@pytest.mark.skipif(not INTERFACE.is_file(), reason="interface PDB fixture missing")
def test_prep_runs_on_fixture_and_logs_failures(tmp_path: Path) -> None:
    pytest.importorskip("openmm")
    pytest.importorskip("pdbfixer")
    from peptide_affinity.prep import prepare_complex, write_prep_manifest

    result = prepare_complex(
        INTERFACE,
        out_dir=tmp_path,
        record_id="PA_CI_PREP",
        peptide_chain="C",
        receptor_chains=("A",),
        use_pdb2pqr=False,
        trim_cutoff_nm=2.0,
    )
    assert result.ok, result.error
    assert result.complex_path and Path(result.complex_path).is_file()
    assert result.receptor_path and Path(result.receptor_path).is_file()
    assert result.peptide_path and Path(result.peptide_path).is_file()

    # Failure path: bogus chain → logged exclusion, not silent patch
    bad = prepare_complex(
        INTERFACE,
        out_dir=tmp_path / "bad",
        record_id="PA_CI_BAD",
        peptide_chain="Z",
        receptor_chains=("A",),
    )
    assert bad.ok is False
    assert bad.error
    write_prep_manifest([result, bad], tmp_path / "prep_manifest.tsv")
    text = (tmp_path / "prep_manifest.tsv").read_text(encoding="utf-8")
    assert "PA_CI_BAD" in text and "0" in text
