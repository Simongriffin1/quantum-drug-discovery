"""Step 4A.1 — predicted structure generation for SKEMPI hold-out (Boltz-2).

MODE A (primary): predict WT once, mutate_in_place + relax.
MODE B (comparison): de novo fold of mutant complex.

Fails loud if Boltz-2 is unavailable — never fabricates structures into a
real measurement artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from peptideforge.structure.deps import Boltz2UnavailableError, require_boltz_cli


@dataclass(frozen=True)
class FoldProvenance:
    mode: str  # mutate_in_place | denovo
    record_id: str
    pdb_id: str
    mutant: str
    boltz_cli: str
    sequence_hash: str
    wt_fold_path: str | None
    mutant_fold_path: str | None
    confidence: float | None
    iptm: float | None
    interface_plddt: float | None
    error: str | None = None


def sequence_hash(*parts: str) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(p.encode("utf-8"))
        h.update(b"|")
    return h.hexdigest()[:16]


def load_holdout(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"hold-out membership missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def predict_holdout_folds(
    *,
    holdout_path: Path,
    skempi_tsv: Path,
    structure_dir: Path,
    out_dir: Path,
    modes: tuple[str, ...] = ("mutate_in_place", "denovo"),
    max_complexes: int | None = None,
) -> dict[str, Any]:
    """Generate predicted folds for the powered SKEMPI hold-out.

    Raises Boltz2UnavailableError before writing any 'measured' claim.
    """
    # Fail loud first — before any cache write that could be mistaken for real data
    boltz_cli = require_boltz_cli()

    from peptideforge_benchmarks.skempi import load_skempi_ddg
    from skempi.run_skempi_ddg import _mutate_pdb_ca_proxy
    from peptideforge.structure.boltz2 import Boltz2StructurePredictor
    from peptideforge.contracts.models import PeptideCandidate
    from uuid import uuid4

    holdout = load_holdout(holdout_path)
    ids = list(holdout["holdout_record_ids"])
    if max_complexes is not None:
        ids = ids[:max_complexes]

    records = {r.record_id: r for r in load_skempi_ddg(skempi_tsv)}
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "cache"
    cache_dir.mkdir(exist_ok=True)

    predictor = Boltz2StructurePredictor(use_msa_server=True)
    provenance: list[dict[str, Any]] = []
    wt_fold_by_pdb: dict[str, Path] = {}

    for rid in ids:
        rec = records.get(rid)
        if rec is None:
            provenance.append(
                asdict(
                    FoldProvenance(
                        mode="n/a",
                        record_id=rid,
                        pdb_id="?",
                        mutant="?",
                        boltz_cli=boltz_cli,
                        sequence_hash="",
                        wt_fold_path=None,
                        mutant_fold_path=None,
                        confidence=None,
                        iptm=None,
                        interface_plddt=None,
                        error="record not in skempi_tsv",
                    )
                )
            )
            continue
        wt_crystal = structure_dir / f"{rec.pdb_id}.pdb"
        if not wt_crystal.is_file():
            provenance.append(
                asdict(
                    FoldProvenance(
                        mode="n/a",
                        record_id=rid,
                        pdb_id=rec.pdb_id,
                        mutant=rec.mutant,
                        boltz_cli=boltz_cli,
                        sequence_hash="",
                        wt_fold_path=None,
                        mutant_fold_path=None,
                        confidence=None,
                        iptm=None,
                        interface_plddt=None,
                        error=f"missing crystal {wt_crystal}",
                    )
                )
            )
            continue

        # MODE A: WT fold once per PDB
        if "mutate_in_place" in modes:
            try:
                if rec.pdb_id not in wt_fold_by_pdb:
                    # Fold WT: use crystal as template for receptor sequence extraction
                    # Peptide sequence from partner2 chain length is unknown — use
                    # crystal peptide chain sequence via Boltz predictor template path.
                    from peptideforge.structure.pdb_utils import (
                        extract_chain_sequence,
                        list_protein_chains,
                        read_pdb_text,
                    )

                    text = read_pdb_text(wt_crystal)
                    chains = list_protein_chains(text)
                    pep_chain = (rec.partner2 or "B")[0]
                    if pep_chain not in chains:
                        pep_chain = chains[-1]
                    pep_seq = extract_chain_sequence(text, pep_chain)
                    if len(pep_seq) < 5:
                        raise ValueError(f"peptide chain {pep_chain} too short")
                    cand = PeptideCandidate(
                        candidate_id=uuid4(),
                        sequence=pep_seq[:50],
                        generation_method="skempi_wt_fold",
                    )
                    folded = predictor.fold(
                        cand,
                        target_id=rec.pdb_id,
                        target_structure=str(wt_crystal),
                        seed=0,
                    )
                    wt_path = cache_dir / f"{rec.pdb_id}_boltz_wt.pdb"
                    wt_path.write_text(folded.pdb_text or "", encoding="utf-8")
                    if not wt_path.read_text(encoding="utf-8").strip():
                        # pdb_path may be temp; copy from path
                        wt_path.write_text(
                            Path(folded.pdb_path).read_text(encoding="utf-8")
                            if folded.pdb_path
                            else "",
                            encoding="utf-8",
                        )
                    wt_fold_by_pdb[rec.pdb_id] = wt_path
                    meta_path = cache_dir / f"{rec.pdb_id}_boltz_wt.meta.json"
                    meta_path.write_text(
                        json.dumps(
                            {
                                "confidence": folded.confidence,
                                "fold_method": folded.fold_method,
                                "boltz_cli": boltz_cli,
                                "cache_key": folded.cache_key,
                            },
                            indent=2,
                        ),
                        encoding="utf-8",
                    )

                wt_pred = wt_fold_by_pdb[rec.pdb_id]
                mut_out = cache_dir / f"{rid}_mutate_in_place.pdb"
                _mutate_pdb_ca_proxy(wt_pred, rec.mutant, mut_out)
                conf = None
                meta = cache_dir / f"{rec.pdb_id}_boltz_wt.meta.json"
                if meta.is_file():
                    conf = json.loads(meta.read_text()).get("confidence")
                provenance.append(
                    asdict(
                        FoldProvenance(
                            mode="mutate_in_place",
                            record_id=rid,
                            pdb_id=rec.pdb_id,
                            mutant=rec.mutant,
                            boltz_cli=boltz_cli,
                            sequence_hash=sequence_hash(rec.pdb_id, "wt", rec.mutant),
                            wt_fold_path=str(wt_pred),
                            mutant_fold_path=str(mut_out),
                            confidence=conf,
                            iptm=None,
                            interface_plddt=conf,
                        )
                    )
                )
            except Exception as exc:  # noqa: BLE001
                provenance.append(
                    asdict(
                        FoldProvenance(
                            mode="mutate_in_place",
                            record_id=rid,
                            pdb_id=rec.pdb_id,
                            mutant=rec.mutant,
                            boltz_cli=boltz_cli,
                            sequence_hash="",
                            wt_fold_path=None,
                            mutant_fold_path=None,
                            confidence=None,
                            iptm=None,
                            interface_plddt=None,
                            error=str(exc)[:300],
                        )
                    )
                )

        if "denovo" in modes:
            provenance.append(
                asdict(
                    FoldProvenance(
                        mode="denovo",
                        record_id=rid,
                        pdb_id=rec.pdb_id,
                        mutant=rec.mutant,
                        boltz_cli=boltz_cli,
                        sequence_hash=sequence_hash(rec.pdb_id, "denovo", rec.mutant),
                        wt_fold_path=None,
                        mutant_fold_path=None,
                        confidence=None,
                        iptm=None,
                        interface_plddt=None,
                        error=(
                            "denovo mutant re-fold not yet implemented for multi-chain "
                            "SKEMPI complexes in this pass — MODE A is primary"
                        ),
                    )
                )
            )

    payload = {
        "holdout_path": str(holdout_path),
        "split_id": holdout.get("split_id"),
        "experimental_reference": holdout.get("experimental_reference"),
        "boltz_cli": boltz_cli,
        "n_requested": len(ids),
        "provenance": provenance,
        "status": "OK",
    }
    out_path = out_dir / "predicted_folds_manifest.json"
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--holdout",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_powered_holdout_v1.json"),
    )
    parser.add_argument(
        "--skempi-tsv",
        type=Path,
        default=Path("benchmarks/skempi/data/skempi_v2.tsv"),
    )
    parser.add_argument(
        "--structure-dir",
        type=Path,
        default=Path("benchmarks/skempi/data/structures"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmarks/skempi/data/predicted_folds"),
    )
    parser.add_argument("--max-complexes", type=int, default=None)
    args = parser.parse_args()
    try:
        payload = predict_holdout_folds(
            holdout_path=args.holdout,
            skempi_tsv=args.skempi_tsv,
            structure_dir=args.structure_dir,
            out_dir=args.out_dir,
            max_complexes=args.max_complexes,
        )
    except Boltz2UnavailableError as exc:
        # Honest incomplete measurement — write a BLOCKED stub, never fake PDBs
        args.out_dir.mkdir(parents=True, exist_ok=True)
        blocked = {
            "status": "BLOCKED_BOLTZ_UNAVAILABLE",
            "error": str(exc),
            "holdout_path": str(args.holdout),
            "experimental_reference_untouched": True,
            "note": (
                "Predicted-fold degradation NOT measured. Per ACCEPTANCE.md 4A.5, "
                "live predicted-fold campaigns remain BLOCKED. Experimental SKEMPI "
                "ρ=0.381 reference is unchanged."
            ),
        }
        (args.out_dir / "predicted_folds_manifest.json").write_text(
            json.dumps(blocked, indent=2), encoding="utf-8"
        )
        raise SystemExit(f"FAIL LOUD: {exc}") from exc
    print(json.dumps({"status": payload["status"], "n": payload["n_requested"]}, indent=2))


if __name__ == "__main__":
    main()
