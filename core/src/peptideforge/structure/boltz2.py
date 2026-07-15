"""Boltz-2 wrapper for peptide–target complex folding.

Invokes the open-source ``boltz predict`` CLI (MIT) on a YAML input built from
the peptide sequence and receptor sequence extracted from ``target_structure``.
Missing CLI or failed runs raise loudly — no fabricated coordinates.

LICENSE: Boltz-2 is MIT-licensed (jwohlwend/boltz). AlphaFold3 weights are not used.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

from peptideforge.contracts.models import ComplexStructure, PeptideCandidate, Provenance
from peptideforge.structure.cache import FoldCache, FoldCacheKey
from peptideforge.structure.deps import Boltz2UnavailableError, require_boltz_cli
from peptideforge.structure.pdb_utils import (
    extract_chain_sequence,
    list_protein_chains,
    mean_plddt_from_pdb,
    read_pdb_text,
)


class Boltz2StructurePredictor:
    """Fold peptide + receptor with Boltz-2; cache by sequence/target hash."""

    fold_method = "boltz2"

    def __init__(
        self,
        *,
        receptor_chain: str | None = None,
        peptide_chain_id: str = "P",
        receptor_chain_id: str = "R",
        cache: FoldCache | None = None,
        use_msa_server: bool = True,
        boltz_extra_args: tuple[str, ...] = (),
    ) -> None:
        self.receptor_chain = receptor_chain
        self.peptide_chain_id = peptide_chain_id
        self.receptor_chain_id = receptor_chain_id
        self.cache = cache
        self.use_msa_server = use_msa_server
        self.boltz_extra_args = boltz_extra_args

    def fold(
        self,
        candidate: PeptideCandidate,
        *,
        target_id: str,
        target_structure: str,
        seed: int | None = None,
    ) -> ComplexStructure:
        key = FoldCacheKey(
            sequence=candidate.sequence,
            target_id=target_id,
            target_structure=target_structure,
        )
        if self.cache is not None:
            hit = self.cache.get(key)
            if hit is not None:
                return hit

        boltz_bin = require_boltz_cli()
        target_path = Path(target_structure)
        if not target_path.is_file():
            raise FileNotFoundError(f"target_structure not found: {target_structure}")

        pdb_text = read_pdb_text(target_path)
        chains = list_protein_chains(pdb_text)
        receptor_chain = self.receptor_chain
        if receptor_chain is None:
            # Heuristic: longest protein chain in the template (typically receptor)
            receptor_chain = max(chains, key=lambda c: len(extract_chain_sequence(pdb_text, c)))
        receptor_seq = extract_chain_sequence(pdb_text, receptor_chain)

        with tempfile.TemporaryDirectory(prefix="peptideforge_boltz_") as tmp:
            tmp_path = Path(tmp)
            yaml_path = tmp_path / "input.yaml"
            out_dir = tmp_path / "out"
            _write_boltz_yaml(
                yaml_path,
                peptide_seq=candidate.sequence,
                peptide_id=self.peptide_chain_id,
                receptor_seq=receptor_seq,
                receptor_id=self.receptor_chain_id,
                template_pdb=target_path.resolve(),
                receptor_template_chain=receptor_chain,
            )
            cmd = [
                boltz_bin,
                "predict",
                str(yaml_path),
                "--out_dir",
                str(out_dir),
            ]
            if self.use_msa_server:
                cmd.append("--use_msa_server")
            if seed is not None:
                cmd.extend(["--seed", str(seed)])
            cmd.extend(self.boltz_extra_args)

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "Boltz-2 prediction failed "
                    f"(exit {proc.returncode}). stderr:\n{proc.stderr[-4000:]}"
                )

            pred_pdb = _find_predicted_pdb(out_dir)
            pred_text = read_pdb_text(pred_pdb)
            confidence = _confidence_from_output(out_dir, pred_text)

        structure = ComplexStructure(
            candidate_id=candidate.candidate_id,
            target_id=target_id,
            sequence=candidate.sequence,
            pdb_path=str(pred_pdb),
            pdb_text=pred_text,
            confidence=confidence,
            fold_method=self.fold_method,
            cache_key=key.digest(),
            provenance=Provenance(
                tool_versions={
                    "boltz_cli": boltz_bin,
                    "receptor_chain": receptor_chain,
                },
            ),
        )
        if self.cache is not None:
            self.cache.put(key, structure)
        return structure


def _write_boltz_yaml(
    path: Path,
    *,
    peptide_seq: str,
    peptide_id: str,
    receptor_seq: str,
    receptor_id: str,
    template_pdb: Path,
    receptor_template_chain: str,
) -> None:
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        f"      id: {peptide_id}",
        f"      sequence: {peptide_seq}",
        "      msa: empty",
        "  - protein:",
        f"      id: {receptor_id}",
        f"      sequence: {receptor_seq}",
        "      msa: empty",
        "templates:",
        "  - pdb: " + str(template_pdb),
        "    chain_id:",
        f"      - {receptor_id}",
        "    template_id:",
        f"      - {receptor_template_chain}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _find_predicted_pdb(out_dir: Path) -> Path:
    pdbs = sorted(out_dir.rglob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"Boltz-2 produced no PDB under {out_dir}")
    # Prefer model_0 if present
    for p in pdbs:
        if "model_0" in p.name or p.name.endswith("_0.pdb"):
            return p
    return pdbs[0]


def _confidence_from_output(out_dir: Path, pdb_text: str) -> float:
    json_files = sorted(out_dir.rglob("*.json"))
    for jpath in json_files:
        try:
            data = json.loads(jpath.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        for key in ("confidence", "confidence_score", "plddt", "complex_plddt"):
            if key in data and isinstance(data[key], (int, float)):
                val = float(data[key])
                return val if val <= 1.0 else val / 100.0
    return mean_plddt_from_pdb(pdb_text)
