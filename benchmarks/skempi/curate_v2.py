"""Phase 5A.2 — true-count curation for the SKEMPI protein–peptide subset.

Goal: enumerate how many *unique PDBs* in SKEMPI v2 qualify for the peptide
within-target ΔΔG test **before** spending any Boltz folding compute.

Criteria (cheap, structure-file based):
- WT crystal structure file exists under --structure-dir as {pdb_id}.pdb
- The mutated partner chain (mutant code like LI18G → chain 'I') is a peptide
  with length in [min_len, max_len] (default 5–50), measured from ATOM records.
- There are ≥K usable mutations on that PDB after the above filters (default K=8).

This does *not* attempt a full geometric interface test (contacts) — it is a
true-count upper bound based on chain length and file availability, suitable for
deciding whether a powered PDB-level evaluation is feasible from SKEMPI alone.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from peptideforge.structure.pdb_utils import extract_chain_sequence, read_pdb_text
from peptideforge_benchmarks.skempi import load_skempi_ddg


@dataclass(frozen=True)
class PdbSummary:
    pdb_id: str
    n_mutations: int
    peptide_chain: str
    peptide_len: int
    example_record_id: str


def _peptide_len_for_mutation_chain(pdb_path: Path, chain_id: str) -> int | None:
    try:
        text = read_pdb_text(pdb_path)
        seq = extract_chain_sequence(text, chain_id)
        return len(seq)
    except Exception:  # noqa: BLE001
        return None


def curate(
    *,
    skempi_tsv: Path,
    structure_dir: Path,
    min_len: int,
    max_len: int,
    min_mutations_per_pdb: int,
) -> dict[str, Any]:
    recs = load_skempi_ddg(skempi_tsv)
    by_pdb: dict[str, list[Any]] = defaultdict(list)
    for r in recs:
        by_pdb[r.pdb_id].append(r)

    qualifying_records: list[str] = []
    summaries: list[PdbSummary] = []
    missing_struct = 0
    chain_len_fail = 0

    for pdb_id, rows in sorted(by_pdb.items()):
        pdb_path = structure_dir / f"{pdb_id}.pdb"
        if not pdb_path.is_file():
            missing_struct += 1
            continue

        kept: list[Any] = []
        pep_chain = None
        pep_len = None
        for r in rows:
            if len(r.mutant) < 2:
                continue
            chain = r.mutant[1]
            L = _peptide_len_for_mutation_chain(pdb_path, chain)
            if L is None or not (min_len <= L <= max_len):
                continue
            kept.append(r)
            pep_chain = chain
            pep_len = L

        if len(kept) < min_mutations_per_pdb:
            if len(rows) > 0:
                chain_len_fail += 1
            continue

        qualifying_records.extend([r.record_id for r in kept])
        summaries.append(
            PdbSummary(
                pdb_id=pdb_id,
                n_mutations=len(kept),
                peptide_chain=str(pep_chain),
                peptide_len=int(pep_len),
                example_record_id=kept[0].record_id,
            )
        )

    summaries.sort(key=lambda s: (-s.n_mutations, s.pdb_id))
    payload = {
        "status": "OK",
        "skempi_tsv": str(skempi_tsv),
        "structure_dir": str(structure_dir),
        "criteria": {
            "peptide_len_range": [min_len, max_len],
            "min_mutations_per_pdb": min_mutations_per_pdb,
            "peptide_chain": "mutation chain (2nd char of SKEMPI mutant code)",
        },
        "counts": {
            "n_records_total": len(recs),
            "n_unique_pdb_total": len(by_pdb),
            "n_unique_pdb_missing_structure": missing_struct,
            "n_unique_pdb_failing_len_or_K": chain_len_fail,
            "n_unique_pdb_qualifying": len(summaries),
            "n_records_qualifying": len(qualifying_records),
        },
        "qualifying_pdbs": [asdict(s) for s in summaries],
        "qualifying_record_ids": qualifying_records,
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
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
    parser.add_argument("--min-len", type=int, default=5)
    parser.add_argument("--max-len", type=int, default=50)
    parser.add_argument("-K", "--min-mutations-per-pdb", type=int, default=8)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("benchmarks/skempi/data/curate_v2_last_run.json"),
    )
    args = parser.parse_args()

    payload = curate(
        skempi_tsv=args.skempi_tsv,
        structure_dir=args.structure_dir,
        min_len=args.min_len,
        max_len=args.max_len,
        min_mutations_per_pdb=args.min_mutations_per_pdb,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=== Phase 5A.2 true-count curation (cheap) ===")
    c = payload["counts"]
    print(
        f"qualifying unique PDBs = {c['n_unique_pdb_qualifying']} "
        f"(records={c['n_records_qualifying']}) "
        f"out of total unique PDBs = {c['n_unique_pdb_total']}"
    )
    top = payload["qualifying_pdbs"][:10]
    if top:
        print("Top qualifying PDBs (by #mutations):")
        for row in top:
            print(
                f"  {row['pdb_id']}: n={row['n_mutations']}, peptide_chain={row['peptide_chain']}, "
                f"len={row['peptide_len']}"
            )
    else:
        print("No qualifying PDBs found under current criteria.")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

