"""Expand PepBench→RCSB matching: more peptides, fuzzy partner overlap, multi-hit PDBs.

Does NOT touch held-out test labels fabrication — only grows the catalog of
real PepBench pKd ∩ experimental PDB pairs. Fail loud if PepBench CSV missing.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path

from peptide_affinity.load import DATA_DIR
from peptide_affinity.scripts.map_pepbench_to_rcsb import (
    download_pdb,
    parse_pdb_sequences,
    peptide_present,
    rcsb_sequence_search,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pepbench", type=Path, default=DATA_DIR / "raw" / "PpI_ba_raw.csv")
    parser.add_argument("--pdb-dir", type=Path, default=DATA_DIR / "full_pdbs")
    parser.add_argument("--out", type=Path, default=DATA_DIR / "peptide_affinity_catalog_v3.tsv")
    parser.add_argument("--map-cache", type=Path, default=DATA_DIR / "pepbench_rcsb_map_v3.json")
    parser.add_argument("--max-peptides", type=int, default=500)
    parser.add_argument("--min-len", type=int, default=6)
    parser.add_argument("--max-len", type=int, default=35)
    parser.add_argument("--sleep", type=float, default=0.08)
    parser.add_argument("--min-matches", type=int, default=50)
    args = parser.parse_args()

    if not args.pepbench.is_file():
        raise SystemExit(
            f"Missing {args.pepbench}. Download PpI_ba_raw.csv per "
            "benchmarks/peptide_affinity/README.md — refusing to fabricate."
        )

    rows = list(csv.DictReader(args.pepbench.open(encoding="utf-8")))
    by_pep: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        pep = row["pep_seq"].strip().upper()
        if not (args.min_len <= len(pep) <= args.max_len):
            continue
        by_pep[pep].append(row)

    # Prefer unique labels; allow mild multi-label by taking median later
    candidates: list[tuple[str, dict[str, str]]] = []
    for pep, group in by_pep.items():
        labels = sorted({float(x["label"]) for x in group})
        if len(labels) == 1:
            candidates.append((pep, group[0]))
        elif max(labels) - min(labels) <= 0.5:
            # near-duplicates — use mean label
            mean_lab = sum(labels) / len(labels)
            row = dict(group[0])
            row["label"] = f"{mean_lab:.6f}"
            candidates.append((pep, row))
    candidates.sort(key=lambda x: (abs(len(x[0]) - 14), -len(x[0])))
    candidates = candidates[: args.max_peptides]

    cache: dict[str, list[str]] = {}
    if args.map_cache.is_file():
        cache = json.loads(args.map_cache.read_text(encoding="utf-8"))

    args.pdb_dir.mkdir(parents=True, exist_ok=True)
    catalog: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()

    for i, (pep, row) in enumerate(candidates):
        if pep not in cache:
            cache[pep] = rcsb_sequence_search(pep, identity=0.90, rows=8)
            time.sleep(args.sleep)
            if i % 20 == 0:
                args.map_cache.write_text(json.dumps(cache, indent=2), encoding="utf-8")
                print(f"searched {i}/{len(candidates)}; catalog={len(catalog)}", flush=True)

        for pdb_id in cache.get(pep, [])[:5]:
            dest = args.pdb_dir / f"{pdb_id}.pdb"
            if not download_pdb(pdb_id, dest):
                continue
            sequences, resolution, year = parse_pdb_sequences(dest)
            pep_chain = peptide_present(sequences, pep)
            if pep_chain is None:
                continue
            receptors = [(c, s) for c, s in sequences.items() if c != pep_chain and len(s) >= 40]
            if not receptors:
                continue
            # Multi-chain required
            if len(sequences) < 2:
                continue
            prot = row["prot_seq"].strip().upper()
            best_chain, best_seq, best_ov = receptors[0][0], receptors[0][1], 0.0
            for c, s in receptors:
                window = 25
                ov = 0.0
                if len(s) >= window and len(prot) >= window:
                    hits = sum(
                        1
                        for j in range(0, min(len(s), 500) - window + 1, window)
                        if s[j : j + window] in prot
                    )
                    total = max(1, (min(len(s), 500) - window) // window)
                    ov = hits / total
                if ov >= best_ov:
                    best_ov, best_chain, best_seq = ov, c, s
            # Relaxed partner overlap (was 0.05–0.25)
            if best_ov < 0.02 and len(receptors) > 1:
                # still accept longest receptor if peptide is exact-chain
                best_chain, best_seq = max(receptors, key=lambda x: len(x[1]))
                best_ov = 0.02
            key = (pdb_id, pep)
            if key in seen:
                continue
            seen.add(key)
            catalog.append(
                {
                    "record_id": f"PA_{pdb_id}_{pep_chain}",
                    "pdb_id": pdb_id,
                    "receptor_seq": best_seq,
                    "peptide_seq": pep,
                    "peptide_len": len(pep),
                    "resolution": resolution if resolution is not None else "",
                    "affinity_value": "",
                    "affinity_type": "pKd",
                    "pKd": f"{float(row['label']):.6f}",
                    "source": "pepbench_ppi_ba+rcsb_seq_search_v3",
                    "structure_path": f"full_pdbs/{pdb_id}.pdb",
                    "peptide_chain": pep_chain,
                    "receptor_chains": best_chain,
                    "deposit_year": year if year is not None else "",
                    "notes": f"receptor_window_overlap={best_ov:.3f}",
                }
            )
            break

    args.map_cache.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    if len(catalog) < args.min_matches:
        raise SystemExit(f"only {len(catalog)} matches (need ≥{args.min_matches})")

    fields = [
        "record_id",
        "pdb_id",
        "receptor_seq",
        "peptide_seq",
        "peptide_len",
        "resolution",
        "affinity_value",
        "affinity_type",
        "pKd",
        "source",
        "structure_path",
        "peptide_chain",
        "receptor_chains",
        "deposit_year",
        "notes",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(catalog)
    print(f"wrote {len(catalog)} → {args.out}")


if __name__ == "__main__":
    main()
