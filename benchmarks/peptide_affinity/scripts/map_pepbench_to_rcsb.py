"""Map PepBench peptides → RCSB structures, download PDBs, write catalog TSV."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

from peptide_affinity.load import DATA_DIR
from peptide_affinity.scripts.match_and_build_catalog import parse_pdb_sequences


def rcsb_sequence_search(seq: str, *, identity: float = 0.95, rows: int = 5) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1.0,
                "identity_cutoff": identity,
                "sequence_type": "protein",
                "value": seq,
            },
        },
        # polymer_entity returns IDs like 1ABC_1; entry return_type often 204s on short peptides
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    req = urllib.request.Request(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(query).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=45) as resp:
            raw = resp.read()
            if not raw:
                return []
            payload = json.loads(raw.decode())
    except urllib.error.HTTPError as exc:
        if exc.code in {204, 400, 404}:
            return []
        raise
    except Exception:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for item in payload.get("result_set", []):
        ident = str(item.get("identifier", "")).upper()
        pdb_id = ident.split("_")[0]
        if len(pdb_id) == 4 and pdb_id not in seen:
            seen.add(pdb_id)
            out.append(pdb_id)
    return out


def download_pdb(pdb_id: str, dest: Path) -> bool:
    if dest.is_file() and dest.stat().st_size > 1000:
        return True
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            data = resp.read()
        if len(data) < 500:
            return False
        dest.write_bytes(data)
        return True
    except Exception:
        return False


def peptide_present(sequences: dict[str, str], pep: str) -> str | None:
    """Return chain ID whose sequence equals pep (separate peptide chain required)."""
    exact = [
        (cid, seq)
        for cid, seq in sequences.items()
        if seq == pep and 5 <= len(seq) <= 50
    ]
    if not exact:
        return None
    exact.sort(key=lambda x: len(x[1]))
    return exact[0][0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pepbench", type=Path, default=DATA_DIR / "raw" / "PpI_ba_raw.csv")
    parser.add_argument("--pdb-dir", type=Path, default=DATA_DIR / "full_pdbs")
    parser.add_argument("--out", type=Path, default=DATA_DIR / "peptide_affinity_catalog_v2.tsv")
    parser.add_argument("--map-cache", type=Path, default=DATA_DIR / "pepbench_rcsb_map.json")
    parser.add_argument("--max-peptides", type=int, default=400)
    parser.add_argument("--sleep", type=float, default=0.15)
    parser.add_argument("--min-matches", type=int, default=30)
    args = parser.parse_args()

    rows = list(csv.DictReader(args.pepbench.open(encoding="utf-8")))
    by_pep: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        pep = row["pep_seq"].strip().upper()
        if not (8 <= len(pep) <= 30):
            continue
        by_pep.setdefault(pep, []).append(row)
    # Prefer mid-length peptides (12–25) — short epitopes often lack RCSB hits
    unique = [(p, v[0]) for p, v in by_pep.items() if len({float(x["label"]) for x in v}) == 1]
    unique.sort(key=lambda x: (abs(len(x[0]) - 16), -len(x[0])))
    unique = unique[: args.max_peptides]

    cache: dict[str, list[str]] = {}
    if args.map_cache.is_file():
        cache = json.loads(args.map_cache.read_text(encoding="utf-8"))

    args.pdb_dir.mkdir(parents=True, exist_ok=True)
    catalog_rows: list[dict[str, object]] = []
    seen_pdb_pep: set[tuple[str, str]] = set()

    for i, (pep, row) in enumerate(unique):
        if pep not in cache:
            try:
                cache[pep] = rcsb_sequence_search(pep)
            except Exception as exc:  # noqa: BLE001
                cache[pep] = []
                print(f"search fail {pep}: {exc}")
            time.sleep(args.sleep)
            if i % 25 == 0:
                args.map_cache.write_text(json.dumps(cache, indent=2), encoding="utf-8")
                print(f"searched {i}/{len(unique)}; catalog={len(catalog_rows)}")

        for pdb_id in cache.get(pep, [])[:3]:
            dest = args.pdb_dir / f"{pdb_id}.pdb"
            if not download_pdb(pdb_id, dest):
                continue
            sequences, resolution, year = parse_pdb_sequences(dest)
            pep_chain = peptide_present(sequences, pep)
            if pep_chain is None:
                continue
            # Receptor = longest other chain (or same-chain partner truncated)
            receptors = [(c, s) for c, s in sequences.items() if c != pep_chain and len(s) >= 40]
            if not receptors:
                # Same polymer: treat partner as receptor sequence from PepBench
                receptors = [("__partner__", row["prot_seq"].strip().upper())]
                # Keep full asymmetric unit as complex; chain annotation = pep_chain + all
                other = [c for c in sequences if c != pep_chain]
                rec_chain_label = ",".join(other) if other else pep_chain
            else:
                rec_chain_label = None
            prot = row["prot_seq"].strip().upper()
            best_chain, best_seq, best_ov = receptors[0][0], receptors[0][1], 0.0
            for c, s in receptors:
                window = 30
                ov = 0.0
                if len(s) >= window and len(prot) >= window:
                    hits = sum(
                        1
                        for j in range(0, min(len(s), 400) - window + 1, window)
                        if s[j : j + window] in prot
                    )
                    total = max(1, (min(len(s), 400) - window) // window)
                    ov = hits / total
                elif s == prot:
                    ov = 1.0
                if ov >= best_ov:
                    best_ov, best_chain, best_seq = ov, c, s
            # Require some receptor evidence OR embedded peptide in a multi-chain PDB
            if best_ov < 0.05 and len(sequences) < 2:
                continue
            key = (pdb_id, pep)
            if key in seen_pdb_pep:
                continue
            seen_pdb_pep.add(key)
            catalog_rows.append(
                {
                    "record_id": f"PA_{pdb_id}_{pep_chain}",
                    "pdb_id": pdb_id,
                    "receptor_seq": best_seq if best_chain != "__partner__" else prot,
                    "peptide_seq": pep,
                    "peptide_len": len(pep),
                    "resolution": resolution if resolution is not None else "",
                    "affinity_value": "",
                    "affinity_type": "pKd",
                    "pKd": f"{float(row['label']):.6f}",
                    "source": "pepbench_ppi_ba+rcsb_seq_search",
                    "structure_path": f"full_pdbs/{pdb_id}.pdb",
                    "peptide_chain": pep_chain,
                    "receptor_chains": rec_chain_label
                    or (best_chain if best_chain != "__partner__" else next(iter(sequences))),
                    "deposit_year": year if year is not None else "",
                    "notes": f"receptor_window_overlap={best_ov:.3f}",
                }
            )
            break  # one structure per peptide

    args.map_cache.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    if len(catalog_rows) < args.min_matches:
        raise SystemExit(
            f"Only matched {len(catalog_rows)} complexes (need ≥{args.min_matches}). "
            "Re-run with higher --max-peptides or check network/RCSB."
        )

    fieldnames = [
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(catalog_rows)
    print(f"wrote {len(catalog_rows)} entries → {args.out}")


if __name__ == "__main__":
    main()
