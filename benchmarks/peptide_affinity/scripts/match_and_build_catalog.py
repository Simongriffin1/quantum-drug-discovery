"""Match RCSB peptide–protein PDBs to PepBench PpI_ba pKd labels and write catalog TSV."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from peptide_affinity.load import DATA_DIR, affinity_to_pkd
from peptide_affinity.models import AffinityType

AA3TO1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_pdb_sequences(pdb_path: Path) -> tuple[dict[str, str], float | None, int | None]:
    """Return chain_id → sequence (standard AA only), resolution, deposit year."""
    chains: dict[str, list[tuple[int, str]]] = defaultdict(list)
    seen: dict[str, set[int]] = defaultdict(set)
    resolution: float | None = None
    year: int | None = None
    seqres: dict[str, list[str]] = defaultdict(list)

    text = pdb_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith("REMARK   2 RESOLUTION") and resolution is None:
            m = re.search(r"RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS", line)
            if m:
                resolution = float(m.group(1))
        if line.startswith("HEADER") and year is None:
            # HEADER .... DD-MMM-YY
            m = re.search(r"(\d{2})-([A-Z]{3})-(\d{2})\s*$", line.strip())
            if m:
                yy = int(m.group(3))
                year = 1900 + yy if yy > 70 else 2000 + yy
        if line.startswith("SEQRES"):
            cid = line[11]
            for tok in line[19:].split():
                if tok in AA3TO1:
                    seqres[cid].append(AA3TO1[tok])
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname not in AA3TO1:
                continue
            cid = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            if resseq in seen[cid]:
                continue
            seen[cid].add(resseq)
            chains[cid].append((resseq, AA3TO1[resname]))

    sequences: dict[str, str] = {}
    for cid, residues in chains.items():
        residues_sorted = sorted(residues, key=lambda x: x[0])
        sequences[cid] = "".join(aa for _, aa in residues_sorted)
    # Prefer SEQRES when available (more complete)
    for cid, aas in seqres.items():
        if len(aas) >= len(sequences.get(cid, "")):
            sequences[cid] = "".join(aas)
    return sequences, resolution, year


def load_pepbench(path: Path) -> dict[str, list[tuple[str, float]]]:
    """Map peptide_seq → list of (receptor_seq, pKd)."""
    if not path.is_file():
        raise FileNotFoundError(
            f"PepBench PpI_ba raw CSV missing: {path}\n"
            "Download:\n"
            "  curl -L -o benchmarks/peptide_affinity/data/raw/PpI_ba_raw.csv \\\n"
            "    https://huggingface.co/datasets/jiahuizhang/PepBenchData_raw/"
            "resolve/main/data/PepPI/nature/PpI_ba/raw.csv\n"
            "See benchmarks/peptide_affinity/README.md"
        )
    by_pep: dict[str, list[tuple[str, float]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pep = row["pep_seq"].strip().upper()
            prot = row["prot_seq"].strip().upper()
            if not (5 <= len(pep) <= 50):
                continue
            by_pep[pep].append((prot, float(row["label"])))
    return by_pep


def receptor_overlap(a: str, b: str, window: int = 40) -> float:
    """Fraction of overlapping windows present — cheap homology score."""
    if not a or not b:
        return 0.0
    short, long = (a, b) if len(a) <= len(b) else (b, a)
    if len(short) < window:
        return 1.0 if short in long else 0.0
    hits = 0
    total = 0
    step = max(1, window // 2)
    for i in range(0, len(short) - window + 1, step):
        total += 1
        if short[i : i + window] in long:
            hits += 1
    return hits / total if total else 0.0


def match_pdb(
    pdb_id: str,
    pdb_path: Path,
    by_pep: dict[str, list[tuple[str, float]]],
    *,
    min_overlap: float = 0.25,
) -> dict[str, object] | None:
    sequences, resolution, year = parse_pdb_sequences(pdb_path)
    peptides = [(cid, seq) for cid, seq in sequences.items() if 5 <= len(seq) <= 50]
    receptors = [(cid, seq) for cid, seq in sequences.items() if len(seq) >= 50]
    if not peptides or not receptors:
        return None

    # Prefer the shortest peptide chain as ligand
    peptides.sort(key=lambda x: len(x[1]))
    pep_chain, pep_seq = peptides[0]
    if pep_seq not in by_pep:
        return None

    # Best receptor by sequence overlap with PepBench partner
    best: tuple[float, str, str, float] | None = None  # overlap, rec_chain, rec_seq, pkd
    for prot, pkd in by_pep[pep_seq]:
        for rec_chain, rec_seq in receptors:
            ov = receptor_overlap(rec_seq, prot)
            if best is None or ov > best[0]:
                best = (ov, rec_chain, rec_seq, pkd)
    if best is None or best[0] < min_overlap:
        return None

    ov, rec_chain, rec_seq, pkd = best
    return {
        "record_id": f"PA_{pdb_id}_{pep_chain}",
        "pdb_id": pdb_id,
        "receptor_seq": rec_seq,
        "peptide_seq": pep_seq,
        "peptide_len": len(pep_seq),
        "resolution": resolution if resolution is not None else "",
        "affinity_value": "",
        "affinity_type": "pKd",
        "pKd": f"{pkd:.6f}",
        "source": "pepbench_ppi_ba+rcsb_match",
        "structure_path": f"full_pdbs/{pdb_id}.pdb",
        "peptide_chain": pep_chain,
        "receptor_chains": rec_chain,
        "deposit_year": year if year is not None else "",
        "notes": f"receptor_overlap={ov:.3f}",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pepbench",
        type=Path,
        default=DATA_DIR / "raw" / "PpI_ba_raw.csv",
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=DATA_DIR / "full_pdbs",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DATA_DIR / "peptide_affinity_catalog_v2.tsv",
    )
    parser.add_argument("--min-overlap", type=float, default=0.25)
    args = parser.parse_args()

    by_pep = load_pepbench(args.pepbench)
    rows: list[dict[str, object]] = []
    seen_keys: set[tuple[str, str]] = set()
    for pdb_path in sorted(args.pdb_dir.glob("*.pdb")):
        pdb_id = pdb_path.stem.upper()
        matched = match_pdb(pdb_id, pdb_path, by_pep, min_overlap=args.min_overlap)
        if matched is None:
            continue
        key = (str(matched["pdb_id"]), str(matched["peptide_seq"]))
        if key in seen_keys:
            continue
        seen_keys.add(key)
        rows.append(matched)

    if len(rows) < 2:
        raise SystemExit(
            f"matched only {len(rows)} complexes under {args.pdb_dir} — "
            "download more candidate PDBs or lower --min-overlap (after QC)."
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
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
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} matched entries → {args.out}")
    _ = affinity_to_pkd  # silence unused in some tooling
    _ = AffinityType


if __name__ == "__main__":
    main()
