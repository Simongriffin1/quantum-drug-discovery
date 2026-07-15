# Peptide–protein affinity catalog (expanded)

This package assembles complexes that have **both** an experimental 3D structure
and an experimentally measured affinity (reported as **pKd = −log10(Kd[M])**),
for peptides of **5–50 residues**.

## Primary sources

| Source | Role | Access |
|---|---|---|
| **PepBenchmark PpI_ba** (Zhang et al. / PepBenchData_raw) | ~1,433 sequence–pKd pairs (PDBbind-derived lineage) | HuggingFace (free) |
| **RCSB PDB** | Experimental coordinates | Public download |
| **Propedia / PepBDB** | Optional structure QC cross-check | Manual |

PepBench PpI_ba provides `prot_seq,pep_seq,label` **without PDB IDs**. We obtain
structures by downloading candidate peptide–protein PDBs from RCSB, extracting
the peptide chain (length 5–50), and **matching `pep_seq` (+ receptor overlap)
to PepBench labels**. Entries are never fabricated.

## Manual download path (fail loud if missing)

```bash
# 1) PepBench affinity labels
mkdir -p benchmarks/peptide_affinity/data/raw
curl -L -o benchmarks/peptide_affinity/data/raw/PpI_ba_raw.csv \
  "https://huggingface.co/datasets/jiahuizhang/PepBenchData_raw/resolve/main/data/PepPI/nature/PpI_ba/raw.csv"

# 2) Candidate PDBs + match + build catalog
python -m peptide_affinity.scripts.match_and_build_catalog

# 3) Prepare structures (PDBFixer / optional PDB2PQR)
python -m peptide_affinity.prep --catalog benchmarks/peptide_affinity/data/peptide_affinity_catalog_v2.tsv
```

If HuggingFace / RCSB is unreachable, loaders raise `FileNotFoundError` with this
README path — they do **not** invent rows.

## CI fixture

`fixtures/peptide_affinity_ci_v1.tsv` is a small schema-valid subset for offline
tests (no network). It does **not** authorize the M2 oracle gate.

## License note

PepBenchData / PDBbind-derived labels: cite Zhang et al. PepBenchmark and PDBbind.
Atomic coordinates: RCSB PDB terms. Do not redistribute full PDBbind dumps without
checking their license.
