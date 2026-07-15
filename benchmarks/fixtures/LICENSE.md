# Benchmark fixtures for PeptideForge
#
# CI uses ONLY these local files — no live network downloads.
#
# Contents:
#   pdbbind_peptide_affinity_v1.tsv  — small protein–peptide affinity subset schema
#   skempi_ddg_v1.tsv                 — small mutation ΔΔG subset schema
#   clusters_v1.tsv                   — precomputed homology clusters for split tests
#
# IMPORTANT: Values in the affinity/ddG fixtures are **format-faithful public
# examples** for loader + harness plumbing. They are NOT a claim that PeptideForge
# has already passed the oracle-validity stage gate. Real oracle Spearman will be
# computed and logged in P3 against an explicitly named subset.
#
# Before redistributing any larger dump of PDBbind-derived labels, verify the
# PDBbind license. Prefer SKEMPI / PDB terms where possible.
#
# contact: see CURSOR_PROJECT_CONTEXT.md
