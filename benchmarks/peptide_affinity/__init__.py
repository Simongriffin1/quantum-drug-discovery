"""Expanded protein–peptide affinity benchmark (experimental structures + pKd)."""

from __future__ import annotations

__all__ = ["PeptideAffinityEntry", "load_peptide_affinity_catalog"]


def __getattr__(name: str):  # type: ignore[no-untyped-def]
    if name == "PeptideAffinityEntry":
        from peptide_affinity.models import PeptideAffinityEntry

        return PeptideAffinityEntry
    if name == "load_peptide_affinity_catalog":
        from peptide_affinity.load import load_peptide_affinity_catalog

        return load_peptide_affinity_catalog
    raise AttributeError(name)
