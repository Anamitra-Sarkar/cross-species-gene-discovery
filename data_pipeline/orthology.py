"""Orthology parsing: OrthoDB / eggNOG TSV -> edge list.

Supports two input formats:
1. Simplified edge list TSV (our fixture format):
   og_id  gene_id_1  species_1  gene_id_2  species_2  score
   where score is orthology confidence (0-1).

2. Native OrthoDB bulk TAB (odb*_OGs.tab):
   Columns vary by version but minimally:
   og_id, gene_id, taxid, ...  (one row per gene-OG membership)
   We expand each OG into pairwise edges (clique per OG).

Also documents the real OrthoDB REST API (not called in sandbox).
Real API: https://www.orthodb.org/?page=api
  GET /ogdetails?id=<OG_ID>
  GET /tab?id=<OG_ID>
  GET /search?query=<gene>&level=<clade>&species=<taxid>
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class OrthologEdge:
    og_id: str
    gene_a: str
    species_a: str
    gene_b: str
    species_b: str
    score: float  # orthology confidence 0-1


def parse_simplified_edge_list(path: str | Path) -> list[OrthologEdge]:
    """Parse simplified 6-column edge list TSV."""
    edges: list[OrthologEdge] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            # Support both TSV and whitespace-splitting fallback, but prefer tab
            # Re-split by tab if needed
            if len(parts) == 1:
                parts = line.split("\t")
            if len(parts) < 6:
                # Try splitting by tab strictly
                parts = line.split("\t")
            if len(parts) < 6:
                raise ValueError(f"Expected 6 columns, got {len(parts)}: {line!r}")
            og_id, gene_a, species_a, gene_b, species_b, score_s = parts[:6]
            edges.append(
                OrthologEdge(
                    og_id=og_id.strip(),
                    gene_a=gene_a.strip(),
                    species_a=species_a.strip(),
                    gene_b=gene_b.strip(),
                    species_b=species_b.strip(),
                    score=float(score_s),
                )
            )
    return edges


def parse_orthodb_tab(path: str | Path) -> list[OrthologEdge]:
    """Parse native OrthoDB *_OGs.tab (one row per gene-OG membership).

    Expected columns (tab-separated, header present):
    e.g. Level, OG_ID, Organism_ID, Gene_ID, UniProt, ... (varies)
    We detect columns by header names case-insensitively.

    Expands each OG's gene list into pairwise edges (complete graph per OG,
    score=1.0 since OrthoDB bulk does not include per-pair confidence).
    """
    og_to_genes: dict[str, list[tuple[str, str]]] = {}  # og_id -> [(gene_id, taxid)]

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty or headerless OrthoDB TAB file")
        # Normalize fieldnames
        lower_fields = [h.lower() for h in reader.fieldnames]
        # Detect OG, gene, and organism/taxon columns
        og_col = _find_col(reader.fieldnames, ["og id", "og_id", "orthodb group", "odbgroup", "og"])
        gene_col = _find_col(reader.fieldnames, ["gene id", "gene_id", "gene", "protein id"])
        taxon_col = _find_col(reader.fieldnames, ["organism", "taxid", "ncbi taxid", "organism id", "species", "tax_id"])

        if og_col is None or gene_col is None:
            raise ValueError(
                f"Cannot detect OG or gene column. Headers: {reader.fieldnames}"
            )

        for row in reader:
            og_id = row[og_col].strip()
            gene_id = row[gene_col].strip()
            taxid = row[taxon_col].strip() if taxon_col else "unknown"
            if not og_id or not gene_id:
                continue
            og_to_genes.setdefault(og_id, []).append((gene_id, taxid))

    edges: list[OrthologEdge] = []
    for og_id, genes in og_to_genes.items():
        # Clique expansion
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                gene_a, tax_a = genes[i]
                gene_b, tax_b = genes[j]
                # Only cross-species edges (orthology is cross-species)
                if tax_a == tax_b:
                    continue
                edges.append(
                    OrthologEdge(
                        og_id=og_id,
                        gene_a=gene_a,
                        species_a=tax_a,
                        gene_b=gene_b,
                        species_b=tax_b,
                        score=1.0,
                    )
                )
    return edges


def _find_col(fieldnames: list[str], candidates: list[str]) -> Optional[str]:
    lower_map = {f.lower().strip(): f for f in fieldnames}
    for cand in candidates:
        if cand in lower_map:
            return lower_map[cand]
    # Substring match
    for cand in candidates:
        for low, orig in lower_map.items():
            if cand in low or low in cand:
                return orig
    return None


def parse_orthology_file(path: str | Path) -> list[OrthologEdge]:
    """Auto-detect format and parse orthology file."""
    path = Path(path)
    # Peek at first non-comment line
    with open(path) as f:
        header = None
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            header = s
            break
    if header is None:
        return []
    # If header contains tab and looks like OG/gene headers, treat as OrthoDB TAB
    low = header.lower()
    if "og" in low and ("gene" in low or "protein" in low) and "\t" in header:
        # Could be either; check column count vs simplified format.
        # Simplified edge list header would be: og_id gene_id_1 ...
        # OrthoDB TAB header contains "Level" or "Organism" etc.
        if "organism" in low or "level" in low or "uniprot" in low:
            return parse_orthodb_tab(path)
    # Try simplified edge list: 6 tab-separated values with numeric last col
    parts = header.split("\t")
    if len(parts) >= 6:
        try:
            float(parts[5])
            return parse_simplified_edge_list(path)
        except ValueError:
            pass
    # Fallback: try OrthoDB TAB
    try:
        return parse_orthodb_tab(path)
    except Exception:
        return parse_simplified_edge_list(path)


def edges_to_gene_set(edges: list[OrthologEdge]) -> dict[str, str]:
    """Return mapping gene_id -> species (taxid string)."""
    gene_species: dict[str, str] = {}
    for e in edges:
        gene_species[e.gene_a] = e.species_a
        gene_species[e.gene_b] = e.species_b
    return gene_species


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse orthology file to edge list JSON/TSV")
    parser.add_argument("--orthology-path", required=True, help="Path to orthology TSV/TAB file")
    parser.add_argument("--output", required=False, default=None, help="Output path (JSON)")
    args = parser.parse_args()

    edges = parse_orthology_file(args.orthology_path)
    print(f"Parsed {len(edges)} ortholog edges from {args.orthology_path}")
    genes = edges_to_gene_set(edges)
    species_counts: dict[str, int] = {}
    for g, s in genes.items():
        species_counts[s] = species_counts.get(s, 0) + 1
    print(f"Genes: {len(genes)}, per species: {species_counts}")

    if args.output:
        import json

        out = [e.__dict__ for e in edges]
        Path(args.output).write_text(json.dumps(out, indent=2))
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
