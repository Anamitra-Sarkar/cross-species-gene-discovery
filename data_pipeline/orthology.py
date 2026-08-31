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
    """Parse simplified 6-column edge list TSV.

    Hardened for real OrthoDB-like edge files:
    - comment lines with optional leading whitespace before '#'
    - BOM stripping, empty/whitespace-only lines
    - both tab and whitespace splitting with consistent handling
    - validates score range [0, 1] and non-empty required fields
    - includes line numbers in errors
    """
    edges: list[OrthologEdge] = []
    with open(path, encoding="utf-8-sig") as f:
        for lineno, raw in enumerate(f, start=1):
            # Strip BOM on first line via utf-8-sig; also handle trailing newline
            stripped = raw.strip()
            if not stripped:
                continue
            # Comment with optional leading whitespace
            if stripped.lstrip().startswith("#"):
                continue
            # Also skip leading-whitespace comment: raw.lstrip startswith '#'
            if raw.lstrip().startswith("#"):
                continue
            # Split: prefer tab if present, else any whitespace
            if "\t" in raw:
                parts = [p.strip() for p in raw.strip().split("\t")]
                # Remove empty parts caused by consecutive tabs / trailing tab
                # but preserve field count: filter only truly-empty from consecutive tabs yields wrong count,
                # so we filter empties only if it would collapse; instead keep as-is and validate below.
                # Remove completely empty strings from split artifacts while preserving real empties between tabs?
                # For robustness, if line had many tabs, split already handles it; just drop empty trailing.
                parts = [p for p in parts if p != "" or len(parts) < 6]
                # Re-split more strictly if still short: try whitespace
                if len(parts) < 6:
                    parts = raw.strip().split()
            else:
                parts = raw.strip().split()
            if len(parts) < 6:
                raise ValueError(f"{path}:{lineno}: Expected 6 columns, got {len(parts)}: {raw.strip()!r}")
            og_id, gene_a, species_a, gene_b, species_b, score_s = [p.strip() for p in parts[:6]]
            if not og_id or not gene_a or not gene_b:
                raise ValueError(f"{path}:{lineno}: og_id and gene IDs must be non-empty: {raw.strip()!r}")
            if not species_a or not species_b:
                raise ValueError(f"{path}:{lineno}: species columns must be non-empty: {raw.strip()!r}")
            try:
                score = float(score_s)
            except ValueError as e:
                raise ValueError(f"{path}:{lineno}: invalid score {score_s!r}: {e}") from e
            if not (0.0 <= score <= 1.0):
                raise ValueError(f"{path}:{lineno}: score {score} out of range [0, 1]: {raw.strip()!r}")
            edges.append(
                OrthologEdge(
                    og_id=og_id,
                    gene_a=gene_a,
                    species_a=species_a,
                    gene_b=gene_b,
                    species_b=species_b,
                    score=score,
                )
            )
    return edges


def parse_orthodb_tab(path: str | Path) -> list[OrthologEdge]:
    """Parse native OrthoDB *_OGs.tab (one row per gene-OG membership).

    Expected columns (tab-separated, header present):
    e.g. Level, OG_ID, Organism_ID, Gene_ID, UniProt, ... (varies)
    We detect columns by header names case-insensitively.

    Hardened for real file quirks:
    - comment lines starting with '#' or '!' before header are skipped
    - BOM stripped, leading/trailing whitespace on headers
    - header detection is case/whitespace-insensitive with substring fallback

    Expands each OG's gene list into pairwise edges (complete graph per OG,
    score=1.0 since OrthoDB bulk does not include per-pair confidence).
    """
    og_to_genes: dict[str, list[tuple[str, str]]] = {}  # og_id -> [(gene_id, taxid)]

    with open(path, encoding="utf-8-sig") as f:
        # Skip leading comment lines before header (real OrthoDB files sometimes have them)
        # We need to find the header line ourselves because csv.DictReader treats first non-empty as header.
        # Peek and skip comment lines.
        first_lines: list[str] = []
        header_line: str | None = None
        header_offset = 0
        for line in f:
            if not line.strip():
                header_offset += 1
                continue
            if line.lstrip().startswith("#") or line.lstrip().startswith("!"):
                header_offset += 1
                continue
            header_line = line
            break
        if header_line is None:
            raise ValueError("Empty or headerless OrthoDB TAB file")
        # Rewind and create DictReader starting at header
        f.seek(0)
        # Skip the comment/blank lines we counted
        for _ in range(header_offset):
            next(f)
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Empty or headerless OrthoDB TAB file")
        # Strip whitespace/BOM from fieldnames
        reader.fieldnames = [h.strip().lstrip("\ufeff") for h in reader.fieldnames if h is not None]
        # Detect OG, gene, and organism/taxon columns
        og_col = _find_col(reader.fieldnames, ["og id", "og_id", "orthodb group", "odbgroup", "og"])
        gene_col = _find_col(reader.fieldnames, ["gene id", "gene_id", "gene", "protein id"])
        taxon_col = _find_col(reader.fieldnames, ["organism", "taxid", "ncbi taxid", "organism id", "species", "tax_id"])

        if og_col is None or gene_col is None:
            raise ValueError(
                f"Cannot detect OG or gene column. Headers: {reader.fieldnames}"
            )

        for row in reader:
            # Skip comment lines that may appear mid-file (real OrthoDB exports occasionally have them)
            # csv.DictReader will parse them as rows; detect via first column starting with # or !
            first_val = next((v for v in row.values() if v is not None), "")
            if isinstance(first_val, str) and first_val.lstrip().startswith(("#", "!")):
                continue
            og_id = (row.get(og_col) or "").strip()
            gene_id = (row.get(gene_col) or "").strip()
            taxid = (row.get(taxon_col) or "").strip() if taxon_col and taxon_col in row else "unknown"
            # Some OrthoDB rows have taxid like "taxon:559292" or plain number; normalize to plain number string
            if taxid.startswith("taxon:"):
                taxid = taxid.split(":", 1)[1]
            taxid = taxid.strip() or "unknown"
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
