"""GO Annotation (GAF 2.2) parsing.

Real GAF 2.2 spec: https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
Columns (17, tab-separated):
1 DB, 2 DB_Object_ID, 3 DB_Object_Symbol, 4 Qualifier, 5 GO_ID,
6 DB:Reference, 7 Evidence Code, 8 With/From, 9 Aspect,
10 DB_Object_Name, 11 DB_Object_Synonym, 12 DB_Object_Type,
13 Taxon, 14 Date, 15 Assigned_By, 16 Annotation_Extension, 17 Gene_Product_Form_ID

Real downloads (no auth):
  https://current.geneontology.org/annotations/sgd.gaf.gz
  https://current.geneontology.org/annotations/mgi.gaf.gz
  https://current.geneontology.org/annotations/pombase.gaf.gz

Prediction target term example (real GO Biological Process):
  GO:0007049 — "cell cycle" (well-studied, hundreds of yeast genes)
Alternatives:
  GO:0006468 — protein phosphorylation
  GO:0006355 — regulation of transcription, DNA-templated
  GO:0007005 — mitochondrion organization
"""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# Known real GO terms (subset) with names
KNOWN_GO_TERMS: dict[str, str] = {
    "GO:0007049": "cell cycle",
    "GO:0006468": "protein phosphorylation",
    "GO:0006355": "regulation of transcription, DNA-templated",
    "GO:0007005": "mitochondrion organization",
    "GO:0008150": "biological_process",
    "GO:0003674": "molecular_function",
    "GO:0005575": "cellular_component",
}


@dataclass
class GOAnnotation:
    db: str
    gene_id: str  # DB:Object_ID combined or Object_ID
    db_object_id: str
    symbol: str
    go_id: str
    evidence_code: str
    aspect: str
    taxon: str  # e.g. taxon:559292
    qualifier: str = ""


def parse_gaf(
    path: str | Path,
    go_term: Optional[str] = None,
    evidence_codes: Optional[set[str]] = None,
    taxon_filter: Optional[str] = None,
) -> list[GOAnnotation]:
    """Parse a GAF 2.2 file (plain or gzipped).

    Args:
        path: path to .gaf or .gaf.gz
        go_term: if set, only keep rows matching this GO ID
        evidence_codes: if set, only keep rows with evidence code in set
        taxon_filter: if set, only keep rows where taxon contains this string
                      e.g. "559292" for S. cerevisiae
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    annotations: list[GOAnnotation] = []
    with opener(path, "rt") as f:  # type: ignore
        for line in f:
            if line.startswith("!"):
                continue
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 13:
                # Skip malformed lines
                continue
            # Pad to 17 if shorter (some files omit trailing empty columns)
            while len(cols) < 17:
                cols.append("")
            db = cols[0]
            db_object_id = cols[1]
            symbol = cols[2]
            qualifier = cols[3]
            go_id = cols[4]
            evidence = cols[6]
            aspect = cols[8]
            taxon = cols[12]

            # Skip NOT qualifier (negated annotations)
            if "NOT" in qualifier:
                continue

            if go_term and go_id != go_term:
                continue
            if evidence_codes and evidence not in evidence_codes:
                continue
            if taxon_filter and taxon_filter not in taxon:
                continue

            # Combine DB and object ID for unique gene key
            gene_id = f"{db}:{db_object_id}" if db else db_object_id

            annotations.append(
                GOAnnotation(
                    db=db,
                    gene_id=gene_id,
                    db_object_id=db_object_id,
                    symbol=symbol,
                    go_id=go_id,
                    evidence_code=evidence,
                    aspect=aspect,
                    taxon=taxon,
                    qualifier=qualifier,
                )
            )
    return annotations


def annotations_to_labels(
    annotations: list[GOAnnotation],
    gene_ids: list[str],
    go_term: str,
) -> dict[str, int]:
    """Create binary label dict for given gene list and GO term.

    Returns dict gene_id -> 0/1.
    Genes with at least one annotation to go_term => 1, else 0.
    Unknown genes (no annotation row at all) => 0.
    """
    positive: set[str] = set()
    for ann in annotations:
        if ann.go_id == go_term:
            positive.add(ann.gene_id)
    return {g: (1 if g in positive else 0) for g in gene_ids}


def load_labels_for_graph(
    gaf_path: str | Path,
    gene_ids: list[str],
    go_term: str = "GO:0007049",
    evidence_codes: Optional[set[str]] = None,
    taxon_filter: Optional[str] = None,
) -> dict[str, int]:
    """Convenience: parse GAF and produce labels for gene_ids."""
    anns = parse_gaf(gaf_path, go_term=go_term, evidence_codes=evidence_codes, taxon_filter=taxon_filter)
    # Need all annotations to determine positives; but if go_term is set,
    # positives are already filtered. For label dict we also need to know
    # which genes are positive vs negative — if filtered by go_term we treat
    # non-positive as 0. So pass full gene list.
    positives = {a.gene_id for a in anns}
    return {g: (1 if g in positives else 0) for g in gene_ids}


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse GO GAF file")
    parser.add_argument("--go-annotations-path", required=True)
    parser.add_argument("--go-term", default="GO:0007049")
    parser.add_argument("--evidence-codes", default=None, help="Comma-separated, e.g. IDA,IMP")
    parser.add_argument("--taxon-filter", default=None)
    args = parser.parse_args()

    ec = set(args.evidence_codes.split(",")) if args.evidence_codes else None
    anns = parse_gaf(args.go_annotations_path, go_term=args.go_term, evidence_codes=ec, taxon_filter=args.taxon_filter)
    print(f"Parsed {len(anns)} annotations for {args.go_term} from {args.go_annotations_path}")
    for a in anns[:5]:
        print(a)


if __name__ == "__main__":
    main()
