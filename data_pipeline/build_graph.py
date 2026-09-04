"""End-to-end CLI: orthology + GO annotations -> cross-species graph.

No downloads in sandbox; documents real-run procedure.
Uses synthetic fixtures for tests (matching real formats).

Real-run example:
  python -m data_pipeline.build_graph \\
    --orthology-path data/odb_Eukaryota_OGs.tab \\
    --go-annotations-path data/sgd.gaf \\
    --go-term GO:0007049 \\
    --evidence-codes IDA,IMP,IGI,IPI \\
    --output data/graph.pt
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .go_annotations import parse_gaf
from .graph import build_graph
from .orthology import parse_orthology_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Build cross-species orthology graph")
    parser.add_argument("--orthology-path", required=True, help="Path to OrthoDB/edge-list TSV")
    parser.add_argument("--go-annotations-path", required=False, default=None, help="Path to GAF 2.2 file (.gaf or .gaf.gz)")
    parser.add_argument("--go-term", default="GO:0007049", help="GO term to use as binary label (default: GO:0007049 cell cycle)")
    parser.add_argument("--evidence-codes", default=None, help="Comma-separated evidence codes to keep, e.g. IDA,IMP")
    parser.add_argument("--taxon-filter", default=None, help="Only keep GAF rows matching this taxon substring")
    parser.add_argument("--target-taxon", default=None, help="Target species taxid (informational)")
    parser.add_argument("--output", required=True, help="Output path (.json or .pt)")
    args = parser.parse_args()

    orth_path = Path(args.orthology_path)
    if not orth_path.exists():
        parser.error(f"orthology file not found: {args.orthology_path}")
    if args.go_annotations_path and not Path(args.go_annotations_path).exists():
        parser.error(f"GO annotations file not found: {args.go_annotations_path}")
    out_path = Path(args.output)
    if out_path.parent != Path(".") and not out_path.parent.exists():
        try:
            out_path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            parser.error(f"cannot create output directory {out_path.parent}: {e}")

    print(f"[build_graph] Parsing orthology: {args.orthology_path}")
    try:
        edges = parse_orthology_file(args.orthology_path)
    except Exception as e:
        parser.error(f"failed to parse orthology file {args.orthology_path}: {e}")
    print(f"[build_graph] Found {len(edges)} ortholog edges")

    # Validate GO term format
    import re

    if not re.match(r"^GO:\d{7}$", args.go_term):
        parser.error(f"invalid GO term format {args.go_term!r}; expected GO:NNNNNNN")

    labels: dict[str, int] = {}
    if args.go_annotations_path:
        ec = set(s.strip() for s in args.evidence_codes.split(",")) if args.evidence_codes else None
        print(f"[build_graph] Parsing GO annotations: {args.go_annotations_path} (term={args.go_term})")
        try:
            anns = parse_gaf(args.go_annotations_path, go_term=args.go_term, evidence_codes=ec, taxon_filter=args.taxon_filter)
        except Exception as e:
            parser.error(f"failed to parse GAF file {args.go_annotations_path}: {e}")
        print(f"[build_graph] Found {len(anns)} annotations for {args.go_term}")
        # Positive gene set. GAF gene ids (e.g. SGD systematic ids like
        # "S000000001") are a DIFFERENT identifier namespace from this
        # module's orthology gene ids (OrthoDB's own "<taxid>_<n>:<seq>"
        # scheme) -- they do not match by string equality, and no
        # cross-reference/bridge between the two is implemented here.
        # Unioning unmatched annotation gene ids directly into the node set
        # (the previous behaviour) silently created isolated nodes with zero
        # orthology edges, so RWR/GNN propagation from them touches nothing
        # -- a graph that "has labels" but is structurally meaningless for
        # the cross-species task. Only genes already present in the real
        # orthology graph are labeled; annotations that don't resolve to an
        # orthology node are counted and reported, never silently added.
        positives = {a.gene_id for a in anns}
        all_genes = set()
        for e in edges:
            all_genes.add(e.gene_a)
            all_genes.add(e.gene_b)
        for g in all_genes:
            labels[g] = 1 if g in positives else 0
        num_pos = sum(labels.values())
        n_unresolved = len(positives - all_genes)
        print(f"[build_graph] Labels: {num_pos} positive / {len(labels)} total "
              f"(within the real orthology graph)")
        if n_unresolved:
            print(
                f"[build_graph] WARNING: {n_unresolved}/{len(positives)} GAF-annotated genes "
                "did not match any orthology-graph gene id (different identifier namespaces, "
                "no cross-reference bridge implemented) and were excluded rather than added "
                "as disconnected nodes. A real fix needs a GAF<->orthology gene-id bridge "
                "(e.g. via OrthoDB's own gene_xrefs UniProt/GOterm cross-references)."
            )
    else:
        # No GAF: all-zero labels (structure-only graph)
        print("[build_graph] No GAF provided — labels all zero (structure-only)")
        for e in edges:
            if e.gene_a not in labels:
                labels[e.gene_a] = 0
            if e.gene_b not in labels:
                labels[e.gene_b] = 0

    graph = build_graph(edges, labels)
    print(f"[build_graph] Graph: {graph.num_nodes} nodes, {graph.num_edges} edges, {len(graph.feature_names)} features")
    print(f"[build_graph] Features: {graph.feature_names}")
    graph.save(args.output)
    print(f"[build_graph] Saved to {args.output}")


if __name__ == "__main__":
    main()
