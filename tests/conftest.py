import tempfile
from pathlib import Path

import pytest

from data_pipeline.go_annotations import parse_gaf
from data_pipeline.graph import build_graph
from data_pipeline.orthology import parse_orthology_file


@pytest.fixture
def fixtures_dir():
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def orthology_edges(fixtures_dir):
    return parse_orthology_file(fixtures_dir / "orthology_edges.tsv")


@pytest.fixture
def go_annotations(fixtures_dir):
    return parse_gaf(fixtures_dir / "go_annotations.gaf")


@pytest.fixture
def synthetic_graph(orthology_edges, fixtures_dir):
    """Graph built from fixture orthology + GO annotations for GO:0007049."""
    gaf_path = fixtures_dir / "go_annotations.gaf"
    anns = parse_gaf(gaf_path, go_term="GO:0007049")
    positives = {a.gene_id for a in anns}
    # Also need all gene ids from orthology
    all_genes: set[str] = set()
    for e in orthology_edges:
        all_genes.add(e.gene_a)
        all_genes.add(e.gene_b)
    for a in anns:
        all_genes.add(a.gene_id)
    labels = {g: (1 if g in positives else 0) for g in all_genes}
    return build_graph(orthology_edges, labels)


@pytest.fixture
def synthetic_two_species_graph():
    """Explicitly synthetic two-species orthology graph with injected signal."""
    from data_pipeline.graph import generate_synthetic_graph

    return generate_synthetic_graph(
        num_source_genes=20,
        num_target_genes=15,
        num_ortholog_pairs=30,
        conservation_rate=0.85,
        seed=42,
    )


@pytest.fixture
def temp_graph_file(synthetic_graph):
    """Save synthetic_graph to a temp file for backend tests."""
    import json
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        obj = {
            "gene_ids": synthetic_graph.gene_ids,
            "species": synthetic_graph.species,
            "labels": synthetic_graph.labels,
            "edges": synthetic_graph.edges,
            "features": synthetic_graph.features,
            "feature_names": synthetic_graph.feature_names,
        }
        json.dump(obj, f)
        path = f.name
    yield path
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
