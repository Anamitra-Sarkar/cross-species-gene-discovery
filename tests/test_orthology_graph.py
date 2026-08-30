"""Tests: orthology parsing and graph construction correctness."""

from pathlib import Path

from data_pipeline.go_annotations import parse_gaf
from data_pipeline.graph import build_graph, generate_synthetic_graph
from data_pipeline.orthology import parse_orthology_file, parse_simplified_edge_list


def test_parse_simplified_edge_list(fixtures_dir=None):
    if fixtures_dir is None:
        fixtures_dir = Path(__file__).parent / "fixtures"
    edges = parse_simplified_edge_list(fixtures_dir / "orthology_edges.tsv")
    assert len(edges) == 14, f"Expected 14 edges from fixture, got {len(edges)}"
    # Check confidence scores are in [0,1]
    for e in edges:
        assert 0.0 <= e.score <= 1.0, f"Score out of range: {e.score}"
    # Check species are expected
    species = {e.species_a for e in edges} | {e.species_b for e in edges}
    assert "559292" in species
    assert "4896" in species


def test_parse_orthology_auto_detect(fixtures_dir=None):
    if fixtures_dir is None:
        fixtures_dir = Path(__file__).parent / "fixtures"
    edges = parse_orthology_file(fixtures_dir / "orthology_edges.tsv")
    assert len(edges) >= 10


def test_graph_construction_node_edge_counts(synthetic_graph):
    g = synthetic_graph
    # Nodes: genes from orthology + labels (SGD + PomBase + MGI)
    assert g.num_nodes >= 20, f"Expected >=20 nodes, got {g.num_nodes}"
    # Edges: fixture has 14 edges, deduped undirected
    assert g.num_edges == 14, f"Expected 14 edges, got {g.num_edges}"
    # All gene_ids are sorted
    assert g.gene_ids == sorted(g.gene_ids)
    # gene_to_idx consistent
    for gid, idx in g.gene_to_idx.items():
        assert g.gene_ids[idx] == gid


def test_graph_weights_preserved(synthetic_graph):
    g = synthetic_graph
    # Check that edge weights are preserved from fixture
    weights = [w for _, _, w in g.edges]
    assert all(0.0 < w <= 1.0 for w in weights)
    # Known high-confidence edge
    # Find edge involving SGD:S000000001
    found = False
    for s, t, w in g.edges:
        if "S000000001" in g.gene_ids[s] or "S000000001" in g.gene_ids[t]:
            found = True
            break
    assert found, "Expected edge involving S000000001"


def test_graph_features_shape(synthetic_graph):
    g = synthetic_graph
    assert g.features is not None
    assert len(g.features) == g.num_nodes
    assert len(g.feature_names) >= 4  # degree_norm, mean_conf, max_conf + species one-hots
    assert "degree_norm" in g.feature_names
    assert "mean_confidence" in g.feature_names
    # Feature dims consistent
    dim = len(g.features[0])
    for row in g.features:
        assert len(row) == dim


def test_graph_undirected_symmetry(synthetic_graph):
    g = synthetic_graph
    # Edges should be stored undirected with small idx first
    for s, t, _ in g.edges:
        assert s < t, f"Edge not normalized: ({s}, {t})"
    # No duplicate edges
    seen = set()
    for s, t, _ in g.edges:
        assert (s, t) not in seen, f"Duplicate edge ({s}, {t})"
        seen.add((s, t))


def test_synthetic_graph_generation():
    g = generate_synthetic_graph(num_source_genes=10, num_target_genes=8, num_ortholog_pairs=12, conservation_rate=0.9, seed=123)
    assert g.num_nodes == 18
    # Should have <= num_ortholog_pairs edges (some may be duplicates)
    assert g.num_edges <= 12
    assert g.num_edges >= 5
    # Species should be two values
    assert set(g.species) == {"559292", "4896"}
    # Labels should have both classes
    assert 0 in g.labels.values()
    assert 1 in g.labels.values()


def test_gaf_parsing(fixtures_dir=None):
    if fixtures_dir is None:
        fixtures_dir = Path(__file__).parent / "fixtures"
    anns = parse_gaf(fixtures_dir / "go_annotations.gaf", go_term="GO:0007049")
    # Fixture has ~11 positives for GO:0007049 across SGD+PomBase+MGI
    assert len(anns) >= 10, f"Expected >=10 annotations for GO:0007049, got {len(anns)}"
    for a in anns:
        assert a.go_id == "GO:0007049"

    # Taxon filter
    yeast_anns = parse_gaf(fixtures_dir / "go_annotations.gaf", go_term="GO:0007049", taxon_filter="559292")
    assert len(yeast_anns) >= 4
    for a in yeast_anns:
        assert "559292" in a.taxon
