"""Tests: RWR / label propagation correctness."""

import numpy as np

from data_pipeline.graph import generate_synthetic_graph
from model.rwr import build_transition_matrix, random_walk_with_restart, rwr_predict


def test_rwr_convergence(synthetic_graph):
    g = synthetic_graph
    scores = random_walk_with_restart(g, alpha=0.3)
    # Scores are non-negative, bounded; sum <=1 (mass may leak to isolated nodes with zero column sum)
    # For connected component with positives, sum close to 1; with isolated nodes, sum <1
    assert len(scores) == g.num_nodes
    assert all(s >= 0 for s in scores)
    assert 0.0 < scores.sum() <= 1.0 + 1e-6
    if scores.sum() > 0:
        # At least some mass retained
        assert scores.sum() > 0.1


def test_rwr_score_ordering_closer_orthologs_higher():
    """Create a chain: source_pos -- target_close (high weight) -- target_far (low weight)
    Closer ortholog should score higher than distant one."""
    from data_pipeline.graph import CrossSpeciesGraph
    from data_pipeline.orthology import OrthologEdge

    # Small graph: S1(positive) -- T1 (0.95) -- T2 (0.1)  (T1 connected to S1, T2 connected via T1? no direct)
    # Actually make: S1 -- T1 (0.95), S1 -- T2 (0.1)  both connected to source
    edges = [
        OrthologEdge("OG1", "SRC_G1", "559292", "TGT_G1", "4896", 0.95),
        OrthologEdge("OG2", "SRC_G1", "559292", "TGT_G2", "4896", 0.10),
    ]
    labels = {"SRC_G1": 1, "TGT_G1": 0, "TGT_G2": 0}
    from data_pipeline.graph import build_graph

    g = build_graph(edges, labels)
    # RWR from SRC_G1
    scores = rwr_predict(g, ["SRC_G1"], alpha=0.3)
    idx_t1 = g.gene_to_idx["TGT_G1"]
    idx_t2 = g.gene_to_idx["TGT_G2"]
    # Weighted RWR should rank higher-confidence ortholog higher
    assert scores[idx_t1] > scores[idx_t2], f"Expected TGT_G1 ({scores[idx_t1]:.4f}) > TGT_G2 ({scores[idx_t2]:.4f}) weighted"


def test_rwr_no_positives_returns_zeros(synthetic_graph):
    from data_pipeline.graph import CrossSpeciesGraph

    # Zero out all labels
    g2 = CrossSpeciesGraph(
        gene_ids=synthetic_graph.gene_ids,
        gene_to_idx=synthetic_graph.gene_to_idx,
        species=synthetic_graph.species,
        labels={gid: 0 for gid in synthetic_graph.gene_ids},
        edges=synthetic_graph.edges,
        features=synthetic_graph.features,
        feature_names=synthetic_graph.feature_names,
    )
    scores = random_walk_with_restart(g2, alpha=0.3)
    assert np.allclose(scores, 0.0)


def test_rwr_isolated_nodes():
    """Isolated target gene (no edges) should get zero score if not a seed."""
    from data_pipeline.graph import CrossSpeciesGraph, build_graph
    from data_pipeline.orthology import OrthologEdge

    edges = [
        OrthologEdge("OG1", "SRC_G1", "559292", "TGT_G1", "4896", 0.9),
    ]
    labels = {"SRC_G1": 1, "TGT_G1": 0, "ISOLATED_G": 0}
    g = build_graph(edges, labels, include_isolated_genes=["ISOLATED_G"])
    scores = rwr_predict(g, ["SRC_G1"], alpha=0.3)
    idx_iso = g.gene_to_idx["ISOLATED_G"]
    idx_tgt = g.gene_to_idx["TGT_G1"]
    assert scores[idx_iso] == 0.0 or scores[idx_iso] < scores[idx_tgt]


def test_transition_matrix_column_stochastic(synthetic_graph):
    W = build_transition_matrix(synthetic_graph, weighted=True)
    # Column sums should be 0 or 1 (isolated nodes 0)
    col_sums = np.array(W.sum(axis=0)).flatten()
    for cs in col_sums:
        assert cs == pytest.approx(0.0, abs=1e-6) or cs == pytest.approx(1.0, abs=1e-6), f"Column sum {cs} not 0 or 1"


def test_rwr_repeatability(synthetic_two_species_graph):
    g = synthetic_two_species_graph
    s1 = random_walk_with_restart(g, alpha=0.3)
    s2 = random_walk_with_restart(g, alpha=0.3)
    np.testing.assert_allclose(s1, s2)


# Need pytest import
import pytest
