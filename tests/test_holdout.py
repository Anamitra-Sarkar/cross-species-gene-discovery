"""Tests: cross-species holdout split correctness (no leakage) and evaluation."""

import numpy as np
import pytest

from data_pipeline.graph import generate_synthetic_graph
from model.evaluation import cross_species_holdout_split, evaluate_rwr_holdout
from model.rwr import rwr_predict


def test_holdout_no_leakage(synthetic_two_species_graph):
    """Ensure training never uses target test labels."""
    g = synthetic_two_species_graph
    source_species = "559292"
    target_species = "4896"

    source_idx, _, target_test_idx = cross_species_holdout_split(
        g, source_species, target_species, test_fraction=1.0, seed=42
    )

    # Source and target sets are disjoint
    assert len(set(source_idx) & set(target_test_idx)) == 0
    # All source belong to source species
    for idx in source_idx:
        assert g.species[idx] == source_species
    # All target test belong to target species
    for idx in target_test_idx:
        assert g.species[idx] == target_species
    # RWR seeds are only source positives
    source_pos_genes = [g.gene_ids[i] for i in source_idx if g.labels[g.gene_ids[i]] == 1]
    # None of the target test genes should be seeds
    for idx in target_test_idx:
        assert g.gene_ids[idx] not in source_pos_genes


def test_holdout_all_target_genes_covered(synthetic_two_species_graph):
    g = synthetic_two_species_graph
    source_idx, hidden_idx, target_test_idx = cross_species_holdout_split(
        g, "559292", "4896", test_fraction=1.0, seed=0
    )
    all_target = [i for i, s in enumerate(g.species) if s == "4896"]
    # With test_fraction=1.0, all target genes should be in test
    assert set(target_test_idx) == set(all_target)
    assert len(hidden_idx) == 0


def test_holdout_rwr_evaluation_runs(synthetic_two_species_graph):
    g = synthetic_two_species_graph
    result = evaluate_rwr_holdout(g, source_species="559292", target_species="4896", alpha=0.3)
    assert 0.0 <= result.auroc <= 1.0
    assert 0.0 <= result.auprc <= 1.0
    assert result.num_test > 0
    assert result.method.startswith("RWR")
    # Baseline should be 0.5 for AUROC
    assert result.baseline_auroc == 0.5
    assert 0.0 <= result.baseline_auprc <= 1.0


def test_injected_signal_recovered():
    """Synthetic graph with high conservation rate should give AUROC > 0.5 (better than random)."""
    g = generate_synthetic_graph(
        num_source_genes=30,
        num_target_genes=20,
        num_ortholog_pairs=50,
        conservation_rate=0.9,  # high conservation
        seed=0,
    )
    result = evaluate_rwr_holdout(g, source_species="559292", target_species="4896", alpha=0.3)
    # With 90% conservation, RWR should beat random
    assert result.auroc > 0.5, f"Expected AUROC > 0.5 with high conservation, got {result.auroc:.3f}"


def test_rwr_uses_only_source_labels(synthetic_two_species_graph):
    """Verify RWR prediction would differ if target positives were incorrectly included as seeds."""
    g = synthetic_two_species_graph
    source_idx = [i for i, s in enumerate(g.species) if s == "559292"]
    target_idx = [i for i, s in enumerate(g.species) if s == "4896"]

    source_pos = [g.gene_ids[i] for i in source_idx if g.labels[g.gene_ids[i]] == 1]
    target_pos = [g.gene_ids[i] for i in target_idx if g.labels[g.gene_ids[i]] == 1]

    scores_source_only = rwr_predict(g, source_pos, alpha=0.3)
    scores_with_leakage = rwr_predict(g, source_pos + target_pos, alpha=0.3)

    # They should differ (leakage changes scores)
    # If no target positives, they would be equal — handle that case
    if target_pos:
        assert not np.allclose(scores_source_only, scores_with_leakage), "Leakage should change RWR scores"
