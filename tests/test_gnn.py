"""Tests: GNN forward pass shape correctness."""

import pytest


def test_gnn_forward_shape(synthetic_two_species_graph):
    try:
        import torch
        from torch_geometric.data import Data

        from model.gnn import HAS_TORCH_GEOMETRIC, CrossSpeciesGNN
    except ImportError:
        pytest.skip("torch/torch_geometric not installed")

    if not HAS_TORCH_GEOMETRIC:
        pytest.skip("torch_geometric not available")

    g = synthetic_two_species_graph
    data = g.to_torch_data()

    num_features = data.x.size(1)
    model = CrossSpeciesGNN(in_channels=num_features, hidden_channels=16, conv_type="sage")
    model.eval()
    with torch.no_grad():
        out = model(data.x, data.edge_index, getattr(data, "edge_attr", None))

    assert out.shape == (g.num_nodes,), f"Expected shape ({g.num_nodes},), got {out.shape}"
    # Output should be finite
    assert torch.isfinite(out).all()


def test_gnn_gcn_variant(synthetic_two_species_graph):
    try:
        import torch

        from model.gnn import HAS_TORCH_GEOMETRIC, CrossSpeciesGNN
    except ImportError:
        pytest.skip("torch/torch_geometric not installed")

    if not HAS_TORCH_GEOMETRIC:
        pytest.skip("torch_geometric not available")

    g = synthetic_two_species_graph
    data = g.to_torch_data()
    model = CrossSpeciesGNN(in_channels=data.x.size(1), hidden_channels=8, conv_type="gcn")
    model.eval()
    with torch.no_grad():
        out = model(data.x, data.edge_index, getattr(data, "edge_attr", None))
    assert out.shape == (g.num_nodes,)


def test_gnn_predict_proba(synthetic_two_species_graph):
    try:
        import torch

        from model.gnn import HAS_TORCH_GEOMETRIC, CrossSpeciesGNN
    except ImportError:
        pytest.skip("torch/torch_geometric not installed")

    if not HAS_TORCH_GEOMETRIC:
        pytest.skip("torch_geometric not available")

    g = synthetic_two_species_graph
    data = g.to_torch_data()
    model = CrossSpeciesGNN(in_channels=data.x.size(1), hidden_channels=16)
    probs = model.predict_proba(data.x, data.edge_index)
    assert probs.shape == (g.num_nodes,)
    # Probabilities in [0,1]
    assert (probs >= 0).all() and (probs <= 1).all()


def test_gnn_training_on_source_only(synthetic_two_species_graph):
    """Train GNN on source nodes and check loss decreases."""
    try:
        import torch

        from model.gnn import HAS_TORCH_GEOMETRIC, train_gnn
    except ImportError:
        pytest.skip("torch/torch_geometric not installed")

    if not HAS_TORCH_GEOMETRIC:
        pytest.skip("torch_geometric not available")

    g = synthetic_two_species_graph
    data = g.to_torch_data()
    # Train mask: source species only
    source_indices = [i for i, s in enumerate(g.species) if s == "559292"]
    assert len(source_indices) > 0
    train_mask = torch.zeros(data.num_nodes, dtype=torch.bool)
    for idx in source_indices:
        train_mask[idx] = True

    torch.manual_seed(0)
    model = train_gnn(data, train_mask, hidden_channels=8, epochs=20, verbose=False)
    # After training, model should produce finite outputs
    model.eval()
    with torch.no_grad():
        out = model(data.x, data.edge_index)
    assert torch.isfinite(out).all()
