"""Cross-species holdout evaluation.

The correct, honest protocol:
  train/tune using only source-species (annotated) genes and graph structure,
  then evaluate on held-out TARGET-species genes whose real GO annotations are known
  (but hidden during training).

Reports AUROC / AUPRC vs naive majority-class baseline.
Supports synthetic two-species graph with injected conservation signal as explicitly
documented synthetic verification (not a real biological finding).

Metrics use scikit-learn if available, otherwise pure numpy fallback.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

try:
    from sklearn.metrics import roc_auc_score, average_precision_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

from data_pipeline.graph import CrossSpeciesGraph
from model.rwr import rwr_predict


def _auroc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute AUROC (requires both classes present; return 0.5 if not)."""
    if len(np.unique(y_true)) < 2:
        return 0.5
    if HAS_SKLEARN:
        return float(roc_auc_score(y_true, y_score))
    # Fallback: Mann-Whitney U / ranking method
    # AUROC = P(score(pos) > score(neg))
    pos_scores = y_score[y_true == 1]
    neg_scores = y_score[y_true == 0]
    if len(pos_scores) == 0 or len(neg_scores) == 0:
        return 0.5
    # Pairwise comparison (O(n*m) — fine for small graphs)
    count = 0.0
    total = len(pos_scores) * len(neg_scores)
    for ps in pos_scores:
        for ns in neg_scores:
            if ps > ns:
                count += 1
            elif ps == ns:
                count += 0.5
    return count / total if total > 0 else 0.5


def _auprc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute AUPRC (average precision). Baseline = positive prevalence."""
    if len(np.unique(y_true)) < 2:
        # If only one class, AUPRC is prevalence or 0.5
        return float(y_true.mean()) if len(y_true) > 0 else 0.5
    if HAS_SKLEARN:
        return float(average_precision_score(y_true, y_score))
    # Fallback: compute average precision via sorted scores
    order = np.argsort(-y_score)
    y_sorted = y_true[order]
    precisions = []
    tp = 0
    for i, label in enumerate(y_sorted):
        if label == 1:
            tp += 1
            precisions.append(tp / (i + 1))
    return float(np.mean(precisions)) if precisions else 0.0


@dataclass
class EvaluationResult:
    auroc: float
    auprc: float
    baseline_auroc: float
    baseline_auprc: float
    num_test: int
    num_positive_test: int
    method: str


def naive_baseline_scores(y_true: np.ndarray, prevalence: float) -> np.ndarray:
    """Naive most-common / prevalence baseline: constant score = prevalence."""
    return np.full_like(y_true, fill_value=prevalence, dtype=float)


def cross_species_holdout_split(
    graph: CrossSpeciesGraph,
    source_species: str,
    target_species: str,
    test_fraction: float = 0.5,
    seed: int = 42,
) -> tuple[list[int], list[int], list[int]]:
    """Split target-species genes into test set; return indices.

    No leakage: training uses only source-species genes + graph structure.

    Returns:
        (source_indices, target_train_hidden_indices, target_test_indices)
        For honest evaluation, target_train_hidden is empty (all target labels hidden).
        We provide it for pipeline verification (if you want to simulate semi-supervised).
        Standard protocol: train on source only, test on all target held-out.

    Args:
        graph: CrossSpeciesGraph
        source_species: taxid string for source (e.g. "559292")
        target_species: taxid string for target (e.g. "4896")
        test_fraction: fraction of target genes to hold out as test (default 0.5 = all evaluated)
        seed: random seed for subsampling test set
    """
    rng = random.Random(seed)
    target_indices = [i for i, s in enumerate(graph.species) if s == target_species]
    source_indices = [i for i, s in enumerate(graph.species) if s == source_species]

    # Shuffle and split target
    rng.shuffle(target_indices)
    n_test = max(1, int(len(target_indices) * test_fraction))
    target_test = target_indices[:n_test]
    target_train_hidden = target_indices[n_test:]

    return source_indices, target_train_hidden, target_test


def evaluate_rwr_holdout(
    graph: CrossSpeciesGraph,
    source_species: str,
    target_species: str,
    alpha: float = 0.3,
    test_fraction: float = 1.0,
    seed: int = 42,
) -> EvaluationResult:
    """Evaluate RWR with cross-species holdout (train on source only, test on target).

    Args:
        graph: CrossSpeciesGraph with labels for all genes (target labels are ground truth for scoring)
        source_species, target_species: taxids
        alpha: RWR restart prob
        test_fraction: fraction of target genes used as test set (1.0 = all target genes)
        seed: random seed

    Returns:
        EvaluationResult with AUROC/AUPRC vs baselines
    """
    source_indices, _, target_test = cross_species_holdout_split(
        graph, source_species, target_species, test_fraction=test_fraction, seed=seed
    )

    # Source positives are seeds for RWR
    source_pos_genes = [graph.gene_ids[i] for i in source_indices if graph.labels.get(graph.gene_ids[i], 0) == 1]

    # Run RWR with only source positives as seeds (no target info)
    scores = rwr_predict(graph, source_pos_genes, alpha=alpha)

    # Evaluate on target test set
    y_true = np.array([graph.labels.get(graph.gene_ids[i], 0) for i in target_test], dtype=int)
    y_score = np.array([scores[i] for i in target_test], dtype=float)

    auroc = _auroc(y_true, y_score)
    auprc = _auprc(y_true, y_score)

    # Baselines
    prevalence = float(y_true.mean()) if len(y_true) > 0 else 0.5
    # AUROC baseline is 0.5 (random), AUPRC baseline is prevalence
    baseline_auroc = 0.5
    baseline_auprc = prevalence

    return EvaluationResult(
        auroc=auroc,
        auprc=auprc,
        baseline_auroc=baseline_auroc,
        baseline_auprc=baseline_auprc,
        num_test=len(y_true),
        num_positive_test=int(y_true.sum()),
        method=f"RWR(alpha={alpha})",
    )


def evaluate_gnn_holdout(
    graph: CrossSpeciesGraph,
    source_species: str,
    target_species: str,
    hidden_channels: int = 16,
    epochs: int = 80,
    test_fraction: float = 1.0,
    seed: int = 42,
) -> EvaluationResult:
    """Evaluate GNN with cross-species holdout.

    Trains GNN only on source-species nodes, tests on target holdout.

    Requires torch + torch_geometric. If not installed, raises ImportError.
    """
    try:
        import torch
        from torch_geometric.data import Data  # noqa: F401

        from model.gnn import HAS_TORCH_GEOMETRIC, train_gnn

        if not HAS_TORCH_GEOMETRIC:
            raise ImportError
    except ImportError:
        raise ImportError("torch and torch_geometric required for GNN evaluation")

    source_indices, _, target_test = cross_species_holdout_split(
        graph, source_species, target_species, test_fraction=test_fraction, seed=seed
    )

    data = graph.to_torch_data()

    # Train mask: only source nodes
    train_mask = torch.zeros(data.num_nodes, dtype=torch.bool)
    for idx in source_indices:
        train_mask[idx] = True

    # Set seeds for reproducibility
    torch.manual_seed(seed)
    np.random.seed(seed)

    model = train_gnn(data, train_mask, hidden_channels=hidden_channels, epochs=epochs, verbose=False)

    # Predict on all nodes
    model.eval()
    with torch.no_grad():
        logits = model(data.x, data.edge_index, getattr(data, "edge_attr", None))
        probs = torch.sigmoid(logits).cpu().numpy()

    y_true = np.array([graph.labels.get(graph.gene_ids[i], 0) for i in target_test], dtype=int)
    y_score = np.array([probs[i] for i in target_test], dtype=float)

    auroc = _auroc(y_true, y_score)
    auprc = _auprc(y_true, y_score)
    prevalence = float(y_true.mean()) if len(y_true) > 0 else 0.5

    return EvaluationResult(
        auroc=auroc,
        auprc=auprc,
        baseline_auroc=0.5,
        baseline_auprc=prevalence,
        num_test=len(y_true),
        num_positive_test=int(y_true.sum()),
        method=f"GNN(hidden={hidden_channels},epochs={epochs})",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Cross-species holdout evaluation")
    parser.add_argument("--graph", required=True, help="Path to graph JSON/PT file")
    parser.add_argument("--source-taxon", default="559292")
    parser.add_argument("--target-taxon", default="4896")
    parser.add_argument("--alpha", type=float, default=0.3)
    parser.add_argument("--target-taxon-param", default=None, help="Alias for --target-taxon")
    args = parser.parse_args()

    target_taxon = args.target_taxon_param or args.target_taxon
    # Support both arg names
    if hasattr(args, "target_taxon"):
        target_taxon = args.target_taxon

    graph = CrossSpeciesGraph.load(args.graph)
    print(f"Graph: {graph.num_nodes} nodes, {graph.num_edges} edges")
    print(f"Species: {set(graph.species)}")

    result = evaluate_rwr_holdout(graph, args.source_taxon, target_taxon, alpha=args.alpha)
    print(json.dumps(result.__dict__, indent=2))
    print(f"RWR AUROC={result.auroc:.3f} (baseline {result.baseline_auroc:.3f}), AUPRC={result.auprc:.3f} (baseline {result.baseline_auprc:.3f})")


if __name__ == "__main__":
    main()
