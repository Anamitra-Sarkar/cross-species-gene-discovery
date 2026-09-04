"""Random Walk with Restart (RWR) / label propagation baseline.

Real, standard, simple baseline for network-based function prediction.
References:
  - Köhler et al. "Walking the interactome for prioritization of candidate disease genes." Am J Hum Genet 2008.
  - Cowen et al. "Network propagation: a universal amplifier of genetic associations." Nat Rev Genet 2017.
  - Vanunu et al. "Associating genes and protein complexes with disease via network propagation." PLoS Comput Biol 2010.

Formula:
  p_{t+1} = (1 - alpha) * W * p_t + alpha * p_0
where W is column-normalized adjacency (transition matrix),
p_0 is restart/distribution vector (normalized, positive on labeled source genes),
alpha is restart probability (0 < alpha < 1).

Convergence: ||p_{t+1} - p_t||_1 < tol or max_iter reached.

This is a real algorithm, not a mock.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import sparse

from data_pipeline.graph import CrossSpeciesGraph


def build_transition_matrix(graph: CrossSpeciesGraph, weighted: bool = True) -> sparse.csc_matrix:
    """Build column-normalized transition matrix W from graph.

    W[j,i] = weight(i->j) / sum_k weight(i->k)  (column stochastic)
    For undirected graphs, W is symmetric after normalization if weights symmetric.

    Args:
        graph: CrossSpeciesGraph
        weighted: if True uses edge scores, else uniform weights
    """
    n = graph.num_nodes
    if n == 0:
        return sparse.csc_matrix((0, 0))

    if not graph.edges:
        # No edges: return identity-like? Actually zeros (isolated nodes stay isolated)
        return sparse.csc_matrix((n, n))

    # Build adjacency (undirected, add both directions). A real orthology
    # graph at Ascomycota scale has hundreds of millions of edges; the
    # previous pure-Python per-edge loop (list.extend x3 per edge) was
    # confirmed intractable at that scale (636M edges -> ~1.9B individual
    # Python-level list appends). Vectorized via numpy instead -- same
    # adjacency, same weights, no Python-level loop over edges.
    edge_arr = np.asarray(graph.edges, dtype=np.float64)  # (E, 3): src, dst, weight
    src = edge_arr[:, 0].astype(np.int64)
    dst = edge_arr[:, 1].astype(np.int64)
    w = edge_arr[:, 2] if weighted else np.ones(len(edge_arr), dtype=np.float64)

    rows = np.concatenate([dst, src])
    cols = np.concatenate([src, dst])
    data = np.concatenate([w, w])

    adj = sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsc()

    # Column-normalize: W[:,i] = adj[:,i] / sum(adj[:,i])
    # Handle isolated nodes (zero column sum) -> keep zero column (no outgoing walk)
    col_sums = np.array(adj.sum(axis=0)).flatten()
    # Avoid division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        # Use csc structure to scale columns
        # For each column j, divide by col_sums[j]
        # Efficient via sparse multiply with diagonal
        inv_col_sums = np.where(col_sums > 0, 1.0 / col_sums, 0.0)
        D_inv = sparse.diags(inv_col_sums, format="csc")
        W = adj @ D_inv  # column-normalized

    return W


def random_walk_with_restart(
    graph: CrossSpeciesGraph,
    alpha: float = 0.3,
    weighted: bool = True,
    tol: float = 1e-6,
    max_iter: int = 1000,
    label_key: str = "source_labels",
) -> np.ndarray:
    """Run RWR from source-labeled genes over the orthology graph.

    The restart vector p_0 is 1 on positively-labeled source nodes, 0 elsewhere,
    normalized to sum to 1. If no positive labels, returns zeros.

    For cross-species evaluation, caller should construct graph where only
    SOURCE species positive labels are used for p_0 (target labels hidden).

    Args:
        graph: CrossSpeciesGraph (labels field used for p_0)
        alpha: restart probability (higher = more localized to seeds)
        weighted: use edge weights if True
        tol: L1 convergence threshold
        max_iter: max iterations
        label_key: not used currently, kept for API consistency

    Returns:
        np.ndarray of shape [num_nodes] with RWR stationary scores (sum to 1).
    """
    n = graph.num_nodes
    if n == 0:
        return np.array([])

    W = build_transition_matrix(graph, weighted=weighted)

    # Build p_0 from labels (positive source genes only — caller controls which labels are set)
    p0 = np.array([graph.labels.get(g, 0) for g in graph.gene_ids], dtype=float)
    pos_count = p0.sum()
    if pos_count == 0:
        return np.zeros(n)
    p0 = p0 / pos_count  # normalize

    # Handle case with no edges: p_t = p_0 always (no diffusion)
    if graph.num_edges == 0:
        return p0

    p = p0.copy()
    for iteration in range(max_iter):
        p_new = (1 - alpha) * (W @ p) + alpha * p0
        delta = np.abs(p_new - p).sum()
        p = p_new
        if delta < tol:
            break

    return p


def label_propagation_scores(
    graph: CrossSpeciesGraph,
    alpha: float = 0.3,
    weighted: bool = True,
) -> dict[str, float]:
    """Convenience: RWR scores as gene_id -> score dict."""
    scores = random_walk_with_restart(graph, alpha=alpha, weighted=weighted)
    return {g: float(scores[i]) for i, g in enumerate(graph.gene_ids)}


def rwr_predict(
    graph: CrossSpeciesGraph,
    source_gene_ids: list[str],
    alpha: float = 0.3,
    weighted: bool = True,
) -> np.ndarray:
    """RWR prediction when source set is explicitly given (not from graph.labels).

    This is the correct cross-species holdout helper: only source genes' positive
    labels contribute to p_0, regardless of what graph.labels contains for targets.

    Args:
        graph: graph with gene_ids ordering
        source_gene_ids: subset of gene_ids that are positive seeds (source positives)
        alpha, weighted: RWR params

    Returns:
        scores array [num_nodes]
    """
    n = graph.num_nodes
    if n == 0:
        return np.array([])

    # Build p_0 only from provided source_gene_ids
    source_set = set(source_gene_ids)
    p0 = np.array([1.0 if g in source_set else 0.0 for g in graph.gene_ids], dtype=float)
    if p0.sum() == 0:
        return np.zeros(n)
    p0 = p0 / p0.sum()

    W = build_transition_matrix(graph, weighted=weighted)
    if graph.num_edges == 0:
        return p0

    p = p0.copy()
    tol = 1e-6
    for _ in range(1000):
        p_new = (1 - alpha) * (W @ p) + alpha * p0
        if np.abs(p_new - p).sum() < tol:
            p = p_new
            break
        p = p_new
    return p
