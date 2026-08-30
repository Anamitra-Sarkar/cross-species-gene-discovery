"""Cross-species orthology graph construction.

Nodes = genes (with species/taxon attribute).
Edges = ortholog pairs (undirected, weighted by orthology confidence).

Node features (real, simple, computable without sequence embeddings):
  - degree (normalized)
  - mean orthology confidence (mean of incident edge scores)
  - max orthology confidence
  - species one-hot (per unique species in graph)
  (No fabricated sequence embeddings; these are documented, real, simple features.)

Graph is represented as:
  - Lightweight Python dict/edge-list for pipeline without torch
  - torch_geometric.data.Data when torch is available (for GNN)

Also provides helpers for synthetic fixture generation (explicitly synthetic).
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .orthology import OrthologEdge


@dataclass
class CrossSpeciesGraph:
    """Lightweight graph representation (no torch dependency for pipeline)."""

    gene_ids: list[str]  # ordered node list
    gene_to_idx: dict[str, int]
    species: list[str]  # per node, taxon string
    labels: dict[str, int]  # gene_id -> 0/1 (binary GO term)
    edges: list[tuple[int, int, float]]  # (src_idx, dst_idx, weight) undirected
    features: Optional[list[list[float]]] = None  # [num_nodes, num_features]
    feature_names: list[str] = field(default_factory=list)

    @property
    def num_nodes(self) -> int:
        return len(self.gene_ids)

    @property
    def num_edges(self) -> int:
        return len(self.edges)

    def to_torch_data(self):
        """Convert to torch_geometric.data.Data (requires torch + torch_geometric)."""
        try:
            import torch
            from torch_geometric.data import Data
        except ImportError as e:
            raise ImportError("torch and torch_geometric required for to_torch_data()") from e

        # Edge index: make bidirectional (undirected graph)
        edge_list = []
        edge_weights = []
        for s, t, w in self.edges:
            edge_list.append((s, t))
            edge_list.append((t, s))
            edge_weights.append(w)
            edge_weights.append(w)

        if edge_list:
            edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
            edge_weight = torch.tensor(edge_weights, dtype=torch.float)
        else:
            edge_index = torch.empty((2, 0), dtype=torch.long)
            edge_weight = torch.empty((0,), dtype=torch.float)

        # Node features
        if self.features is not None:
            x = torch.tensor(self.features, dtype=torch.float)
        else:
            # Fallback: identity / degree-based
            x = torch.eye(len(self.gene_ids), dtype=torch.float)

        # Labels
        y_list = [self.labels.get(g, 0) for g in self.gene_ids]
        y = torch.tensor(y_list, dtype=torch.float)

        data = Data(x=x, edge_index=edge_index, edge_attr=edge_weight, y=y)
        data.gene_ids = self.gene_ids
        data.species = self.species
        return data

    def save(self, path: str | Path) -> None:
        """Save graph to JSON (lightweight) or .pt (torch)."""
        path = Path(path)
        if path.suffix == ".pt":
            data = self.to_torch_data()
            try:
                import torch

                torch.save(data, path)
            except ImportError:
                # Fallback to JSON
                path = path.with_suffix(".json")
                self._save_json(path)
        else:
            self._save_json(path)

    def _save_json(self, path: Path) -> None:
        obj = {
            "gene_ids": self.gene_ids,
            "species": self.species,
            "labels": self.labels,
            "edges": self.edges,
            "features": self.features,
            "feature_names": self.feature_names,
        }
        path.write_text(json.dumps(obj, indent=2))

    @staticmethod
    def load(path: str | Path) -> "CrossSpeciesGraph":
        path = Path(path)
        if path.suffix == ".pt":
            try:
                import torch

                data = torch.load(path, weights_only=False)
                # Reconstruct CrossSpeciesGraph from Data
                gene_ids = getattr(data, "gene_ids", [f"gene_{i}" for i in range(data.num_nodes)])
                species = getattr(data, "species", ["unknown"] * data.num_nodes)
                # Edges: deduplicate undirected
                edge_index = data.edge_index.numpy() if hasattr(data.edge_index, "numpy") else []
                edges = []
                seen = set()
                for i in range(edge_index.shape[1] if len(edge_index) else 0):
                    s, t = int(edge_index[0, i]), int(edge_index[1, i])
                    if s > t:
                        s, t = t, s
                    if (s, t) in seen:
                        continue
                    seen.add((s, t))
                    w = float(data.edge_attr[i]) if hasattr(data, "edge_attr") and data.edge_attr is not None else 1.0
                    edges.append((s, t, w))
                y = data.y.tolist() if hasattr(data, "y") else [0] * len(gene_ids)
                labels = {g: int(y[i]) for i, g in enumerate(gene_ids)}
                features = data.x.tolist() if hasattr(data, "x") else None
                return CrossSpeciesGraph(
                    gene_ids=gene_ids,
                    gene_to_idx={g: i for i, g in enumerate(gene_ids)},
                    species=species,
                    labels=labels,
                    edges=edges,
                    features=features,
                )
            except ImportError:
                pass
        # JSON fallback
        obj = json.loads(path.read_text())
        return CrossSpeciesGraph(
            gene_ids=obj["gene_ids"],
            gene_to_idx={g: i for i, g in enumerate(obj["gene_ids"])},
            species=obj["species"],
            labels={k: int(v) for k, v in obj["labels"].items()},
            edges=[tuple(e) for e in obj["edges"]],
            features=obj.get("features"),
            feature_names=obj.get("feature_names", []),
        )


def build_graph(
    ortholog_edges: list[OrthologEdge],
    labels: dict[str, int],
    include_isolated_genes: Optional[list[str]] = None,
) -> CrossSpeciesGraph:
    """Build cross-species graph from ortholog edges + labels.

    Args:
        ortholog_edges: list of OrthologEdge
        labels: gene_id -> 0/1
        include_isolated_genes: extra genes with no edges (still become nodes)
    """
    # Collect all genes
    gene_species: dict[str, str] = {}
    for e in ortholog_edges:
        if e.gene_a not in gene_species:
            gene_species[e.gene_a] = e.species_a
        if e.gene_b not in gene_species:
            gene_species[e.gene_b] = e.species_b

    # Add labels-only genes (genes that appear in GAF but not orthology)
    for g in labels:
        if g not in gene_species:
            gene_species[g] = "unknown"
    # Add isolated genes
    if include_isolated_genes:
        for g in include_isolated_genes:
            if g not in gene_species:
                gene_species[g] = "unknown"

    # Stable ordering: sorted by gene_id
    gene_ids = sorted(gene_species.keys())
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    species = [gene_species[g] for g in gene_ids]

    # Deduplicate edges (undirected)
    seen: set[tuple[int, int]] = set()
    edges: list[tuple[int, int, float]] = []
    for e in ortholog_edges:
        # Ensure both genes are in mapping (they are)
        a = gene_to_idx.get(e.gene_a)
        b = gene_to_idx.get(e.gene_b)
        if a is None or b is None:
            continue
        if a == b:
            continue
        if a > b:
            a, b = b, a
        if (a, b) in seen:
            # Keep max score for duplicate pairs
            for idx, (s, t, w) in enumerate(edges):
                if s == a and t == b:
                    if e.score > w:
                        edges[idx] = (a, b, e.score)
                    break
            continue
        seen.add((a, b))
        edges.append((a, b, e.score))

    # Build features
    features, feature_names = _build_features(gene_ids, gene_to_idx, species, edges)

    return CrossSpeciesGraph(
        gene_ids=gene_ids,
        gene_to_idx=gene_to_idx,
        species=species,
        labels={g: int(labels.get(g, 0)) for g in gene_ids},
        edges=edges,
        features=features,
        feature_names=feature_names,
    )


def _build_features(
    gene_ids: list[str],
    gene_to_idx: dict[str, int],
    species: list[str],
    edges: list[tuple[int, int, float]],
) -> tuple[list[list[float]], list[str]]:
    """Build real, simple node features.

    Features:
      1. degree (count, then normalized by max degree)
      2. mean incident orthology confidence
      3. max incident orthology confidence
      4. species one-hot (one column per unique species)

    All features are real and computable from orthology data alone.
    """
    n = len(gene_ids)
    # Degree + confidence accumulators
    degree = [0] * n
    conf_sum = [0.0] * n
    conf_max = [0.0] * n

    for s, t, w in edges:
        degree[s] += 1
        degree[t] += 1
        conf_sum[s] += w
        conf_sum[t] += w
        conf_max[s] = max(conf_max[s], w)
        conf_max[t] = max(conf_max[t], w)

    max_deg = max(degree) if degree and max(degree) > 0 else 1
    unique_species = sorted(set(species))
    species_to_idx = {s: i for i, s in enumerate(unique_species)}

    feature_names = ["degree_norm", "mean_confidence", "max_confidence"] + [
        f"species_{s}" for s in unique_species
    ]

    features: list[list[float]] = []
    for i in range(n):
        deg_norm = degree[i] / max_deg
        mean_conf = (conf_sum[i] / degree[i]) if degree[i] > 0 else 0.0
        max_conf = conf_max[i]
        one_hot = [0.0] * len(unique_species)
        one_hot[species_to_idx[species[i]]] = 1.0
        row = [deg_norm, mean_conf, max_conf] + one_hot
        features.append(row)

    return features, feature_names


def generate_synthetic_graph(
    num_source_genes: int = 20,
    num_target_genes: int = 15,
    num_ortholog_pairs: int = 25,
    go_term: str = "GO:0007049",
    conservation_rate: float = 0.8,
    seed: int = 42,
) -> CrossSpeciesGraph:
    """Generate an explicitly synthetic two-species orthology graph with injected signal.

    This is NOT real biological data — it is a synthetic correctness test fixture.
    Most ortholog pairs share function (conserved), some diverged pairs do not,
    allowing the pipeline's ability to recover transferable signal to be verified.

    Args:
        num_source_genes: number of source-species genes
        num_target_genes: number of target-species genes
        num_ortholog_pairs: number of ortholog edges
        go_term: GO term name (informational)
        conservation_rate: fraction of ortholog edges where function is conserved (shared label)
        seed: random seed

    Returns:
        CrossSpeciesGraph with synthetic labels and edges
    """
    rng = random.Random(seed)
    source_genes = [f"SOURCE_GENE_{i:04d}" for i in range(num_source_genes)]
    target_genes = [f"TARGET_GENE_{i:04d}" for i in range(num_target_genes)]

    # Species taxids (synthetic, using real taxid-like strings)
    source_species = "559292"  # S. cerevisiae
    target_species = "4896"  # S. pombe (example target)

    # Assign source genes: ~40% positive for GO term
    source_labels: dict[str, int] = {}
    for g in source_genes:
        source_labels[g] = 1 if rng.random() < 0.4 else 0

    # Create ortholog pairs
    ortholog_edges: list[OrthologEdge] = []
    target_labels: dict[str, int] = {}

    # First, create pairs ensuring coverage
    pairs_made: set[tuple[str, str]] = set()
    for _ in range(num_ortholog_pairs):
        sg = rng.choice(source_genes)
        tg = rng.choice(target_genes)
        if (sg, tg) in pairs_made:
            continue
        pairs_made.add((sg, tg))

        # Injected signal: with conservation_rate, target inherits source label
        if rng.random() < conservation_rate:
            target_labels[tg] = source_labels[sg]
        else:
            # Diverged: flip
            target_labels[tg] = 1 - source_labels[sg]

        score = rng.uniform(0.6, 0.99)
        ortholog_edges.append(
            OrthologEdge(
                og_id=f"OG_SYNTH_{rng.randint(1, 100):04d}",
                gene_a=sg,
                species_a=source_species,
                gene_b=tg,
                species_b=target_species,
                score=score,
            )
        )

    # Unpaired target genes: random label
    for tg in target_genes:
        if tg not in target_labels:
            target_labels[tg] = 1 if rng.random() < 0.3 else 0

    labels = {**source_labels, **target_labels}

    # Build graph with explicit species mapping so isolated nodes keep correct species
    # (build_graph alone would set "unknown" for isolated label-only genes)
    known_species: dict[str, str] = {}
    for g in source_genes:
        known_species[g] = source_species
    for g in target_genes:
        known_species[g] = target_species

    # Collect genes preserving known species
    gene_species: dict[str, str] = {}
    for e in ortholog_edges:
        if e.gene_a not in gene_species:
            gene_species[e.gene_a] = e.species_a
        if e.gene_b not in gene_species:
            gene_species[e.gene_b] = e.species_b
    for g in labels:
        if g not in gene_species:
            gene_species[g] = known_species.get(g, "unknown")

    gene_ids = sorted(gene_species.keys())
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    species_list = [gene_species[g] for g in gene_ids]

    # Deduplicate edges
    seen: set[tuple[int, int]] = set()
    edges_idx: list[tuple[int, int, float]] = []
    for e in ortholog_edges:
        a = gene_to_idx[e.gene_a]
        b = gene_to_idx[e.gene_b]
        if a == b:
            continue
        if a > b:
            a, b = b, a
        if (a, b) in seen:
            continue
        seen.add((a, b))
        edges_idx.append((a, b, e.score))

    features, feature_names = _build_features(gene_ids, gene_to_idx, species_list, edges_idx)
    return CrossSpeciesGraph(
        gene_ids=gene_ids,
        gene_to_idx=gene_to_idx,
        species=species_list,
        labels={g: int(labels.get(g, 0)) for g in gene_ids},
        edges=edges_idx,
        features=features,
        feature_names=feature_names,
    )
