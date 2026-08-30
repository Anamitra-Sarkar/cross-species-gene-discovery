"""GNN for cross-species gene function prediction.

Real GraphSAGE / GCN using torch_geometric.
Architecture: 2-layer GraphSAGE (or GCN) + ReLU + dropout + final linear for binary GO-term classification.
Trained with BCEWithLogitsLoss on source-species nodes only (cross-species holdout).

Node features are real orthology-derived scalars (degree, confidence, species one-hot) —
not fabricated sequence embeddings. If ESM/protein embeddings are available later,
they concatenate as additional x columns.

References:
  - Hamilton et al. "Inductive Representation Learning on Large Graphs (GraphSAGE)." NeurIPS 2017.
  - Kipf & Welling "Semi-Supervised Classification with Graph Convolutional Networks." ICLR 2017.
  - torch_geometric documentation: https://pytorch-geometric.readthedocs.io/
"""

from __future__ import annotations

from typing import Optional

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import Data
    from torch_geometric.nn import SAGEConv, GCNConv

    HAS_TORCH_GEOMETRIC = True
except ImportError:
    HAS_TORCH_GEOMETRIC = False
    torch = None  # type: ignore
    nn = None  # type: ignore
    F = None  # type: ignore
    Data = object  # type: ignore
    SAGEConv = None  # type: ignore
    GCNConv = None  # type: ignore


if HAS_TORCH_GEOMETRIC:

    class CrossSpeciesGNN(nn.Module):  # type: ignore
        """2-layer GraphSAGE/GCN for binary GO-term prediction.

        Args:
            in_channels: number of input node features
            hidden_channels: hidden dimension
            out_channels: 1 for binary classification
            num_layers: 2 (default)
            conv_type: "sage" or "gcn"
            dropout: dropout rate
        """

        def __init__(
            self,
            in_channels: int,
            hidden_channels: int = 32,
            out_channels: int = 1,
            num_layers: int = 2,
            conv_type: str = "sage",
            dropout: float = 0.3,
        ):
            super().__init__()
            self.dropout = dropout
            self.num_layers = num_layers
            ConvClass = SAGEConv if conv_type == "sage" else GCNConv

            self.convs = nn.ModuleList()
            self.convs.append(ConvClass(in_channels, hidden_channels))
            for _ in range(num_layers - 2):
                self.convs.append(ConvClass(hidden_channels, hidden_channels))
            # For 2 layers, we have 1 hidden + 1 output; for >2 layers, hidden layers in between
            if num_layers >= 2:
                self.convs.append(ConvClass(hidden_channels, hidden_channels))

            self.lin = nn.Linear(hidden_channels, out_channels)

        def forward(self, x: torch.Tensor, edge_index: torch.Tensor, edge_weight: Optional[torch.Tensor] = None) -> torch.Tensor:  # type: ignore
            for conv in self.convs:
                # SAGEConv and GCNConv handle edge_index; edge_weight only for GCN
                if isinstance(conv, SAGEConv):
                    x = conv(x, edge_index)
                else:
                    x = conv(x, edge_index, edge_weight=edge_weight)
                x = F.relu(x)
                x = F.dropout(x, p=self.dropout, training=self.training)
            x = self.lin(x)
            return x.squeeze(-1)  # [num_nodes]

        def predict_proba(self, x: torch.Tensor, edge_index: torch.Tensor, edge_weight: Optional[torch.Tensor] = None) -> torch.Tensor:  # type: ignore
            self.eval()
            with torch.no_grad():
                logits = self.forward(x, edge_index, edge_weight)
                return torch.sigmoid(logits)

else:

    class CrossSpeciesGNN:  # type: ignore
        def __init__(self, *args, **kwargs):
            raise ImportError("torch and torch_geometric are required for CrossSpeciesGNN. Install with: pip install torch torch_geometric")

        def forward(self, *args, **kwargs):
            raise ImportError("torch and torch_geometric not installed")


def build_gnn(
    in_channels: int,
    hidden_channels: int = 32,
    conv_type: str = "sage",
    dropout: float = 0.3,
) -> "CrossSpeciesGNN":
    """Factory for CrossSpeciesGNN."""
    return CrossSpeciesGNN(
        in_channels=in_channels,
        hidden_channels=hidden_channels,
        conv_type=conv_type,
        dropout=dropout,
    )


def train_gnn(
    data: "Data",
    train_mask: "torch.Tensor",
    hidden_channels: int = 32,
    conv_type: str = "sage",
    lr: float = 0.01,
    weight_decay: float = 5e-4,
    epochs: int = 100,
    dropout: float = 0.3,
    verbose: bool = False,
) -> "CrossSpeciesGNN":
    """Train GNN on source-species nodes (train_mask).

    Args:
        data: torch_geometric Data with x, edge_index, edge_weight, y
        train_mask: boolean tensor [num_nodes], True for source training nodes
        hidden_channels, conv_type, lr, weight_decay, epochs, dropout: hyperparams

    Returns:
        trained CrossSpeciesGNN model
    """
    if not HAS_TORCH_GEOMETRIC:
        raise ImportError("torch and torch_geometric required for training")

    in_channels = data.x.size(1)
    model = CrossSpeciesGNN(in_channels=in_channels, hidden_channels=hidden_channels, conv_type=conv_type, dropout=dropout)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = nn.BCEWithLogitsLoss()

    model.train()
    for epoch in range(epochs):
        optimizer.zero_grad()
        logits = model(data.x, data.edge_index, getattr(data, "edge_attr", None))
        loss = criterion(logits[train_mask], data.y[train_mask].float())
        loss.backward()
        optimizer.step()
        if verbose and (epoch + 1) % 20 == 0:
            print(f"Epoch {epoch+1}/{epochs} loss={loss.item():.4f}")

    return model
