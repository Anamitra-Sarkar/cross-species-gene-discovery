# Architecture

## Overview

Evolution-aware graph learning for cross-species gene function discovery. The system transfers Gene Ontology (GO) functional annotations from well-annotated source species (e.g. *S. cerevisiae* yeast, *M. musculus* mouse) to a less-annotated target species via a cross-species orthology graph.

```
Real data sources (no auth)          Pipeline                    Models                    Serving
┌─────────────────────┐     ┌─────────────────────┐     ┌──────────────────┐     ┌──────────────┐
│ OrthoDB / eggNOG     │────▶│ data_pipeline/      │────▶│ model/           │────▶│ backend/     │
│ orthology groups     │     │ - orthology.py      │     │ - rwr.py (RWR)   │     │ FastAPI +    │
│ GO / GOA (GAF)       │     │ - go_annotations.py │     │ - gnn.py (GNN)   │     │ release gate │
│ annotations          │     │ - graph.py          │     │ - evaluation.py  │     │ + auth stub  │
│                     │     │ - build_graph.py    │     │                  │     │ frontend/    │
└─────────────────────┘     └─────────────────────┘     └──────────────────┘     │ React/Vite   │
                                                                                 └──────────────┘
```

## Components

### 1. `data_pipeline/` — Data ingestion and graph construction

- **`orthology.py`** — Parses OrthoDB/eggNOG orthology data. Supports two input formats:
  - *Simplified edge list* (TSV: `og_id  gene_id_1  species_1  gene_id_2  species_2  score`) — used by fixtures and bulk post-processing.
  - *Native OrthoDB TAB* (`*_OGs.tab` with columns `og_id, gene_id, taxid, ...`) — expands ortholog groups into pairwise edges.
  CLI: `python -m data_pipeline.orthology --orthology-path <tsv> --output <json>`

- **`go_annotations.py`** — Parses GO Annotation File (GAF 2.2) format (17-column TSV, see `docs/data_sources.md`). Filters by GO term, evidence codes, taxon, aspect. Produces binary label vectors per gene.
  CLI flags: `--go-annotations-path`, `--go-term`, `--evidence-codes`, `--taxon-filter`

- **`graph.py`** — Builds a cross-species orthology graph as a `torch_geometric.data.Data` object (or pure Python/NetworkX fallback for pipeline without torch). Nodes = genes (with species attribute), edges = ortholog pairs (weighted by confidence score). Node features: degree, orthology-confidence-derived scalar(s), species one-hot — real, simple, computable without sequence embeddings.
  CLI: `python -m data_pipeline.build_graph --orthology-path ... --go-annotations-path ... --output data/graph.pt`

- **`build_graph.py`** — End-to-end CLI that wires orthology + GO parsing + graph construction. No downloads in sandbox; real-run procedure documented in `docs/data_sources.md`.

### 2. `model/` — Graph learning

- **`rwr.py`** — Random Walk with Restart (RWR) / label propagation baseline. A real, standard, simple baseline for network-based function prediction.
  Formula: `p_{t+1} = (1 - alpha) * W * p_t + alpha * p_0`
  where `W` is the column-normalized adjacency matrix, `p_0` is the restart vector (1 on positively-labeled source genes, 0 elsewhere, normalized), `alpha` is restart probability. Iterated to convergence (L1 delta < tol). Implemented with NumPy/SciPy only (no torch required for baseline).

- **`gnn.py`** — GraphSAGE / GCN using `torch_geometric`. Real GNN that learns from graph structure + node features. Architecture: 2-layer GraphSAGE (or GCN) with ReLU + dropout, final linear for binary GO-term classification. Trained with BCE loss on source-species nodes only. Node features are real orthology-derived scalars (see `graph.py`).

- **`evaluation.py`** — **Cross-species holdout** evaluation. The honest protocol:
  1. Split target-species genes into train (hidden) / test (held-out). Training uses *only* source-species labels + graph structure.
  2. Evaluate predictions on held-out target-species genes against their real GO annotations (ground truth exists but is hidden during training).
  3. Report AUROC, AUPRC vs naive majority-class baseline.
  4. If no real dual-annotated pair is wired, the pipeline supports a **synthetic two-species graph with injected conservation signal** (most orthologs share function, some diverged) — explicitly documented as synthetic verification, not a biological finding.

### 3. `backend/` — FastAPI serving

- **`app.py`** — FastAPI application.
- **Release gate (fail-closed):** Model artifacts are NOT loaded/served unless `MODEL_RELEASE_APPROVED=true` AND `APPROVED_ARTIFACT_REVISION` is set and points to a valid artifact. Health/readiness endpoints (`/health`, `/ready`) honestly reflect whether a real approved model is loaded. Predictions return 503 with `model not yet released` when gate is closed.
- **Auth stub:** Firebase-auth-shaped bearer-token verification. Reads a JSON service-account path from `FIREBASE_SERVICE_ACCOUNT_JSON` env var. If no file is present (sandbox), it runs in permissive mock mode (documented). Includes a real `verify_bearer_token(token)` function that validates JWT structure and is unit-tested with mocked verifier.
- **Endpoints:**
  - `GET /health` — liveness, always 200.
  - `GET /ready` — readiness, 200 only if model loaded via release gate.
  - `GET /genes/search?q=<query>` — search target-species genes.
  - `GET /genes/{gene_id}/predict` — predicted GO-term score + ortholog explanation (which source orthologs drove prediction, with confidence).

### 4. `frontend/` — React + Vite + TypeScript

- Tooling: Vite + TypeScript + React.
- **Real, non-boilerplate UI:**
  - Target-species gene search box (calls `GET /genes/search`).
  - Gene detail card: predicted GO term / function score, explanation list (source orthologs with scores), species badges.
  - Honest abstention banner: shows "Model not yet released — predictions unavailable" when `GET /ready` reports not-ready, matching backend release gate.
  - Clean, modern scientific-dashboard design (CSS, no heavy UI framework dependency).
- Dev: `npm install && npm run dev` (Vite). Build: `npm run build`.

### 5. `tests/` — Pytest

- Small synthetic fixture graphs (10-50 genes, 2-3 synthetic species) exercising:
  - Orthology graph construction correctness (node/edge counts, weights, undirected symmetry)
  - RWR / label propagation correctness (convergence, score ordering: closer orthologs score higher)
  - GNN forward pass shape correctness (output shape = [num_nodes, 1] or [num_nodes] for binary classification)
  - Cross-species holdout split correctness (no leakage: target test genes never seen during training)
  - Backend API tests (release-gate on/off, auth stub mocking)
- Fixtures in `tests/fixtures/` match real file formats (OrthoDB-like TSV, GAF 2.2).

## Data Flow (real run)

```bash
# 1. Download (real, no auth) — documented, not executed in sandbox
curl -L https://www.orthodb.org/download/odb11v0_Eukaryota_OGs.tab.gz -o data/odb.tab.gz
curl -L https://current.geneontology.org/annotations/sgd.gaf.gz -o data/sgd.gaf.gz

# 2. Build graph
python -m data_pipeline.build_graph \
  --orthology-path data/odb.tab \
  --go-annotations-path data/sgd.gaf \
  --go-term GO:0007049 \
  --output data/graph.pt

# 3. Train & evaluate
python -m model.evaluation --graph data/graph.pt --target-taxon 4896

# 4. Serve (requires explicit release approval)
MODEL_RELEASE_APPROVED=true APPROVED_ARTIFACT_REVISION=v1.0 uvicorn backend.app:app --reload
```

## Design Decisions

- **OrthoDB over eggNOG as primary:** OrthoDB has a real REST API with per-gene lookup (useful for interactive serving) and clear bulk TAB format. Both are cited; eggNOG remains a documented alternative.
- **RWR + GNN (both):** RWR is a well-established, interpretable baseline (Köhler et al. 2008; Cowen et al. 2017 for network propagation in biology). GNN (GraphSAGE/GCN) is the learned counterpart. Reporting both is standard.
- **Node features are honest/real:** Orthology-confidence-derived scalars + degree + species one-hot. No fabricated sequence embeddings. If sequence embeddings (e.g. ESM) are available later, they plug in as additional `x` columns.
- **Fail-closed release gate:** Prevents serving stale/unapproved models. Pattern: env-var-gated artifact loading, `/ready` reflects truth, predictions 503 when not ready.
- **Synthetic fixtures are explicit:** Tests use synthetic graphs with injected signal, clearly labeled as verification, not biological claims.
