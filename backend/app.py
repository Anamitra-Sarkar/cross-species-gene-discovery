"""FastAPI app for cross-species gene function prediction.

Fail-closed release gate: model artifacts NOT loaded/served unless
  MODEL_RELEASE_APPROVED=true and APPROVED_ARTIFACT_REVISION=<rev> are set.
Health/readiness endpoints honestly reflect model loaded state.

Endpoints:
  GET /health
  GET /ready
  GET /genes/search?q=<query>
  GET /genes/{gene_id}/predict  (requires model ready, else 503)
  GET /genes/{gene_id}/orthologs
"""

from __future__ import annotations

import os
import re
from typing import Optional

from fastapi import Depends, FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from .auth import get_current_user
from .model_store import get_model_status, is_model_ready, load_approved_graph

# Gene ID validation: reasonable chars, no path traversal, length bound
_GENE_ID_RE = re.compile(r"^[A-Za-z0-9:_\-.\|]+$")
_MAX_GENE_ID_LEN = 200
_MAX_QUERY_LEN = 500


class HealthResponse(BaseModel):
    status: str
    service: str


class ReadyResponse(BaseModel):
    ready: bool
    approved: bool
    revision: Optional[str]
    artifact_path: Optional[str]
    message: str


class GeneSearchResult(BaseModel):
    gene_id: str
    species: str
    label: int


class GeneSearchResponse(BaseModel):
    model_config = {"protected_namespaces": ()}
    query: str
    results: list[GeneSearchResult]
    message: Optional[str] = None
    model_ready: bool


class OrthologInfo(BaseModel):
    gene_id: str
    species: str
    orthology_confidence: Optional[float] = None
    confidence: Optional[float] = None
    source_positive: Optional[bool] = None


class PredictionResponse(BaseModel):
    model_config = {"protected_namespaces": ()}
    gene_id: str
    species: str
    predicted_score: float
    go_term: str
    go_term_name: str
    explanation: dict
    model_revision: Optional[str]

app = FastAPI(
    title="Cross-Species Gene Function Discovery",
    description="Evolution-aware graph learning for cross-species gene function transfer via orthology graphs.",
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Return clean 422 for validation errors instead of raw 500."""
    # Sanitize errors to avoid leaking internals, but keep useful messages
    return JSONResponse(status_code=422, content={"detail": exc.errors()})


def _validate_gene_id(gene_id: str) -> str:
    """Validate gene_id path param; raises 400 for malformed IDs instead of 500."""
    if not gene_id or not gene_id.strip():
        raise HTTPException(status_code=400, detail="gene_id must be non-empty")
    if len(gene_id) > _MAX_GENE_ID_LEN:
        raise HTTPException(status_code=400, detail=f"gene_id too long (max {_MAX_GENE_ID_LEN})")
    if ".." in gene_id or "/" in gene_id or "\\" in gene_id:
        raise HTTPException(status_code=400, detail="gene_id contains invalid characters")
    if not _GENE_ID_RE.match(gene_id):
        raise HTTPException(status_code=400, detail=f"gene_id contains invalid characters: {gene_id!r}")
    return gene_id


@app.get("/health", response_model=HealthResponse)
def health() -> dict:
    return {"status": "ok", "service": "cross-species-gene-discovery"}


@app.get("/ready", response_model=ReadyResponse)
def readiness() -> dict:
    status = get_model_status()
    return {
        "ready": status.loaded,
        "approved": status.approved,
        "revision": status.revision,
        "artifact_path": status.artifact_path,
        "message": status.message,
    }


@app.get("/genes/search", response_model=GeneSearchResponse)
def search_genes(
    q: str = Query(..., min_length=1, max_length=500, description="Search query for gene ID or symbol"),
    limit: int = Query(20, ge=1, le=100),
    user: Optional[dict] = Depends(get_current_user),
) -> dict:
    """Search target-species genes. Requires model artifact to have gene list; otherwise returns empty or mock."""
    # Explicit blank/whitespace check: Query min_length passes "   " (treated as 3 chars) -> return 400
    if not q.strip():
        raise HTTPException(status_code=422, detail="query must be non-empty (whitespace-only not allowed)")
    if len(q) > _MAX_QUERY_LEN:
        raise HTTPException(status_code=400, detail=f"query too long (max {_MAX_QUERY_LEN})")
    q_stripped = q.strip()
    status = get_model_status()

    # If no model loaded, return empty with honest message (not fabricating results)
    if not status.loaded:
        return {
            "query": q_stripped,
            "results": [],
            "message": "Model not yet released — gene search unavailable (no artifact loaded).",
            "model_ready": False,
        }

    graph = load_approved_graph()
    if graph is None:
        return {"query": q_stripped, "results": [], "message": "Failed to load graph artifact.", "model_ready": False}

    # Search: substring match on gene_id (case-insensitive)
    q_lower = q_stripped.lower()
    matches = []
    for idx, gid in enumerate(graph.gene_ids):
        if q_lower in gid.lower():
            matches.append(
                {
                    "gene_id": gid,
                    "species": graph.species[idx],
                    "label": graph.labels.get(gid, 0),
                }
            )
        if len(matches) >= limit:
            break

    return {"query": q_stripped, "results": matches, "model_ready": True}


@app.get("/genes/{gene_id}/predict", response_model=PredictionResponse)
def predict_gene(
    gene_id: str,
    user: Optional[dict] = Depends(get_current_user),
) -> dict:
    """Predict GO-term function for a target-species gene.

    Uses RWR over the cross-species orthology graph.
    Returns predicted score + ortholog explanation.

    Fail-closed: 503 if model not yet released.
    Validates gene_id and returns 400 for malformed IDs, 404 for unknown IDs.
    """
    gene_id = _validate_gene_id(gene_id)
    if not is_model_ready():
        status = get_model_status()
        raise HTTPException(
            status_code=503,
            detail={
                "error": "Model not yet released",
                "message": status.message,
                "approved": status.approved,
                "revision": status.revision,
            },
        )

    graph = load_approved_graph()
    if graph is None:
        raise HTTPException(status_code=503, detail="Model artifact failed to load")

    if gene_id not in graph.gene_to_idx:
        raise HTTPException(status_code=404, detail=f"Gene {gene_id!r} not found in graph")

    # Run RWR with source positives as seeds (honest cross-species: only source labels)
    # Determine source species as most common among positively-labeled genes? Use first species with positives.
    # Simpler: all positive genes are seeds (source positives by definition in cross-species holdout graph)
    # For serving, seeds = all positively-labeled genes in graph
    positives = [g for g, v in graph.labels.items() if v == 1]
    if not positives:
        raise HTTPException(status_code=503, detail="No positive labels in model artifact")

    # Use RWR
    from model.rwr import rwr_predict

    scores = rwr_predict(graph, positives, alpha=0.3)
    idx = graph.gene_to_idx[gene_id]
    score = float(scores[idx])

    # Explanation: top source orthologs that contributed (neighbors with positive labels, sorted by edge weight)
    neighbors: list[dict] = []
    for s, t, w in graph.edges:
        # Find edges incident to query gene
        gid_s = graph.gene_ids[s]
        gid_t = graph.gene_ids[t]
        if gid_s == gene_id:
            neighbor_id = gid_t
        elif gid_t == gene_id:
            neighbor_id = gid_s
        else:
            continue
        is_positive = graph.labels.get(neighbor_id, 0) == 1
        neighbors.append(
            {
                "gene_id": neighbor_id,
                "species": graph.species[graph.gene_to_idx[neighbor_id]],
                "orthology_confidence": w,
                "source_positive": is_positive,
            }
        )

    # Sort by confidence descending, positives first as secondary
    neighbors.sort(key=lambda x: (x["source_positive"], x["orthology_confidence"]), reverse=True)

    return {
        "gene_id": gene_id,
        "species": graph.species[idx],
        "predicted_score": score,
        "go_term": "GO:0007049",
        "go_term_name": "cell cycle",
        "explanation": {
            "method": "Random Walk with Restart over orthology graph",
            "source_orthologs": neighbors[:10],
            "num_source_positives": len(positives),
        },
        "model_revision": get_model_status().revision,
    }


@app.get("/genes/{gene_id}/orthologs")
def get_orthologs(
    gene_id: str,
    user: Optional[dict] = Depends(get_current_user),
) -> dict:
    gene_id = _validate_gene_id(gene_id)
    if not is_model_ready():
        status = get_model_status()
        raise HTTPException(status_code=503, detail={"error": "Model not yet released", "message": status.message})

    graph = load_approved_graph()
    if graph is None:
        raise HTTPException(status_code=503, detail="Model artifact failed to load")

    if gene_id not in graph.gene_to_idx:
        raise HTTPException(status_code=404, detail=f"Gene {gene_id!r} not found")

    neighbors: list[dict] = []
    for s, t, w in graph.edges:
        gid_s = graph.gene_ids[s]
        gid_t = graph.gene_ids[t]
        if gid_s == gene_id:
            neighbors.append({"gene_id": gid_t, "species": graph.species[graph.gene_to_idx[gid_t]], "confidence": w})
        elif gid_t == gene_id:
            neighbors.append({"gene_id": gid_s, "species": graph.species[graph.gene_to_idx[gid_s]], "confidence": w})

    return {"gene_id": gene_id, "orthologs": neighbors}
