"""Model artifact store with fail-closed release gate.

Model artifacts are NOT loaded/served unless explicit env vars are set:
  MODEL_RELEASE_APPROVED=true
  APPROVED_ARTIFACT_REVISION=<revision string>

This is the developer's standard fail-closed release gate pattern.
Health/readiness endpoints honestly reflect whether a real approved model is loaded,
never fabricate a loaded state.

Artifact lookup:
  - If MODEL_ARTIFACT_PATH is set, load graph/model from that path
  - Otherwise, no model is loaded (even if gate is open but no artifact exists)

The loaded model is cached in memory after successful gate check.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class ModelStatus:
    approved: bool
    revision: Optional[str]
    loaded: bool
    artifact_path: Optional[str]
    message: str


# Module-level cache for loaded graph / model
_cached_graph = None
_cached_status: Optional[ModelStatus] = None


def check_release_gate() -> ModelStatus:
    """Check the fail-closed release gate and return status.

    Gate is OPEN only if:
      MODEL_RELEASE_APPROVED == "true" (case-insensitive)
      AND APPROVED_ARTIFACT_REVISION is non-empty

    Even when gate is open, loaded=False if no artifact file exists at MODEL_ARTIFACT_PATH.
    """
    approved_flag = os.environ.get("MODEL_RELEASE_APPROVED", "").lower() == "true"
    revision = os.environ.get("APPROVED_ARTIFACT_REVISION", "").strip() or None
    artifact_path = os.environ.get("MODEL_ARTIFACT_PATH", "").strip() or None

    if not approved_flag:
        return ModelStatus(
            approved=False,
            revision=revision,
            loaded=False,
            artifact_path=artifact_path,
            message="Model not yet released — MODEL_RELEASE_APPROVED != true",
        )

    if not revision:
        return ModelStatus(
            approved=False,
            revision=None,
            loaded=False,
            artifact_path=artifact_path,
            message="Model not yet released — APPROVED_ARTIFACT_REVISION not set",
        )

    # Gate open, but check artifact exists
    if not artifact_path:
        return ModelStatus(
            approved=True,
            revision=revision,
            loaded=False,
            artifact_path=None,
            message=f"Release approved (rev={revision}) but MODEL_ARTIFACT_PATH not set — no artifact to load",
        )

    if not Path(artifact_path).exists():
        return ModelStatus(
            approved=True,
            revision=revision,
            loaded=False,
            artifact_path=artifact_path,
            message=f"Release approved (rev={revision}) but artifact not found at {artifact_path}",
        )

    return ModelStatus(
        approved=True,
        revision=revision,
        loaded=True,
        artifact_path=artifact_path,
        message=f"Model loaded (rev={revision}) from {artifact_path}",
    )


def get_model_status() -> ModelStatus:
    """Get current model status (checks gate each call, no stale cache for status)."""
    return check_release_gate()


def is_model_ready() -> bool:
    """Whether a real approved model is loaded and ready to serve."""
    return get_model_status().loaded


def load_approved_graph():
    """Load the approved graph artifact if gate is open and artifact exists.

    Returns:
        CrossSpeciesGraph if loaded, None otherwise.
    """
    status = get_model_status()
    if not status.loaded or not status.artifact_path:
        return None

    global _cached_graph
    # Simple cache keyed by artifact path + revision
    # Invalidate if path changes
    cache_key = (status.artifact_path, status.revision)
    if _cached_graph is not None and getattr(_cached_graph, "_cache_key", None) == cache_key:
        return _cached_graph

    # Load graph
    try:
        from data_pipeline.graph import CrossSpeciesGraph

        graph = CrossSpeciesGraph.load(status.artifact_path)
        graph._cache_key = cache_key  # type: ignore
        _cached_graph = graph
        return graph
    except Exception:
        return None


def clear_cache() -> None:
    """Clear cached graph (for testing)."""
    global _cached_graph
    _cached_graph = None
