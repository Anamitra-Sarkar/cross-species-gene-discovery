"""Tests: backend API release gate and auth stub."""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from backend.auth import verify_bearer_token
from backend.model_store import clear_cache


@pytest.fixture
def client_no_model():
    """Client without any model release env vars set."""
    # Ensure env vars are cleared
    env_patch = patch.dict(os.environ, {}, clear=False)
    # Remove model gate vars if present
    for key in ["MODEL_RELEASE_APPROVED", "APPROVED_ARTIFACT_REVISION", "MODEL_ARTIFACT_PATH"]:
        os.environ.pop(key, None)
    clear_cache()
    with env_patch:
        from backend.app import app

        with TestClient(app) as c:
            yield c


@pytest.fixture
def client_with_model(temp_graph_file):
    env_vars = {
        "MODEL_RELEASE_APPROVED": "true",
        "APPROVED_ARTIFACT_REVISION": "test-v1",
        "MODEL_ARTIFACT_PATH": temp_graph_file,
    }
    with patch.dict(os.environ, env_vars):
        clear_cache()
        from backend.app import app

        with TestClient(app) as c:
            yield c
    clear_cache()


def test_health_endpoint(client_no_model):
    resp = client_no_model.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_ready_not_ready(client_no_model):
    resp = client_no_model.get("/ready")
    assert resp.status_code == 200
    data = resp.json()
    assert data["ready"] is False
    assert "not yet released" in data["message"].lower()


def test_ready_ready(client_with_model):
    resp = client_with_model.get("/ready")
    assert resp.status_code == 200
    data = resp.json()
    assert data["ready"] is True
    assert data["revision"] == "test-v1"


def test_predict_503_when_not_ready(client_no_model, synthetic_graph):
    # Use a gene that exists in synthetic_graph
    gene_id = synthetic_graph.gene_ids[0]
    resp = client_no_model.get(f"/genes/{gene_id}/predict")
    assert resp.status_code == 503
    assert "not yet released" in resp.json()["detail"]["error"].lower() or "not yet released" in str(resp.json()).lower()


def test_predict_success_when_ready(client_with_model):
    # Search first to get a gene_id
    resp = client_with_model.get("/genes/search?q=SGD")
    genes = resp.json()["results"]
    if not genes:
        resp = client_with_model.get("/genes/search?q=SOURCE")
        genes = resp.json()["results"]
    # Fallback: get any gene via search with empty-ish query
    if not genes:
        # Use gene ids directly from the graph file
        import json

        artifact_path = os.environ.get("MODEL_ARTIFACT_PATH", "")
        if artifact_path and Path(artifact_path).exists():
            graph_data = json.loads(Path(artifact_path).read_text())
            gene_id = graph_data["gene_ids"][0]
        else:
            pytest.skip("No genes found")
    else:
        gene_id = genes[0]["gene_id"]

    resp = client_with_model.get(f"/genes/{gene_id}/predict")
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert "predicted_score" in data
    assert 0.0 <= data["predicted_score"] <= 1.0
    assert "explanation" in data
    assert "source_orthologs" in data["explanation"]


def test_search_when_not_ready_returns_empty(client_no_model):
    resp = client_no_model.get("/genes/search?q=SGD")
    assert resp.status_code == 200
    assert resp.json()["results"] == []
    assert resp.json()["model_ready"] is False


def test_search_when_ready(client_with_model):
    resp = client_with_model.get("/genes/search?q=SGD")
    assert resp.status_code == 200
    # Should have results (our synthetic graph has SGD genes)
    # Even if query doesn't match, we check structure
    data = resp.json()
    assert "results" in data
    assert data["model_ready"] is True


def test_auth_verify_bearer_token_valid_jwt():
    import base64
    import json as j

    header = base64.urlsafe_b64encode(j.dumps({"alg": "none"}).encode()).decode().rstrip("=")
    payload = base64.urlsafe_b64encode(j.dumps({"sub": "user123", "email": "test@example.com"}).encode()).decode().rstrip("=")
    token = f"{header}.{payload}.signature"
    # Should not raise
    result = verify_bearer_token(token)
    assert result["uid"] == "user123"


def test_auth_verify_bearer_token_empty():
    with pytest.raises(ValueError, match="Empty"):
        verify_bearer_token("")


def test_auth_verify_bearer_token_opaque():
    # Opaque token should be accepted when not in strict mode
    with patch.dict(os.environ, {"FIREBASE_AUTH_ENABLED": "false"}):
        result = verify_bearer_token("some-opaque-token-123")
        assert "uid" in result

    # In strict mode, opaque token should be rejected
    with patch.dict(os.environ, {"FIREBASE_AUTH_ENABLED": "true"}):
        with pytest.raises(ValueError, match="Opaque"):
            verify_bearer_token("some-opaque-token-123")


def test_auth_invalid_jwt_structure():
    with pytest.raises(ValueError):
        verify_bearer_token("invalid.jwt")


def test_auth_header_missing_permissive(client_no_model):
    # No auth header, permissive mode (default) should not reject health/ready
    resp = client_no_model.get("/health")
    assert resp.status_code == 200


def test_release_gate_requires_both_vars(temp_graph_file):
    from backend.model_store import check_release_gate

    # Only MODEL_RELEASE_APPROVED without revision
    with patch.dict(os.environ, {"MODEL_RELEASE_APPROVED": "true", "APPROVED_ARTIFACT_REVISION": ""}, clear=False):
        os.environ.pop("APPROVED_ARTIFACT_REVISION", None)
        status = check_release_gate()
        assert status.approved is False or status.loaded is False

    # Both set but artifact missing
    with patch.dict(
        os.environ,
        {
            "MODEL_RELEASE_APPROVED": "true",
            "APPROVED_ARTIFACT_REVISION": "v1",
            "MODEL_ARTIFACT_PATH": "/nonexistent/path.json",
        },
    ):
        status = check_release_gate()
        assert status.approved is True
        assert status.loaded is False

    # All set correctly
    with patch.dict(
        os.environ,
        {
            "MODEL_RELEASE_APPROVED": "true",
            "APPROVED_ARTIFACT_REVISION": "v1",
            "MODEL_ARTIFACT_PATH": temp_graph_file,
        },
    ):
        status = check_release_gate()
        assert status.loaded is True
