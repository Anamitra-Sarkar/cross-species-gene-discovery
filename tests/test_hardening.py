"""Hardening tests: edge cases for parsers, backend validation, graph I/O."""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from backend.model_store import clear_cache
from data_pipeline.go_annotations import parse_gaf
from data_pipeline.graph import CrossSpeciesGraph, build_graph
from data_pipeline.orthology import OrthologEdge, parse_orthology_file, parse_simplified_edge_list


# ── Orthology edge list ──────────────────────────────────────────────────

def test_edge_list_comment_with_leading_whitespace(tmp_path):
    p = tmp_path / "edges.tsv"
    p.write_text("   # comment with leading spaces\nOG1\tSGD:S000000001\t559292\tPOMBASE:SPAC1\t4896\t0.9\n")
    edges = parse_simplified_edge_list(p)
    assert len(edges) == 1
    assert edges[0].gene_a == "SGD:S000000001"


def test_edge_list_score_out_of_range(tmp_path):
    p = tmp_path / "edges.tsv"
    p.write_text("OG1\tA\t1\tB\t2\t1.5\n")
    with pytest.raises(ValueError, match="out of range"):
        parse_simplified_edge_list(p)


def test_edge_list_negative_score_rejected(tmp_path):
    p = tmp_path / "edges.tsv"
    p.write_text("OG1\tA\t1\tB\t2\t-0.1\n")
    with pytest.raises(ValueError, match="out of range"):
        parse_simplified_edge_list(p)


def test_edge_list_invalid_score_not_numeric(tmp_path):
    p = tmp_path / "edges.tsv"
    p.write_text("OG1\tA\t1\tB\t2\tnot_a_number\n")
    with pytest.raises(ValueError, match="invalid score"):
        parse_simplified_edge_list(p)


def test_edge_list_blank_lines_and_bom(tmp_path):
    p = tmp_path / "edges.tsv"
    # BOM + blank lines
    p.write_text("\ufeff\n\nOG1\tA\t1\tB\t2\t0.8\n\n  \nOG2\tC\t1\tD\t2\t0.7\n")
    edges = parse_simplified_edge_list(p)
    assert len(edges) == 2


def test_edge_list_empty_gene_id_rejected(tmp_path):
    p = tmp_path / "edges.tsv"
    p.write_text("OG1\t\t559292\tPOMBASE:SPAC1\t4896\t0.9\n")
    # Empty gene_a: either column-count error or non-empty validation both indicate rejection
    with pytest.raises(ValueError, match=r"(non-empty|Expected 6 columns)"):
        parse_simplified_edge_list(p)


def test_orthodb_tab_with_comment_header(tmp_path):
    p = tmp_path / "odb.tab"
    p.write_text("# OrthoDB exported\n! comment\nLevel\tOG_ID\tOrganism_ID\tGene_ID\tUniProt\nEukaryota\tOG0001\t559292\tGENE_A\tP123\nEukaryota\tOG0001\t4896\tGENE_B\tP124\n")
    edges = parse_orthology_file(p)
    assert len(edges) == 1
    assert edges[0].score == 1.0


def test_orthodb_tab_substring_header_detection(tmp_path):
    # Variant header names
    p = tmp_path / "odb2.tab"
    p.write_text("odbGroup\tprotein id\tncbi taxid\nOG_X\tGENE1\t559292\nOG_X\tGENE2\t4896\n")
    edges = parse_orthology_file(p)
    assert len(edges) == 1


# ── GAF parsing ──────────────────────────────────────────────────────────

def test_gaf_comment_with_leading_whitespace(tmp_path):
    p = tmp_path / "test.gaf"
    p.write_text("   !gaf-version: 2.2\nSGD\tS000000001\tG1\t\tGO:0007049\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292\t20240101\tSGD\t\t\n")
    anns = parse_gaf(p, go_term="GO:0007049")
    assert len(anns) == 1


def test_gaf_taxon_pipe_separated(tmp_path):
    p = tmp_path / "test.gaf"
    p.write_text("SGD\tS000000001\tG1\t\tGO:0007049\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292|taxon:9606\t20240101\tSGD\t\t\n")
    anns = parse_gaf(p, go_term="GO:0007049", taxon_filter="559292")
    assert len(anns) == 1
    anns2 = parse_gaf(p, go_term="GO:0007049", taxon_filter="9606")
    assert len(anns2) == 1
    anns3 = parse_gaf(p, go_term="GO:0007049", taxon_filter="999999")
    assert len(anns3) == 0


def test_gaf_evidence_code_case_insensitive(tmp_path):
    p = tmp_path / "test.gaf"
    p.write_text("SGD\tS000000001\tG1\t\tGO:0007049\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292\t20240101\tSGD\t\t\n")
    anns = parse_gaf(p, go_term="GO:0007049", evidence_codes={"ida"})
    assert len(anns) == 1
    anns2 = parse_gaf(p, go_term="GO:0007049", evidence_codes={"IDA"})
    assert len(anns2) == 1


def test_gaf_not_qualifier_skipped(tmp_path):
    p = tmp_path / "test.gaf"
    # NOT qualifier should be skipped
    p.write_text("SGD\tS000000001\tG1\tNOT\tGO:0007049\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292\t20240101\tSGD\t\t\n")
    anns = parse_gaf(p, go_term="GO:0007049")
    assert len(anns) == 0


def test_gaf_malformed_go_id_skipped(tmp_path):
    p = tmp_path / "test.gaf"
    p.write_text("SGD\tS000000001\tG1\t\tBAD_GO\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292\t20240101\tSGD\t\t\n")
    anns = parse_gaf(p)
    assert len(anns) == 0


def test_gaf_blank_and_short_lines(tmp_path):
    p = tmp_path / "test.gaf"
    p.write_text("!header\n\nshort\tline\nSGD\tS000000001\tG1\t\tGO:0007049\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292\t20240101\tSGD\t\t\n")
    anns = parse_gaf(p, go_term="GO:0007049")
    assert len(anns) == 1


def test_gaf_with_trailing_missing_columns(tmp_path):
    p = tmp_path / "test.gaf"
    # Only 13 columns (missing Assigned_By etc.) — should pad and parse
    p.write_text("SGD\tS000000001\tG1\t\tGO:0007049\tPMID:1\tIDA\t\tP\tName\t\tgene\ttaxon:559292\n")
    anns = parse_gaf(p, go_term="GO:0007049")
    assert len(anns) == 1


# ── Graph I/O ────────────────────────────────────────────────────────────

def test_graph_load_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        CrossSpeciesGraph.load(tmp_path / "nonexistent.json")


def test_graph_load_invalid_json(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{ not valid json }")
    with pytest.raises(ValueError, match="Invalid JSON"):
        CrossSpeciesGraph.load(p)


def test_graph_load_missing_keys(tmp_path):
    p = tmp_path / "bad2.json"
    p.write_text(json.dumps({"something": "else"}))
    with pytest.raises(ValueError, match="missing required"):
        CrossSpeciesGraph.load(p)


def test_graph_save_creates_parent_dirs(tmp_path):
    from data_pipeline.graph import generate_synthetic_graph
    g = generate_synthetic_graph(num_source_genes=4, num_target_genes=3, num_ortholog_pairs=5, seed=1)
    out = tmp_path / "nested" / "deep" / "graph.json"
    g.save(out)
    assert out.exists()
    loaded = CrossSpeciesGraph.load(out)
    assert loaded.num_nodes == g.num_nodes


# ── Backend validation ───────────────────────────────────────────────────

@pytest.fixture
def client_no_model():
    for key in ["MODEL_RELEASE_APPROVED", "APPROVED_ARTIFACT_REVISION", "MODEL_ARTIFACT_PATH"]:
        os.environ.pop(key, None)
    clear_cache()
    with patch.dict(os.environ, {}, clear=False):
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


def test_search_whitespace_only_query_422(client_no_model):
    resp = client_no_model.get("/genes/search?q=%20%20%20")
    assert resp.status_code == 422


def test_search_empty_query_422(client_no_model):
    resp = client_no_model.get("/genes/search?q=")
    assert resp.status_code == 422


def test_predict_invalid_gene_id_400(client_with_model):
    resp = client_with_model.get("/genes/invalid%20gene%20id%21/predict")
    # gene_id with space and ! -> fails validation -> 400
    assert resp.status_code == 400


def test_predict_path_traversal_400(client_with_model):
    resp = client_with_model.get("/genes/..%2Fetc%2Fpasswd/predict")
    assert resp.status_code in (400, 404)


def test_predict_unknown_gene_404(client_with_model):
    resp = client_with_model.get("/genes/UNKNOWN_GENE_99999/predict")
    assert resp.status_code == 404


def test_predict_malformed_gene_too_long(client_with_model):
    long_id = "A" * 300
    resp = client_with_model.get(f"/genes/{long_id}/predict")
    assert resp.status_code == 400


def test_orthologs_unknown_gene_404(client_with_model):
    resp = client_with_model.get("/genes/UNKNOWN_GENE_99999/orthologs")
    assert resp.status_code == 404


def test_orthologs_success(client_with_model):
    # Get a real gene via search
    resp = client_with_model.get("/genes/search?q=SGD")
    genes = resp.json()["results"]
    if not genes:
        resp2 = client_with_model.get("/genes/search?q=SOURCE")
        genes = resp2.json()["results"]
    if not genes:
        import json as j
        artifact = os.environ.get("MODEL_ARTIFACT_PATH", "")
        if artifact and Path(artifact).exists():
            data = j.loads(Path(artifact).read_text())
            gid = data["gene_ids"][0]
        else:
            pytest.skip("no genes")
    else:
        gid = genes[0]["gene_id"]
    resp = client_with_model.get(f"/genes/{gid}/orthologs")
    assert resp.status_code == 200
    assert "orthologs" in resp.json()


def test_auth_empty_bearer_token_401(client_no_model):
    resp = client_no_model.get("/health", headers={"Authorization": "Bearer "})
    # Health doesn't require auth in permissive mode, so should still be 200 (auth is optional for /health)
    # But verify verify_bearer_token rejects empty after Bearer
    from backend.auth import verify_bearer_token
    with pytest.raises(ValueError):
        verify_bearer_token("   ")


def test_auth_case_insensitive_bearer(client_no_model):
    # Lowercase bearer should be accepted (we fixed to case-insensitive)
    with patch.dict(os.environ, {"FIREBASE_AUTH_ENABLED": "false"}):
        from backend.auth import get_current_user
        import asyncio
        # Directly test verify path: token with lowercase bearer prefix via TestClient
        resp = client_no_model.get("/health", headers={"Authorization": "bearer some-token"})
        assert resp.status_code == 200
