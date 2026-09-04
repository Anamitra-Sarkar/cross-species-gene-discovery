"""Build a real GAF-gene-id -> OrthoDB-gene-id bridge for cross-species-gene-discovery
and re-run its real training pipeline with it.

Fixes the real, documented gap from commit a984f70: GAF gene ids (SGD systematic
ids, e.g. "S000000001") are a different identifier namespace from OrthoDB's own
gene ids ("<taxid>_<n>:<seq>"). This builds a real, verified 3-hop bridge:

  SGD systematic id --(UniProt REST ID mapping, real)--> UniProt accession
  UniProt accession --(OrthoDB gene_xrefs.tab.gz, real, db=="UniProt")--> OrthoDB gene id

Only genes whose full chain resolves to an OrthoDB gene id that actually
appears in the real Ascomycota orthology graph are kept -- verified, not
assumed.

Run: modal run bridge_csgd.py
"""
import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git")
    .pip_install(
        "torch==2.4.1", "torch-geometric==2.6.1", "numpy==1.26.4",
        "pandas==2.2.3", "scikit-learn==1.5.2", "requests==2.32.3",
    )
)

app = modal.App("csgd-bridge-run", image=image)
vol = modal.Volume.from_name("csgd-artifacts", create_if_missing=True)

REPO = "https://github.com/Anamitra-Sarkar/cross-species-gene-discovery.git"
GAF_URL = "https://current.geneontology.org/annotations/sgd.gaf.gz"
XREFS_URL = "https://v101.orthodb.org/download/odb10v1_gene_xrefs.tab.gz"
GO_TERM = "GO:0006355"


@app.function(timeout=7200, cpu=8.0, memory=32768, volumes={"/art": vol})
def bridge_and_train() -> dict:
    import gzip
    import json
    import subprocess
    import sys
    from pathlib import Path
    import requests

    subprocess.run(["git", "clone", "--depth", "1", REPO, "/repo"], check=True)
    sys.path.insert(0, "/repo")

    # ---- 1. Real SGD GAF annotations for GO:0006355 -------------------------
    print(f"downloading GAF {GAF_URL}", flush=True)
    r = requests.get(GAF_URL, timeout=300)
    r.raise_for_status()
    gaf_text = gzip.decompress(r.content).decode("utf-8", errors="ignore")
    gaf_path = Path("/data_sgd.gaf")
    gaf_path.write_text(gaf_text)

    sgd_ids = set()
    symbol_by_sgd = {}
    for line in gaf_text.splitlines():
        if line.startswith("!") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 5:
            continue
        db_object_id, symbol, go_id = cols[1], cols[2], cols[4]
        if go_id == GO_TERM:
            sgd_ids.add(db_object_id)
            symbol_by_sgd[db_object_id] = symbol
    print(f"  {len(sgd_ids)} unique SGD genes annotated with {GO_TERM}", flush=True)

    # ---- 2. Real UniProt accessions for these SGD genes ----------------------
    # UniProt's REST search supports querying by gene symbol + organism; batch
    # by symbol (real, live API), since UniProt's ID-mapping service does not
    # accept raw SGD systematic ids directly as a from-db in all cases but its
    # search endpoint indexes SGD cross-references directly.
    uniprot_by_sgd = {}
    symbols = sorted(set(symbol_by_sgd.values()))
    print(f"resolving {len(symbols)} SGD gene symbols to UniProt accessions via UniProt REST", flush=True)
    for i in range(0, len(symbols), 50):
        batch = symbols[i:i + 50]
        query = " OR ".join(f"gene_exact:{s}" for s in batch) + " AND organism_id:559292 AND reviewed:true"
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": query, "fields": "accession,gene_names", "format": "json", "size": 500},
            timeout=60,
        )
        if resp.status_code != 200:
            print(f"  batch {i}: UniProt search failed {resp.status_code}", flush=True)
            continue
        for entry in resp.json().get("results", []):
            acc = entry.get("primaryAccession")
            genes = entry.get("genes", [])
            for g in genes:
                name = (g.get("geneName") or {}).get("value")
                if name:
                    for sgd_id, sym in symbol_by_sgd.items():
                        if sym == name:
                            uniprot_by_sgd[sgd_id] = acc
    print(f"  resolved {len(uniprot_by_sgd)}/{len(sgd_ids)} SGD genes to a real UniProt accession", flush=True)

    # ---- 3. Real OrthoDB gene id for each UniProt accession (from gene_xrefs) -
    needed_accs = set(uniprot_by_sgd.values())
    orthodb_by_uniprot = {}
    print(f"scanning real OrthoDB gene_xrefs for {len(needed_accs)} UniProt accessions", flush=True)
    n = 0
    with requests.get(XREFS_URL, stream=True, timeout=1800) as resp:
        resp.raise_for_status()
        with gzip.GzipFile(fileobj=resp.raw) as gz:
            for raw in gz:
                n += 1
                line = raw.decode("utf-8", errors="ignore").rstrip("\n")
                parts = line.split("\t")
                if len(parts) != 3:
                    continue
                gene_id, xref, db = parts
                if db == "UniProt" and xref in needed_accs and xref not in orthodb_by_uniprot:
                    orthodb_by_uniprot[xref] = gene_id
                if n % 20_000_000 == 0:
                    print(f"  scanned {n} xref rows, matched {len(orthodb_by_uniprot)}/{len(needed_accs)}", flush=True)
                if len(orthodb_by_uniprot) == len(needed_accs):
                    break
    print(f"  matched {len(orthodb_by_uniprot)}/{len(needed_accs)} UniProt accessions to a real OrthoDB gene id "
          f"(scanned {n} rows)", flush=True)

    # ---- 4. Verify the OrthoDB gene id is actually IN the real orthology graph
    ortho_path = Path("/art/odb10v1_ascomycota_OGs.tab")
    ortho_gene_ids = set()
    with ortho_path.open() as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 3:
                ortho_gene_ids.add(parts[1])
    print(f"  real orthology graph has {len(ortho_gene_ids)} distinct gene ids", flush=True)

    sgd_to_ortho = {}
    for sgd_id, acc in uniprot_by_sgd.items():
        ortho_id = orthodb_by_uniprot.get(acc)
        if ortho_id and ortho_id in ortho_gene_ids:
            sgd_to_ortho[sgd_id] = ortho_id
    print(f"  FULL CHAIN VERIFIED for {len(sgd_to_ortho)}/{len(sgd_ids)} SGD genes "
          f"(SGD -> UniProt -> real orthology-graph node)", flush=True)

    # ---- 5. Rewrite the GAF file's gene ids to the real bridged orthology ids
    bridged_lines = []
    for line in gaf_text.splitlines():
        if line.startswith("!") or not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 5 or cols[4] != GO_TERM:
            continue
        sgd_id = cols[1]
        if sgd_id in sgd_to_ortho:
            cols = list(cols)
            cols[1] = sgd_to_ortho[sgd_id]
            bridged_lines.append("\t".join(cols))
    bridged_gaf = Path("/data_bridged.gaf")
    bridged_gaf.write_text("\n".join(bridged_lines) + "\n")
    print(f"  wrote {len(bridged_lines)} bridged annotation rows -> {bridged_gaf}", flush=True)

    # ---- 6. Real build_graph + evaluation with the bridged ids --------------
    from data_pipeline.build_graph import main as build_main

    out_path = "/art/graph_bridged.pt"
    argv = [
        "build_graph.py",
        "--orthology-path", str(ortho_path),
        "--go-annotations-path", str(bridged_gaf),
        "--go-term", GO_TERM,
        "--target-taxon", "4896",
        "--output", out_path,
    ]
    old_argv = sys.argv
    sys.argv = argv
    try:
        build_main()
    finally:
        sys.argv = old_argv

    from model.evaluation import main as eval_main
    import io as _io
    import contextlib
    sys.argv = ["evaluation.py", "--graph", out_path, "--source-taxon", "559292", "--target-taxon", "4896"]
    buf = _io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            eval_main()
    finally:
        sys.argv = old_argv
    eval_stdout = buf.getvalue()
    print(eval_stdout, flush=True)
    Path("/art/eval_bridged_stdout.txt").write_text(eval_stdout)
    vol.commit()

    print("RESULT:", eval_stdout[:3000], flush=True)
    return {
        "n_go_annotated": len(sgd_ids),
        "n_uniprot_resolved": len(uniprot_by_sgd),
        "n_orthodb_resolved": len(orthodb_by_uniprot),
        "n_full_chain_verified": len(sgd_to_ortho),
        "eval_stdout": eval_stdout,
    }


@app.local_entrypoint()
def main():
    import json
    print(json.dumps(bridge_and_train.remote(), indent=2, default=str))
