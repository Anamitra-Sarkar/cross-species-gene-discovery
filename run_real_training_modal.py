"""Run cross-species-gene-discovery's own pipeline on real OrthoDB + GO data.

Real, public, no-auth sources:
  - OrthoDB v10.1 bulk download server (v101.orthodb.org/download/) -- the URL the
    repo's own docs point at (orthodb.org's newer download page) is a JS-rendered
    SPA with no static link, and the documented odb11v0_* file pattern 404s. v10.1's
    static file server is real, reachable, and the parser (data_pipeline/orthology.py)
    accepts any OrthoDB TAB-format export.
  - odb10v1_OGs.tab.gz: og_id, level_taxid, description (113MB).
  - odb10v1_OG2genes.tab.gz: og_id, gene_id (gene_id embeds its NCBI taxid as a
    prefix, e.g. "559292_0:001234" -> taxid 559292; this is OrthoDB's own
    documented gene-id convention, not an invented mapping). 1.27GB total across
    all of life -- filtered here to Ascomycota (taxid 4890, the fungal phylum
    containing both the source and target species below) rather than pulling
    all of Eukaryota, since that is the orthology scope this task actually needs.
  - SGD (yeast) real GO annotations: current.geneontology.org/annotations/sgd.gaf.gz.
  - Target: GO:0006355 (regulation of transcription, DNA-templated -- GO:0007049 "cell cycle" is not obsolete but genuinely has zero DIRECT annotations in the real current SGD GAF, a real GO-annotation-specificity characteristic, not a bug; GO:0006355 is the repo's own documented alternative with the best real annotation count, 295), source taxon 559292 (S. cerevisiae), target
    taxon 4896 (S. pombe) -- the repo's own documented defaults.

Run: modal run train_csgd.py
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

app = modal.App("csgd-real-run", image=image)
vol = modal.Volume.from_name("csgd-artifacts", create_if_missing=True)

REPO = "https://github.com/Anamitra-Sarkar/cross-species-gene-discovery.git"
OGS_URL = "https://v101.orthodb.org/download/odb10v1_OGs.tab.gz"
OG2GENES_URL = "https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz"
GAF_URL = "https://current.geneontology.org/annotations/sgd.gaf.gz"
ASCOMYCOTA_LEVEL = "4890"


@app.function(timeout=7200, cpu=8.0, memory=32768, volumes={"/art": vol})
def train() -> dict:
    import gzip
    import json
    import subprocess
    import sys
    from pathlib import Path
    import requests

    subprocess.run(["git", "clone", "--depth", "1", REPO, "/repo"], check=True)
    sys.path.insert(0, "/repo")

    data = Path("/data"); data.mkdir(exist_ok=True)

    print(f"downloading OGs {OGS_URL}", flush=True)
    r = requests.get(OGS_URL, timeout=600)
    r.raise_for_status()
    ogs_text = gzip.decompress(r.content).decode("utf-8", errors="ignore")
    print(f"  {len(r.content)} bytes compressed, {len(ogs_text)} bytes decompressed", flush=True)

    ascomycota_ogs = set()
    for line in ogs_text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].strip() == ASCOMYCOTA_LEVEL:
            ascomycota_ogs.add(parts[0].strip())
    print(f"  {len(ascomycota_ogs)} Ascomycota-level (taxid {ASCOMYCOTA_LEVEL}) OGs", flush=True)
    if len(ascomycota_ogs) < 100:
        raise RuntimeError(f"Only {len(ascomycota_ogs)} Ascomycota OGs found -- level filter likely wrong.")

    cached_path = Path("/art/odb10v1_ascomycota_OGs.tab")
    filtered_path = data / "odb10v1_ascomycota_OGs.tab"
    if cached_path.exists():
        print(f"reusing cached filtered orthology file from a prior run: {cached_path} "
              f"({cached_path.stat().st_size} bytes)", flush=True)
        filtered_path = cached_path
        n_kept = sum(1 for _ in cached_path.open()) - 1
    else:
        print(f"streaming OG2genes {OG2GENES_URL} (1.27GB compressed, filtering to Ascomycota)", flush=True)
        n_rows = n_kept = 0
        with requests.get(OG2GENES_URL, stream=True, timeout=1800) as resp:
            resp.raise_for_status()
            with gzip.GzipFile(fileobj=resp.raw) as gz, filtered_path.open("w") as out:
                # Header the repo's own parser (data_pipeline/orthology.py) looks for.
                out.write("og_id\tgene_id\ttaxid\n")
                for raw_line in gz:
                    n_rows += 1
                    line = raw_line.decode("utf-8", errors="ignore").rstrip("\n")
                    og_id, _, gene_id = line.partition("\t")
                    if og_id not in ascomycota_ogs or not gene_id:
                        continue
                    # OrthoDB gene ids are documented as "<ncbi_taxid>_<n>:<seq>"; the
                    # leading numeric taxid is the real, stable part.
                    taxid = gene_id.split("_", 1)[0]
                    if not taxid.isdigit():
                        continue
                    out.write(f"{og_id}\t{gene_id}\t{taxid}\n")
                    n_kept += 1
                    if n_rows % 5_000_000 == 0:
                        print(f"  scanned {n_rows} rows, kept {n_kept}", flush=True)
        print(f"  done: scanned {n_rows} total OG2genes rows, kept {n_kept} Ascomycota rows -> {filtered_path}", flush=True)
        if n_kept < 1000:
            raise RuntimeError(f"Only kept {n_kept} Ascomycota gene rows -- something is wrong upstream.")
        # Cache for a possible re-run (e.g. trying a different --go-term).
        import shutil
        shutil.copy(filtered_path, cached_path)
        vol.commit()

    print(f"downloading GAF {GAF_URL}", flush=True)
    r = requests.get(GAF_URL, timeout=300)
    r.raise_for_status()
    gaf_path = data / "sgd.gaf"
    gaf_path.write_bytes(gzip.decompress(r.content))
    print(f"  {gaf_path} {gaf_path.stat().st_size} bytes", flush=True)

    from data_pipeline.build_graph import main as build_main

    out_path = "/art/graph.pt"
    argv = [
        "build_graph.py",
        "--orthology-path", str(filtered_path),
        "--go-annotations-path", str(gaf_path),
        "--go-term", "GO:0006355",
        "--target-taxon", "4896",
        "--output", out_path,
    ]
    old_argv = sys.argv
    sys.argv = argv
    try:
        build_main()
    finally:
        sys.argv = old_argv
    print(f"graph built: {Path(out_path).stat().st_size} bytes", flush=True)

    from model.evaluation import main as eval_main
    sys.argv = ["evaluation.py", "--graph", out_path, "--source-taxon", "559292", "--target-taxon", "4896"]
    import io as _io
    import contextlib
    buf = _io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            eval_main()
    finally:
        sys.argv = old_argv
    eval_stdout = buf.getvalue()
    print(eval_stdout, flush=True)
    Path("/art/eval_stdout.txt").write_text(eval_stdout)

    vol.commit()
    produced = sorted(str(p.relative_to("/art")) for p in Path("/art").rglob("*") if p.is_file())
    print("ARTIFACTS:", produced, flush=True)
    print("RESULT:", eval_stdout[:3000], flush=True)
    return {"artifacts": produced, "n_ascomycota_ogs": len(ascomycota_ogs), "n_ortholog_rows": n_kept, "eval_stdout": eval_stdout}


@app.local_entrypoint()
def main():
    import json
    print(json.dumps(train.remote(), indent=2, default=str))
