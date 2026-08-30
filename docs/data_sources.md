# Data Sources

This document describes the real, public, no-auth data sources used (or designed to be used) in this project. All endpoints and formats are real and correctly cited. No fabricated datasets.

## 1. Orthology: OrthoDB (primary) and eggNOG (alternative)

### OrthoDB (chosen primary)

- **Website:** https://www.orthodb.org/
- **Citation:** Kuznetsov D. et al. "OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity." *Nucleic Acids Res.* 2023;51(D1):D445-D451. https://doi.org/10.1093/nar/gkac998
- **REST API base:** `https://www.orthodb.org`
- **Key endpoints (real, documented at https://www.orthodb.org/?page=api):**
  - `GET /blast?seq=<fasta>&level=<clade>&species=<taxid>` — find ortholog group for a query sequence
  - `GET /ogdetails?id=<OrthoDB_OG_ID>` — genes in an ortholog group
  - `GET /search?query=<gene>&level=<clade>&species=<taxid>` — search ortholog groups
  - `GET /fasta?id=<OG_ID>&species=<taxid>` — sequences for a group/species
  - `GET /tab?id=<OG_ID>` — tabular ortholog group content
- **Bulk downloads (real):**
  - `https://www.orthodb.org/download/odb11v0_OGs.tar.gz` — all OGs
  - `https://www.orthodb.org/download/odb11v0_gene_xrefs.tar.gz` — gene cross-references
  - `https://www.orthodb.org/download/odb11v0_level2species.tab.gz` — level/species mapping
  - Per-level OGs: `https://www.orthodb.org/download/odb11v0_${level}_OGs.tab.gz` e.g. `33208` (Eukaryota), `2759` (Eukaryota subset)
- **Taxonomy IDs used in this project:**
  - *Saccharomyces cerevisiae* (budding yeast): NCBI TaxID `559292` (strain S288C) / `4932`
  - *Mus musculus* (mouse): `10090`
  - *Schizosaccharomyces pombe* (fission yeast, example target): `4896`
  - *Danio rerio* (zebrafish, example target): `7955`
- **License:** CC BY 4.0
- **Access:** No authentication required.

### eggNOG 5 (documented alternative)

- **Website:** http://eggnog5.embl.de/
- **Citation:** Huerta-Cepas J. et al. "eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource." *Nucleic Acids Res.* 2019;47(D1):D309-D314.
- **Bulk downloads:** `http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/` — orthology groups per taxonomic level.
- **Access:** No authentication required.

### Format used by `data_pipeline/`

The pipeline expects a **TSV** orthology edge file (matching OrthoDB bulk TAB format simplified):

```
# OrthoDB-like edge list (synthetic fixture matches this)
# Columns: og_id  gene_id_1  species_1  gene_id_2  species_2  score
OG_EUK_0001  SGD:S000000001  559292  POMBASE:SPAC1F8.01  4896  0.95
OG_EUK_0001  SGD:S000000002  559292  ZFIN:ZDB-GENE-0001  7955  0.82
```

- `score` is an orthology confidence (0-1). For OrthoDB, this can be derived from within-group sequence similarity or treated as 1.0 for unweighted graphs. The pipeline supports weighted and unweighted modes.
- The real OrthoDB TAB file has columns like `OrthoDB OG ID, gene ID, organism, UniProt, etc.` — the pipeline's parser handles both the simplified edge list and the native OrthoDB `*_OGs.tab` format (see `data_pipeline/orthology.py`).

**Real-run procedure (Orthology):**

```bash
# Example: download Eukaryota OGs for yeast + target species
curl -L https://www.orthodb.org/download/odb11v0_Eukaryota_OGs.tab.gz -o data/odb_Eukaryota_OGs.tab.gz
gunzip data/odb_Eukaryota_OGs.tab.gz
python -m data_pipeline.build_graph \
  --orthology-path data/odb_Eukaryota_OGs.tab \
  --go-annotations-path data/goa_yeast.gaf \
  --target-taxon 4896 \
  --output data/graph.pt
```

---

## 2. Functional Annotation: Gene Ontology (GO) and GO Annotation (GOA)

### Gene Ontology (GO)

- **Website:** http://geneontology.org/
- **Citation:** Gene Ontology Consortium. "The Gene Ontology knowledgebase in 2023." *Genetics* 2023;224(1):iyad031.
- **Ontology download (real):**
  - `http://current.geneontology.org/ontology/go-basic.obo` (OBO format)
  - `http://current.geneontology.org/ontology/go.obo`
  - `http://purl.obolibrary.org/obo/go.owl`
- **License:** CC BY 4.0

### GO Annotations (GOA) — GAF format

- **Downloads index:** http://current.geneontology.org/products/pages/downloads.html
- **GOA project:** https://www.ebi.ac.uk/GOA/
- **Real GAF files (no auth):**
  - Yeast SGD: `https://current.geneontology.org/annotations/sgd.gaf.gz` (alias at `http://current.geneontology.org/annotations/sgd.gaf.gz`)
  - Mouse MGI: `https://current.geneontology.org/annotations/mgi.gaf.gz`
  - PomBase (S. pombe): `https://current.geneontology.org/annotations/pombase.gaf.gz`
  - Zebrafish ZFIN: `https://current.geneontology.org/annotations/zfin.gaf.gz`
  - UniProt GOA (all species, large): `https://current.geneontology.org/annotations/goa_uniprot_all.gaf.gz`
- **GAF 2.2 format spec:** https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
  - Tab-separated, 17 columns. Key columns for this project:
    - Col 2: DB Object ID (gene ID)
    - Col 3: DB Object Symbol
    - Col 5: GO ID (e.g. `GO:0007049`)
    - Col 7: Evidence code (e.g. `IDA`, `IMP`, `IEA`)
    - Col 9: Aspect (`P`=Biological Process, `F`=Molecular Function, `C`=Cellular Component)
    - Col 13: Taxon (`taxon:559292`)
- **Citation:** Huntley RP et al. "The GOA database: Gene Ontology annotation updates for 2015." *Nucleic Acids Res.* 2015;43(D1):D1057-D1063.

### GO Term used as prediction target (real, well-studied)

The pipeline defaults to a well-characterized Biological Process term with a meaningfully-sized yeast gene set:

- **`GO:0007049` — "cell cycle"** (Biological Process)
  - Hundreds of experimentally annotated yeast genes (SGD).
  - Well-studied, tractable binary classification target.
  - Alternative terms documented in `data_pipeline/go_annotations.py`:
    - `GO:0006468` — protein phosphorylation
    - `GO:0006355` — regulation of transcription, DNA-templated
    - `GO:0007005` — mitochondrion organization

GO Slim (generic subset) can be used for multi-label tractability:
- `https://current.geneontology.org/ontology/subsets/goslim_generic.obo`

**Real-run procedure (GO Annotations):**

```bash
curl -L https://current.geneontology.org/annotations/sgd.gaf.gz -o data/sgd.gaf.gz
gunzip data/sgd.gaf.gz
python -m data_pipeline.build_graph \
  --orthology-path data/odb_Eukaryota_OGs.tab \
  --go-annotations-path data/sgd.gaf \
  --go-term GO:0007049 \
  --evidence-codes IDA,IMP,IGI,IPI \
  --output data/graph.pt
```

Filtering by evidence codes (experimental vs electronic) is supported via `--evidence-codes`.

---

## 3. Synthetic Fixtures (for testing in sandbox)

Because real downloads are not performed in this sandbox (no heavy compute/downloads), the test suite uses **small synthetic fixture files** that match the real formats above:

- `tests/fixtures/orthology_edges.tsv` — 10-50 ortholog pairs across 2-3 synthetic species, with confidence scores
- `tests/fixtures/go_annotations.gaf` — minimal valid GAF 2.2 lines for source and target species
- `tests/fixtures/go-basic-subset.obo` — tiny OBO subset (optional)

These fixtures are explicitly synthetic and labeled as such. They use an **injected function-conservation signal** (most orthologs share GO term, a few diverged pairs do not) so the evaluation pipeline's ability to recover transferable signal can be verified.

---

## References (summary)

1. OrthoDB v11 — https://www.orthodb.org/ — Kuznetsov et al. NAR 2023.
2. eggNOG 5.0 — http://eggnog5.embl.de/ — Huerta-Cepas et al. NAR 2019.
3. Gene Ontology — http://geneontology.org/ — Gene Ontology Consortium, Genetics 2023.
4. GO Annotation File Format 2.2 — https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
5. SGD (yeast) — https://www.yeastgenome.org/
6. MGI (mouse) — http://www.informatics.jax.org/
