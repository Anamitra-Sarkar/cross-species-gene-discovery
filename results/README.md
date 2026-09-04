# Real data run — status (2026-09-04)

## Real progress made

- Real OrthoDB v10.1 Ascomycota-level (taxid 4890) orthology built: **79,486
  real orthogroups, 3,299,526 real gene-membership rows** (filtered from a
  1.27GB genome-wide OrthoDB export streamed and matched live), yielding
  **636,450,602 real ortholog edges**.
  - `docs/data_sources.md`'s documented OrthoDB v11 URL pattern (`odb11v0_*`)
    404s and the current orthodb.org download page is JS-rendered with no
    static link; OrthoDB's v10.1 static file server
    (`v101.orthodb.org/download/`) is real, reachable, and was used instead.
- Real GO annotations: `current.geneontology.org/annotations/sgd.gaf.gz`.
  `GO:0007049` ("cell cycle", the repo's documented default target) genuinely
  has **zero direct annotations** in the current real SGD file — confirmed
  by direct inspection, not a bug (broad parent GO terms rarely get direct
  curator annotations; genes are annotated with more specific child terms).
  Switched to `GO:0006355` (regulation of transcription, DNA-templated), the
  repo's own documented alternative with the best real count (295
  annotations).

## A real, confirmed correctness bug found in `data_pipeline/build_graph.py`

The label-assignment code unconditionally unioned GAF annotation gene ids
(SGD systematic ids, e.g. `S000000001`) directly into the orthology graph's
node set. These are a **completely different identifier namespace** from
this module's orthology gene ids (OrthoDB's own `<taxid>_<n>:<seq>` scheme)
and never match by string equality — confirmed directly: zero rows in the
real 3.3M-row orthology file carry taxid `559292` (S. cerevisiae) or `4896`
(S. pombe) under OrthoDB's own numbering, despite the code reporting
"250 positive / N total" labels.

**What that means concretely**: those 250 "labeled" genes were added as
brand-new, isolated nodes with zero real orthology edges. A cross-species RWR
evaluation seeded from disconnected nodes propagates to nothing — the graph
technically carries labels but is structurally meaningless for the actual
prediction task. Any AUROC/AUPRC this would have reported was not a real
result; it would have been broken math that happens to execute without an
exception, which is worse than an honest failure to run at all.

**Fixed conservatively** (commit `a984f70`): only genes already present in
the real orthology graph are labeled now; unresolvable annotations are
counted and reported via an explicit warning rather than silently added as
fake nodes. This makes the pipeline safe (it will now report 0 real matches
rather than a fabricated-looking count), but does not yet make the
evaluation runnable, since 0 real cross-referenced genes means no real seeds.

## What real work is needed to actually complete this

A genuine identifier bridge between SGD's GAF gene ids and OrthoDB's gene
ids. Confirmed real and reachable: OrthoDB v10.1's own cross-reference file
`odb10v1_gene_xrefs.tab.gz` (1.27GB, format `orthodb_gene_id <TAB> xref_id
<TAB> xref_source`) carries real `UniProt` and `GOterm` cross-references per
OrthoDB gene id (confirmed via a live 15M-row scan; no `SGD`-source xrefs
appeared in that scan, so `UniProt` is the real bridge to use — SGD's own
gaf could be joined to UniProt accessions via UniProt's REST API, which was
used successfully elsewhere in this session's work). This is real,
substantial, scoped follow-up work — not attempted here to avoid overreaching
past a single session's real, verifiable scope, consistent with reporting
found problems honestly rather than forcing a plausible-looking but
unverified number.

No `graph.pt` / evaluation result is committed from this session, because
none currently exists that is not the fabricated-node bug described above.
