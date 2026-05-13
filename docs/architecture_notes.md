# Gene Search Architecture and Biological Choices

This document captures the design decisions behind BIL's gene-search prototype. The goal is to make those choices explicit so that future reviewers and maintainers can understand *why* the system is structured the way it is — not just *what* it does.

## Definitions

**Gene symbol.** A short, human-readable identifier for a gene (e.g., `TP53`, `Adra1a`). Symbols are set by species-specific nomenclature committees: HGNC for human, MGI for mouse, RGD for rat. Symbols can and do change over time as gene function is refined or family relationships are formalized.

**Alias / previous symbol.** Older or alternative names by which a gene is also known. When a gene is renamed, the prior names persist in the database as aliases so legacy literature stays findable. Example: `Lhfp` was the original name for the gene now called `Lhfpl6`.

**NCBI Entrez gene ID.** A stable, numeric identifier maintained by NCBI for each curated gene record. One Entrez ID per gene per species. Entrez IDs do not change over time and are suitable as a canonical primary key. Example: `7157` is the Entrez ID for human TP53.

**Ensembl gene ID.** A stable, alphanumeric identifier maintained by Ensembl / EMBL-EBI, of the form `ENSG00000141510` (human), `ENSMUSG00000059552` (mouse), `ENSMMUG00000045427` (rhesus macaque). Ensembl IDs are also stable, but **a single biological gene can have multiple Ensembl IDs** when it appears on more than one path through the assembled genome — alternate haplotype contigs, assembly patches, or pseudoautosomal regions all produce additional Ensembl gene entries for the same biological gene.

**HomoloGene cluster ID.** An NCBI-curated identifier that groups orthologous genes across species into a single cluster. For example, human TP53, mouse Trp53, and rhesus macaque TP53 share HomoloGene cluster `460`. HomoloGene has not been actively curated since approximately 2014, so coverage of recently catalogued genes (especially long non-coding RNAs and non-mammalian orthologs) is incomplete.

**Ortholog.** Genes in different species derived from a common ancestral gene through speciation, typically retaining similar function. Orthologs are what one wants to compare across species in cross-species transcriptomic studies.

**Surface form.** Any string a user might type to refer to a gene: a current symbol, an alias, a previous symbol, an Ensembl ID, an Entrez ID. The alias-resolution index maps every surface form to its canonical gene record.

**Alt-haplotype / patch sequence.** A region of the genome that exists in multiple versions across different individuals, included in the reference assembly as separate sequence contigs. Genes located in these regions receive a separate Ensembl ID per haplotype/patch. The most extreme example in human is the MHC region on chromosome 6, where HLA genes can have eight or more Ensembl IDs each.

**mygene** / mygene.info. A Python client (and underlying web service, maintained by the BioThings team at Scripps Research) that aggregates gene annotations from NCBI Gene, Ensembl, UniProt, HGNC, MGI, and other primary sources into a unified query API. Given any common surface form — symbol, alias, Ensembl ID, Entrez ID — it returns the merged gene record with canonical identifiers, aliases, and (where available) HomoloGene cluster membership. 

## Canonical identifier choice

### Internal canonical key: NCBI Entrez gene ID

This decision drives the rest of the architecture. Reasoning:

1. **Uniqueness.** Each biological gene maps to exactly one Entrez ID within a species. By contrast, a single gene can have multiple Ensembl IDs (alt haplotypes, patches, pseudoautosomal regions). Using Entrez as the join key gives a clean one-to-one mapping between rows of the resolved-gene table and biological entities.

2. **Stability.** Both Entrez and Ensembl IDs are stable, but Entrez tends to retain its ID across reassembly events that can cause Ensembl IDs to be retired or split.

3. **Cross-reference availability.** Essentially every well-curated Ensembl gene record carries an Entrez cross-reference, and the reverse is also typically true, making the two systems interoperable.

### Displayed canonical ID: the primary Ensembl gene ID

The user-facing display preferentially shows the primary Ensembl ID — the one located on the primary genome assembly (not an alt haplotype) and cross-referenced from NCBI. Reasoning: spatial transcriptomics workflows (cellxgene, gene panels from Vizgen / 10x Xenium / Slide-seq) overwhelmingly use Ensembl IDs as the primary gene identifier, so users expect to see them.

The architecture cleanly separates these concerns:

| Layer | Identifier | Why |
|---|---|---|
| Database join | Entrez | One stable scalar per gene |
| User-facing display | Primary Ensembl ID | Matches what users see in transcriptomics files |
| Search input | Any surface form | Maximum forgiveness for what the user types |

## Merging rules

When mygene returns multiple records for a single query, the resolution pipeline classifies the result and decides whether to merge them, keep them as a single record, or flag them as ambiguous.

### Mergeable case 1 — NCBI / Ensembl unmerged records ("soft duplicate")

**Pattern.** All returned records agree on symbol, name, and species. Each record carries only one of `{ensembl, entrezgene}` — one is Ensembl-sourced, the other is NCBI-sourced, and the cross-reference linking them is missing.

**Example.** `CASC6` (long non-coding RNA on chromosome 14) — two records, one keyed on `ENSG00000224944` (Ensembl only), one on Entrez `101929083` (NCBI only).

**Decision.** Merge into one record. This is a database integration artifact, not biological ambiguity. The merged record carries the union of identifiers. Eight such soft duplicates were collapsed during the initial resolution pass.

### Mergeable case 2 — Multi-Ensembl-ID genes (alt haplotype)

**Pattern.** All records agree on symbol, name, and species, but each carries a *different* Ensembl gene ID. One record typically also carries the Entrez ID; the others do not.

**Example.** `AKT3` — two Ensembl IDs (`ENSG00000117020` on the primary chromosome 1, `ENSG00000275199` on alt-haplotype contig `HSCHR1_3_CTG32_1`), both sharing Entrez `10000`. The same pattern affects `TRBC1`, `IGHA1`, the HLA family, and approximately 59 other human genes residing in haplotype-rich regions.

**Decision.** Merge into one record. The merged record's `ensembl` field is a list containing all Ensembl IDs; the `entrez` field holds the single Entrez ID. The primary Ensembl ID (the one paired with the Entrez record) is exposed as the canonical display ID; the alt-haplotype IDs are searchable but not displayed as separate options in the autocomplete.

### Non-mergeable case — real biological ambiguity

**Pattern.** Returned records disagree on either symbol or name. The query is an alias for two or more genuinely distinct genes.

**Example.** `Lhfp` (mouse) — historically the name for the gene now called `Lhfpl6`, but also a lingering alias on the paralog `Lhfpl1` due to family-naming history (when the LHFPL tetraspan family was formalized, the founder gene was renamed and its paralogs were named `-like 1`, `-like 2`, etc., creating a cross-alias that persists in MGI).

**Decision.** Do not merge. Flag the query as ambiguous; log it; do not add it to the canonical index. For user-typed input the UI presents both candidates and requires an explicit choice. For ambiguous data found in deposited BIL datasets, the data depositor is the authoritative source — they alone know which gene their probe actually targets.

In BIL's current data, exactly one such truly ambiguous case remained after the merge passes (`Lhfp` in a mouse dataset) and at the stage of this prototype was omitted.

## Species-aware presence index

The presence matrix is keyed on `(Entrez ID, dataset ID)`. A cell is `True` only when **both** of the following hold:

1. The original data deposited at BIL marks this gene as present in this dataset.
2. The species of the dataset matches the species of the gene record.

The second condition matters because BIL's original data is keyed on raw symbol strings, and some symbols overlap loosely across species. The species-aware index ensures that, for example, mouse `Trp53` (Entrez `22059`) is only ever attributed to mouse datasets — even if some human dataset happened to contain a raw string that loosely mapped to mouse Trp53 through aliases.

## Cross-species search via HomoloGene

For users interested in a gene across all available species in BIL, the system surfaces an aggregate suggestion keyed on the gene's HomoloGene cluster ID. When the user types any surface form that resolves to at least one gene with a known HomoloGene cluster, the system additionally suggests:

> `[all species] {Symbol_set} · {species_in_cluster} · {total_datasets} datasets`

This works even when the typed string only directly matches one species — for example, typing `TP53` only matches the human gene directly, but the human gene's HomoloGene cluster pulls mouse `Trp53` and macaque `TP53` into the aggregate suggestion. Clicking the aggregate suggestion shows the full per-species breakdown.

### Coverage caveats

- HomoloGene was last comprehensively updated by NCBI around 2014. Recently catalogued genes — especially long non-coding RNAs and many macaque genes — may lack a HomoloGene cluster assignment.
- In BIL's current data, approximately 14% of macaque genes (42 of 300) have no HomoloGene assignment. These genes are reachable only by direct symbol or ID search, not via cross-species aggregation.
- see Future work for ideas to make it more robust.

## Two-pass resolution strategy

Gene resolution uses two passes to ensure that exact symbol matches always take precedence over alias matches.

1. **Pass 1 — strict symbol scope.** Query mygene with `scopes="symbol"`. If the query string is a current official symbol for a gene, this pass resolves it cleanly. The majority of queries resolve here.

2. **Pass 2 — alias and cross-reference scopes.** For queries that miss pass 1, query with `scopes="alias,ensembl.gene,entrezgene,name"`. This handles older symbols, Ensembl IDs, and Entrez IDs.

Splitting into two passes prevents the failure mode where an alias-field match outscores a symbol-field match in mygene's relevance ranking. This is the source of the original `Adra1a` / `Adra1d` issue that motivated this architecture: querying `Adra1a` with all scopes simultaneously returned `Adra1d` as the top hit (because `Adra1a` is on Adra1d's alias list for historical-naming reasons), even though `Adra1a` is the canonical current symbol for a different gene entirely. The two-pass approach guarantees that current-symbol matches always win.

## Macaque-specific notes

*Macaca mulatta* has substantially less curated gene annotation than human or mouse, with several practical consequences:

- **Alias coverage is sparse.** Most macaque genes in BIL's resolved index carry zero aliases, or just the symbol echoed back as its own alias. **Implication for UX**: users searching for a macaque gene by anything other than its exact symbol, Ensembl ID, or Entrez ID will not find it via direct prefix matching. The HomoloGene cluster mechanism is the principal fallback.
- **Entrez cross-references are present for all 300 macaque genes** in BIL's current data — that surface remained fully populated even though mygene's macaque records lack the alias depth of human or mouse.
- **HomoloGene coverage is 86%** (258 of 300 macaque genes have a cluster ID). The remaining 14% are reachable only by direct symbol or ID search.

## Data flow summary

```
BIL raw gene-by-dataset presence (rows = raw symbols, mixed across species)
  │
  ├── resolve_genes_for_species()          [per species, calls mygene.info]
  │      │
  │      └── genes_resolved.pkl            (one row per gene per species)
  │
  ├── build_presence_by_entrez()
  │      │
  │      └── presence.pkl                  (rows = Entrez ID, cols = dataset ID)
  │
  └── BIL dataset metadata
         │
         └── datasets.pkl                  (rows = dataset ID, cols = title, species, url)

At app startup:
  load_genes()      ─┐
  load_presence()    ├── build_alias_index()  →  in-memory autocomplete index
  load_datasets()   ─┘
```

## What is *not* in scope for this prototype

- **Production-scale search infrastructure** (Elasticsearch, Typesense, SQLite FTS). The Pandas-based prefix search is adequate for the ~15,000 surface-form rows in current BIL data, but the artifacts produced (alias index, presence matrix) are directly portable to whichever full-text search backend the web team chooses.
- **Real-time updates.** The pipeline is offline batch; new BIL deposits would need to trigger recreation of data artefacts.

## Future work

- **Ensembl Compara** as a replacement for HomoloGene to improve coverage of recently catalogued and non-mammalian orthologs.
- **HCOP** (HGNC Comparison of Orthology Predictions) for human-centric, manually curated orthology data.
- **Schema requirement at deposit time.** Future BIL spatial transcriptomics deposits could require Ensembl gene IDs in the gene metadata (matching cellxgene's submission schema). This would eliminate the symbol-to-canonical resolution step entirely for newly submitted datasets, leaving only the legacy backlog to handle with mygene-based resolution.
