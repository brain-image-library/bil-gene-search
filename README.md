# BIL Gene Search Prototype

A prototype for browsing the [Brain Image Library (BIL)](https://www.brainimagelibrary.org/) spatial transcriptomics datasets by gene. Type any of the following into the search box and matching datasets across all available species will surface:

- Current official gene symbol (e.g., `Ebf1`)
- Historical or alternative names (aliases, previous symbols)
- NCBI Entrez gene ID (e.g., `7157`)
- Ensembl gene ID (e.g., `ENSG00000141510`, `ENSMUSG00000048332`)

Results are presented in the same dropdown:

- **Per-species hits** — direct matches, one entry per species where the gene is catalogued in BIL data.
- **Cross-species homolog cluster** — a single aggregated entry that spans every species sharing the gene's HomoloGene cluster, for cross-species comparison work.

This is a design prototype intended for UX validation and handoff. 
Data artefacts are based on the preliminary search for the spatial transcriptomics datasets in BIL.
Production implementation will run against a different stack, but the data artifacts produced by this repo are portable.

## Quick start

This project uses [uv](https://docs.astral.sh/uv/) for environment management.

```bash
# Install dependencies from the lockfile
uv sync

# Launch the app
uv run streamlit run app.py
```

The app opens at `http://localhost:8501`. 


## Regenerating the data artifacts

The app expects three pickle files in the working directory:

| File | Contents |
|---|---|
| `genes_resolved.pkl` | Gene resolution table — one row per gene per species, with canonical Entrez/Ensembl IDs, aliases, and HomoloGene cluster IDs. |
| `presence.pkl` | Gene-by-dataset presence matrix, keyed on Entrez ID. |
| `datasets.pkl` | Per-dataset metadata: title, species, location (anything that we want to display for the user). |

To regenerate from BIL's source data use `create_data_artefacts.ipynb`.

## Status

Prototype. Suitable for UX review and handoff. Not production code.
