"""
BIL gene-search demo — unified search, Entrez-keyed presence.

Run with:
    pip install streamlit streamlit-searchbox pandas
    streamlit run app.py

Expected files in the working directory (pickle by default; change loaders to
your preferred format):

  genes_resolved.pkl      output of resolve_genes_for_species() concatenated
                          across species. Columns: species, gene_query,
                          ensembl (list[str]), entrez (list[str]),
                          aliases (list[str]), homolog_id (int|None).

  presence.pkl            output of build_presence_by_entrez().
                          index = Entrez ID (as str), columns = dataset_id,
                          values = bool.

  datasets.pkl            per-dataset metadata. index = dataset_id,
                          columns = title, species, url (extra cols ignored).
"""
from __future__ import annotations

import json
from typing import Any

import pandas as pd
import streamlit as st
from streamlit_searchbox import st_searchbox


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

@st.cache_data
def load_genes() -> pd.DataFrame:
    return pd.read_pickle("./data/genes_resolved.pkl")


@st.cache_data
def load_presence() -> pd.DataFrame:
    """Entrez-keyed presence matrix from build_presence_by_entrez()."""
    df = pd.read_pickle("./data/presence_by_entrez.pkl")
    # Normalize the index to string so JSON-encoded payloads round-trip safely
    df.index = df.index.astype(str)
    return df


@st.cache_data
def load_datasets() -> pd.DataFrame:
    return pd.read_pickle("./data/datasets.pkl")


# -----------------------------------------------------------------------------
# Alias-resolution index
#
# One row per (surface_form, gene). Surface forms include the symbol, every
# alias, every Ensembl ID, and the Entrez ID. Each row carries its canonical
# Entrez ID so the dataset lookup is a single hop from a search hit.
# -----------------------------------------------------------------------------

@st.cache_data
def build_alias_index(genes: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for _, g in genes.iterrows():
        canonical_entrez = g["entrez"][0] if g["entrez"] else None
        if canonical_entrez is None:
            # Skip genes with no Entrez — they can't be joined to the
            # entrez-keyed presence table. (Logged in your earlier QC pass.)
            continue
        base = {
            "species": g["species"],
            "gene_query": g["gene_query"],
            "display_symbol": g["gene_query"],
            "canonical_entrez": str(canonical_entrez),
            "homolog_id": g["homolog_id"],
        }

        def emit(form: str, mt: str) -> None:
            rows.append({**base, "surface_form": form, "match_type": mt})

        emit(g["gene_query"], "symbol")
        for a in (g["aliases"] or []):
            if a != g["gene_query"]:
                emit(a, "alias")
        for e in (g["ensembl"] or []):
            emit(e, "ensembl")
        for e in (g["entrez"] or []):
            emit(str(e), "entrez")

    df = pd.DataFrame(rows)
    df["surface_form_lower"] = df["surface_form"].str.lower()
    return df


@st.cache_data
def dataset_count_per_entrez(presence: pd.DataFrame) -> dict[str, int]:
    return presence.fillna(False).astype(bool).sum(axis=1).to_dict()


# -----------------------------------------------------------------------------
# Unified search — per-species hits plus an ortholog-cluster suggestion when
# ≥2 species in BIL share the gene's homolog cluster.
# Called per keystroke; do NOT cache.
# -----------------------------------------------------------------------------

def _per_species_label(row: pd.Series, n_datasets: int) -> str:
    sym = row["display_symbol"]
    sp = row["species"]
    if row["match_type"] == "symbol":
        return f"{sym} · {sp} · {n_datasets} datasets"
    return f"{row['surface_form']} → {sym} · {sp} · {n_datasets} datasets"


def search_unified(
    query: str,
    alias_idx: pd.DataFrame,
    counts: dict[str, int],
    limit: int = 10,
) -> list[tuple[str, str]]:
    if not query or len(query) < 2:
        return []
    q = query.lower()
    hits = alias_idx[alias_idx["surface_form_lower"].str.startswith(q)]
    if hits.empty:
        return []

    # Per-species hits + collect every homolog_id we encounter
    seen: set[tuple[str, str]] = set()
    per_species_hits: list[pd.Series] = []
    matched_homolog_ids: set[int] = set()
    for _, row in hits.iterrows():
        key = (row["species"], row["gene_query"])
        if key in seen:
            continue
        seen.add(key)
        per_species_hits.append(row)
        h = row["homolog_id"]
        if pd.notna(h):
            matched_homolog_ids.add(int(h))

    results: list[tuple[str, str]] = []
    for row in per_species_hits:
        n = counts.get(row["canonical_entrez"], 0)
        label = _per_species_label(row, n)
        value = json.dumps({
            "mode": "species",
            "species": row["species"],
            "gene_query": row["gene_query"],
        })
        results.append((label, value))

    # For every homolog_id we touched, look up ALL its members in BIL
    # (NOT just those matching the prefix). Show a cluster row if ≥2 species.
    for h in matched_homolog_ids:
        members = alias_idx[alias_idx["homolog_id"] == h].drop_duplicates(
            subset=["species", "gene_query"]
        )
        species_set = sorted(members["species"].unique())
        if len(species_set) < 2:
            continue
        total = sum(counts.get(r["canonical_entrez"], 0) for _, r in members.iterrows())
        symbol_label = " / ".join(sorted(members["display_symbol"].unique()))
        label = (f"[all species] {symbol_label} · "
                 f"{', '.join(species_set)} · {total} datasets")
        value = json.dumps({"mode": "homolog", "homolog_id": h})
        results.append((label, value))

    return results[:limit]


# -----------------------------------------------------------------------------
# Result rendering — all dataset lookups go through Entrez
# -----------------------------------------------------------------------------

def _list_datasets(ds_ids: list[str], datasets: pd.DataFrame) -> None:
    for ds_id in ds_ids:
        if ds_id in datasets.index:
            row = datasets.loc[ds_id]
            title = row.get("title", ds_id)
            url = row.get("url", None)
            if pd.notna(url):
                st.markdown(f"- **{title}** — location: {url}")
            else:
                st.markdown(f"- **{title}**")
        else:
            st.markdown(f"- {ds_id}")


def render_gene(
    species: str,
    gene_query: str,
    genes: pd.DataFrame,
    presence: pd.DataFrame,
    datasets: pd.DataFrame,
) -> None:
    rec_q = genes[(genes["species"] == species) & (genes["gene_query"] == gene_query)]
    if rec_q.empty:
        st.warning(f"No record for {gene_query} in {species}.")
        return
    rec = rec_q.iloc[0]
    st.subheader(f"{rec['gene_query']} · {species}")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Ensembl ID(s)**")
        for e in rec["ensembl"]:
            st.code(e)
    with c2:
        st.markdown("**Entrez ID**")
        for e in rec["entrez"]:
            st.code(e)
    if rec["aliases"]:
        st.markdown("**Aliases**: " + ", ".join(rec["aliases"]))
    if pd.notna(rec["homolog_id"]):
        st.markdown(f"**HomoloGene cluster**: `{int(rec['homolog_id'])}`")

    # Dataset lookup goes through ENTREZ now
    if not rec["entrez"]:
        st.info("No Entrez ID for this gene — cannot resolve datasets.")
        return
    entrez = str(rec["entrez"][0])
    if entrez in presence.index:
        ds_bool = presence.loc[entrez].fillna(False).astype(bool)
        ds_ids = ds_bool[ds_bool].index.tolist()
        st.markdown(f"### Found in {len(ds_ids)} dataset(s)")
        _list_datasets(ds_ids, datasets)
    else:
        st.info(f"Entrez `{entrez}` not present in any dataset.")


def render_homolog_group(
    homolog_id: int,
    genes: pd.DataFrame,
    presence: pd.DataFrame,
    datasets: pd.DataFrame,
) -> None:
    members = genes[genes["homolog_id"] == homolog_id]
    if members.empty:
        st.warning(f"No genes in ortholog cluster {homolog_id}.")
        return

    species_list = sorted(members["species"].unique())
    symbols = sorted(members["gene_query"].unique())
    st.subheader(f"Ortholog cluster — {' / '.join(symbols)}")
    st.caption(f"HomoloGene `{homolog_id}` · "
               f"{len(species_list)} species: {', '.join(species_list)}")

    # OR across all members' entrez rows for the combined dataset count
    member_entrez = [str(r["entrez"][0]) for _, r in members.iterrows() if r["entrez"]]
    available = [e for e in member_entrez if e in presence.index]
    if available:
        combined = presence.loc[available].any(axis=0)
        combined_ds = combined[combined].index.tolist()
        st.markdown(f"### Found in {len(combined_ds)} dataset(s) across the cluster")
    else:
        combined_ds = []

    # Per-species breakdown
    for sp in species_list:
        sp_members = members[members["species"] == sp]
        sp_entrez = [str(r["entrez"][0]) for _, r in sp_members.iterrows() if r["entrez"]]
        sp_entrez = [e for e in sp_entrez if e in presence.index]
        if sp_entrez:
            sp_bool = presence.loc[sp_entrez].any(axis=0)
            sp_ds_ids = sp_bool[sp_bool].index.tolist()
        else:
            sp_ds_ids = []
        with st.expander(f"{sp} — {len(sp_ds_ids)} dataset(s)", expanded=True):
            for _, r in sp_members.iterrows():
                ens = ", ".join(r["ensembl"]) if r["ensembl"] else "—"
                ent = r["entrez"][0] if r["entrez"] else "—"
                st.markdown(f"**{r['gene_query']}** — Entrez `{ent}` · Ensembl: {ens}")
            _list_datasets(sp_ds_ids, datasets)


# -----------------------------------------------------------------------------
# Main app
# -----------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(page_title="BIL Gene Search Prototype", layout="wide")
    st.title("BIL Spatial Transcriptomics — Gene Search Prototype")
    st.caption("Type a gene symbol, alias, Ensembl ID, or Entrez ID. "
               "Results span all species; cross-species homolog clusters are "
               "shown as a separate '[all species]' option when applicable.")

    genes = load_genes()
    presence = load_presence()
    datasets = load_datasets()
    alias_idx = build_alias_index(genes)
    counts = dataset_count_per_entrez(presence)

    selected: Any = st_searchbox(
        lambda q: search_unified(q, alias_idx, counts),
        placeholder="e.g. AKT3, p53, ENSG00000117020, 10000 …",
        key="gene_search",
    )

    if selected:
        try:
            payload = json.loads(selected)
        except (json.JSONDecodeError, TypeError):
            st.error("Unexpected selection payload.")
            return
        if payload.get("mode") == "species":
            render_gene(payload["species"], payload["gene_query"],
                        genes, presence, datasets)
        elif payload.get("mode") == "homolog":
            render_homolog_group(int(payload["homolog_id"]),
                                 genes, presence, datasets)


if __name__ == "__main__":
    main()
