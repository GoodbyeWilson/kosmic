"""
Query CellxGene REST API for COVID-19 PBMC scRNA-seq datasets
suitable for pseudobulk meta-analysis.

No special dependencies -- uses urllib + json only.

Usage:
    conda activate kirk
    python simulations/query_covid_pbmc_datasets.py

Output:
    simulations/results/covid_pbmc_datasets.csv
    simulations/results/covid_pbmc_suitable.csv
"""

import sys
import os
import json
import urllib.request
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>

import pandas as pd

# -- Filtering criteria --
MIN_CASES    = 6
MIN_CONTROLS = 6
ASSAY_FILTER = ["10x 3' v2", "10x 3' v3", "10x 5' v1", "10x 5' v2",
                "10x 3' transcription profiling"]

API_BASE = "https://api.cellxgene.cziscience.com"


def _extract_labels(field):
    """Safely extract labels from a field that may be str, list of dicts, or list of str."""
    if isinstance(field, str):
        return field
    if isinstance(field, list):
        labels = []
        for item in field:
            if isinstance(item, dict):
                labels.append(item.get("label", str(item)))
            else:
                labels.append(str(item))
        return "; ".join(labels)
    return str(field)


def fetch_json(url, timeout=60):
    """Fetch JSON from URL with retry."""
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            if attempt < 2:
                print(f"  Retry {attempt+1}/3: {e}", flush=True)
                time.sleep(2)
            else:
                raise


def query_collections():
    """Fetch all collections from CellxGene."""
    print("Fetching collections from CellxGene API...", flush=True)
    url = f"{API_BASE}/curation/v1/collections"
    collections = fetch_json(url)
    print(f"  Found {len(collections)} collections", flush=True)
    return collections


def query_datasets():
    """Fetch all datasets from CellxGene."""
    print("Fetching datasets from CellxGene API...", flush=True)
    url = f"{API_BASE}/curation/v1/datasets"
    datasets = fetch_json(url, timeout=120)
    print(f"  Found {len(datasets)} datasets", flush=True)
    return datasets


def filter_covid_blood(datasets):
    """Filter datasets to COVID-19 + blood/PBMC."""
    print("Filtering for COVID-19 blood/PBMC datasets...", flush=True)

    blood_terms = {"blood", "peripheral blood", "pbmc",
                   "peripheral blood mononuclear cell"}
    covid_terms = {"COVID-19", "SARS-CoV-2"}

    results = []
    for ds in datasets:
        # Check disease
        diseases = set()
        for d in ds.get("disease", []):
            label = d.get("label", "")
            diseases.add(label)

        has_covid = bool(diseases & covid_terms)
        has_normal = "normal" in diseases

        if not has_covid and not has_normal:
            continue

        # Check tissue
        tissues = set()
        for t in ds.get("tissue", []):
            label = t.get("label", "").lower()
            tissues.add(label)

        is_blood = bool(tissues & blood_terms)
        if not is_blood:
            continue

        # Extract metadata
        assays = [a.get("label", "") for a in ds.get("assay", [])]
        cell_types = [c.get("label", "") for c in ds.get("cell_type", [])]
        organisms = [o.get("label", "") for o in ds.get("organism", [])]

        if "Homo sapiens" not in organisms:
            continue

        # Donor counts from dataset metadata
        # CellxGene API doesn't give per-condition donor counts directly,
        # so we use cell_count and disease labels as proxy
        results.append({
            "dataset_id": ds.get("dataset_id", ""),
            "collection_id": ds.get("collection_id", ""),
            "title": ds.get("title", ""),
            "assay": "; ".join(assays),
            "diseases": "; ".join(sorted(diseases)),
            "tissues": "; ".join(sorted(tissues)),
            "cell_count": ds.get("cell_count", 0),
            "n_cell_types": len(cell_types),
            "donor_id_count": len(ds.get("donor_id", [])),
            "has_covid": has_covid,
            "has_normal": has_normal,
            "has_both": has_covid and has_normal,
            "suspension_type": _extract_labels(ds.get("suspension_type", "")),
            "explorer_url": ds.get("explorer_url", ""),
        })

    print(f"  COVID/blood datasets: {len(results)}", flush=True)
    return pd.DataFrame(results)


def enrich_with_collection_info(df, collections):
    """Add collection title and DOI."""
    coll_map = {c.get("collection_id", ""): c for c in collections}

    coll_titles = []
    coll_dois = []
    for cid in df["collection_id"]:
        coll = coll_map.get(cid, {})
        coll_titles.append(coll.get("name", ""))
        dois = coll.get("doi", "")
        if not dois:
            links = coll.get("links", [])
            for link in links:
                if link.get("link_type") == "DOI":
                    dois = link.get("link_url", "")
                    break
        coll_dois.append(dois)

    df["collection_title"] = coll_titles
    df["doi"] = coll_dois
    return df


def apply_suitability_filter(df):
    """Filter to datasets likely suitable for pseudobulk MA."""
    # Must have both COVID and normal in same dataset
    suitable = df[df["has_both"]].copy()

    # Minimum cell count as proxy for having enough donors
    # (6 donors x ~1000 cells each = ~6000 minimum)
    suitable = suitable[suitable["cell_count"] >= 5000].copy()

    # Flag assays
    suitable["assay_flag"] = suitable["assay"].apply(
        lambda a: "OK" if any(x in a for x in ASSAY_FILTER) else "CHECK"
    )

    suitable = suitable.sort_values("cell_count", ascending=False).reset_index(drop=True)
    return suitable


def print_summary(suitable):
    print("\n" + "=" * 70)
    print(f"SUITABLE DATASETS (both COVID + normal, >= 5000 cells)")
    print("=" * 70)
    print(f"Found {len(suitable)} datasets\n")

    for i, row in suitable.iterrows():
        title = row["title"][:60] + "..." if len(row["title"]) > 60 else row["title"]
        print(f"  {i+1}. {title}")
        print(f"     Cells: {row['cell_count']:,}  |  "
              f"Cell types: {row['n_cell_types']}  |  "
              f"Assay: {row['assay']}  [{row['assay_flag']}]")
        print(f"     Diseases: {row['diseases']}")
        if row["collection_title"]:
            print(f"     Collection: {row['collection_title']}")
        if row["doi"]:
            print(f"     DOI: {row['doi']}")
        print()

    total_cells = suitable["cell_count"].sum()
    print(f"Total cells across suitable datasets: {total_cells:,}")


def main():
    print("CellxGene COVID-19 PBMC Dataset Query (REST API)")
    print("=" * 70)
    print()

    # Fetch data
    datasets = query_datasets()
    collections = query_collections()

    # Filter
    df = filter_covid_blood(datasets)

    if len(df) == 0:
        print("No matching datasets found.")
        return

    # Enrich with collection info
    df = enrich_with_collection_info(df, collections)

    # Save full inventory
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)

    df.to_csv(os.path.join(out_dir, "covid_pbmc_datasets.csv"), index=False)
    print(f"\nFull inventory: results/covid_pbmc_datasets.csv ({len(df)} datasets)")

    # Filter suitable
    suitable = apply_suitability_filter(df)
    suitable.to_csv(os.path.join(out_dir, "covid_pbmc_suitable.csv"), index=False)
    print(f"Suitable datasets: results/covid_pbmc_suitable.csv ({len(suitable)} datasets)")

    print_summary(suitable)

    print("\nNext steps:")
    print("  1. Review results/covid_pbmc_suitable.csv")
    print("  2. Cross-reference with SumRank paper's 16 datasets")
    print("  3. Download h5ad files from explorer_url links")
    print("  4. Run KOSMIC pseudobulk DE + meta-analysis pipeline")


if __name__ == "__main__":
    main()
