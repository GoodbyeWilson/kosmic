"""
Fetch metadata for SumRank 16 COVID-19 PBMC datasets from GEO
and merge with CellxGene Census results.

No special dependencies beyond requests + pandas.

Usage:
    conda activate kirk
    python simulations/query_covid_pbmc_datasets.py   # produces results/covid_pbmc_suitable.csv
    python simulations/merge_covid_datasets.py         # merges and deduplicates

Output:
    simulations/results/covid_pbmc_geo.csv
    simulations/results/covid_pbmc_combined.csv
"""

import sys
import os
import time
import json
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>

import pandas as pd

# -- SumRank 16 COVID-19 PBMC datasets --
# Source: Nakatsuka et al. Nature Communications 2025, Supplementary Data
SUMRANK_DATASETS = [
    {"reference": "Su_Cell2020", "geo": "GSE150728", "first_author": "Su",
     "journal": "Cell", "year": 2020, "covid_donors": 129, "control_donors": 16,
     "assay": "10x Genomics (multi-ome)", "note": "Multi-ome, flagged as lower quality by SumRank"},
    {"reference": "SchulteSchrepping_Cell2020", "geo": "EGAS00001004571",
     "first_author": "Schulte-Schrepping", "journal": "Cell", "year": 2020,
     "covid_donors": 27, "control_donors": 38, "assay": "10x Genomics",
     "note": "EGA controlled access"},
    {"reference": "Yu_CellResearch2020", "geo": "GSE155673", "first_author": "Yu",
     "journal": "Cell Research", "year": 2020, "covid_donors": 7, "control_donors": 3,
     "assay": "10x Genomics", "note": "Small N"},
    {"reference": "Zhu_Immunity2020", "geo": "GSE150861", "first_author": "Zhu",
     "journal": "Immunity", "year": 2020, "covid_donors": 5, "control_donors": 3,
     "assay": "DNBelab/MGI", "note": "Small N, non-10x assay"},
    {"reference": "Liao_NatMed2020", "geo": "GSE145926", "first_author": "Liao",
     "journal": "Nature Medicine", "year": 2020, "covid_donors": 9, "control_donors": 4,
     "assay": "10x Genomics", "note": "Small N"},
    {"reference": "Trump_NatBiotech2020", "geo": "GSE149689", "first_author": "Trump",
     "journal": "Nature Biotechnology", "year": 2020, "covid_donors": 32, "control_donors": 16,
     "assay": "10x Genomics", "note": ""},
    {"reference": "Wen_CellDiscovery2020", "geo": "GSE158055", "first_author": "Wen",
     "journal": "Cell Discovery", "year": 2020, "covid_donors": 10, "control_donors": 5,
     "assay": "10x Genomics", "note": "Small N"},
    {"reference": "Lee_SciImmunology2020", "geo": "GSE149689", "first_author": "Lee",
     "journal": "Science Immunology", "year": 2020, "covid_donors": 11, "control_donors": 4,
     "assay": "10x Genomics", "note": "Small N"},
    {"reference": "Wilk_NatMed2020", "geo": "GSE150728", "first_author": "Wilk",
     "journal": "Nature Medicine", "year": 2020, "covid_donors": 7, "control_donors": 6,
     "assay": "Seq-Well", "note": "Non-10x assay"},
    {"reference": "Arunachalam_Science2020", "geo": "GSE155673", "first_author": "Arunachalam",
     "journal": "Science", "year": 2020, "covid_donors": 7, "control_donors": 5,
     "assay": "10x CITE-seq", "note": "Small N"},
    {"reference": "Combes_Nature2021", "geo": "GSE174072", "first_author": "Combes",
     "journal": "Nature", "year": 2021, "covid_donors": 20, "control_donors": 14,
     "assay": "10x CITE-seq", "note": ""},
    {"reference": "Stephenson_NatMed2021", "geo": "GSE174072", "first_author": "Stephenson",
     "journal": "Nature Medicine", "year": 2021, "covid_donors": 86, "control_donors": 23,
     "assay": "10x CITE-seq", "note": "Largest dataset"},
    {"reference": "Bacher_Immunity2020", "geo": "GSE161918", "first_author": "Bacher",
     "journal": "Immunity", "year": 2020, "covid_donors": 14, "control_donors": 6,
     "assay": "10x CITE-seq", "note": ""},
    {"reference": "Chua_NatBiotech2020", "geo": "GSE145926", "first_author": "Chua",
     "journal": "Nature Biotechnology", "year": 2020, "covid_donors": 19, "control_donors": 5,
     "assay": "10x CITE-seq", "note": "Small control N"},
    {"reference": "Kusnadi_SciImmunology2021", "geo": "GSE161918", "first_author": "Kusnadi",
     "journal": "Science Immunology", "year": 2021, "covid_donors": 37, "control_donors": 9,
     "assay": "10x CITE-seq", "note": ""},
    {"reference": "Meckiff_Cell2020", "geo": "GSE161918", "first_author": "Meckiff",
     "journal": "Cell", "year": 2020, "covid_donors": 37, "control_donors": 9,
     "assay": "10x CITE-seq", "note": ""},
]

MIN_CASES = 6
MIN_CONTROLS = 6


def fetch_geo_title(geo_acc):
    """Fetch dataset title from NCBI GEO via Entrez API."""
    if not geo_acc.startswith("GSE"):
        return ""
    try:
        url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
               f"?db=gds&term={geo_acc}[Accession]&retmode=json")
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            ids = json.loads(resp.read().decode()).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return ""

        time.sleep(0.35)

        url2 = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=gds&id={ids[0]}&retmode=json")
        req2 = urllib.request.Request(url2)
        with urllib.request.urlopen(req2, timeout=10) as resp2:
            result = json.loads(resp2.read().decode()).get("result", {})
        entry = result.get(ids[0], {})
        return entry.get("title", "")
    except Exception as e:
        print(f"    GEO fetch failed for {geo_acc}: {e}")
        return ""


def build_geo_table():
    """Build summary table from SumRank 16 hardcoded metadata."""
    print("Building SumRank 16 dataset table...")
    rows = []
    for d in SUMRANK_DATASETS:
        rows.append({
            "source": "GEO/EGA (SumRank 16)",
            "reference": d["reference"],
            "first_author": d["first_author"],
            "journal": d["journal"],
            "year": d["year"],
            "accession": d["geo"],
            "covid_donors": d["covid_donors"],
            "control_donors": d["control_donors"],
            "assay": d["assay"],
            "note": d["note"],
            "title": "",
        })

    df = pd.DataFrame(rows)

    print("Fetching GEO titles (~30 seconds due to rate limiting)...")
    seen = {}
    for i, row in df.iterrows():
        acc = row["accession"]
        if acc not in seen:
            print(f"  Fetching {acc}...", flush=True)
            seen[acc] = fetch_geo_title(acc)
            time.sleep(0.35)
        df.at[i, "title"] = seen[acc]

    return df


def load_cellxgene_results():
    """Load CellxGene results if available."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'results', 'covid_pbmc_suitable.csv')
    if not os.path.exists(path):
        print(f"  {path} not found -- skipping CellxGene merge.")
        print("  Run query_covid_pbmc_datasets.py first to include CellxGene datasets.")
        return None

    cxg = pd.read_csv(path)
    cxg["source"] = "CellxGene"
    cxg["accession"] = cxg["dataset_id"]
    cxg["reference"] = cxg.get("title", cxg["dataset_id"])
    cxg["first_author"] = ""
    cxg["journal"] = ""
    cxg["year"] = ""
    cxg["note"] = cxg.get("assay_flag", "")
    return cxg


def flag_suitability(df):
    df = df.copy()
    df["suitable"] = (
        (df["covid_donors"] >= MIN_CASES) &
        (df["control_donors"] >= MIN_CONTROLS)
    )
    reasons = []
    for _, row in df.iterrows():
        r = []
        if row["covid_donors"] < MIN_CASES:
            r.append(f"only {row['covid_donors']} COVID donors")
        if row["control_donors"] < MIN_CONTROLS:
            r.append(f"only {row['control_donors']} controls")
        reasons.append("; ".join(r) if r else "")
    df["exclusion_reason"] = reasons
    return df


def print_combined_summary(df):
    print("\n" + "=" * 70)
    print("COMBINED DATASET INVENTORY")
    print("=" * 70)

    suitable = df[df["suitable"]]
    excluded = df[~df["suitable"]]

    print(f"\nSUITABLE (>={MIN_CASES} COVID, >={MIN_CONTROLS} controls): "
          f"{len(suitable)} datasets")
    print("-" * 70)
    cols = ["source", "first_author", "year", "covid_donors",
            "control_donors", "assay", "note"]
    cols = [c for c in cols if c in suitable.columns]
    print(suitable[cols].to_string(index=False))

    print(f"\nEXCLUDED (below threshold): {len(excluded)} datasets")
    print("-" * 70)
    ecols = ["source", "first_author", "year",
             "covid_donors", "control_donors", "exclusion_reason"]
    ecols = [c for c in ecols if c in excluded.columns]
    print(excluded[ecols].to_string(index=False))

    print(f"\nTotals across suitable datasets:")
    print(f"  COVID-19 donors: {suitable['covid_donors'].sum()}")
    print(f"  Control donors:  {suitable['control_donors'].sum()}")
    print(f"  Unique datasets: {len(suitable)}")


def main():
    print("COVID-19 PBMC Dataset Inventory")
    print("=" * 70)
    print(f"Merging SumRank 16 (GEO/EGA) with CellxGene Census results\n")

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)

    # Build GEO table
    geo_df = build_geo_table()
    geo_path = os.path.join(out_dir, "covid_pbmc_geo.csv")
    geo_df.to_csv(geo_path, index=False)
    print(f"\nGEO table saved: {geo_path}")

    # Load CellxGene results
    cxg_df = load_cellxgene_results()

    # Merge
    if cxg_df is not None:
        keep_cols = ["source", "reference", "first_author", "journal",
                     "year", "accession", "covid_donors", "control_donors",
                     "assay", "note", "title"]
        geo_keep = [c for c in keep_cols if c in geo_df.columns]
        cxg_keep = [c for c in keep_cols if c in cxg_df.columns]

        combined = pd.concat(
            [geo_df[geo_keep], cxg_df[cxg_keep]],
            ignore_index=True
        )
        print(f"Combined: {len(geo_df)} GEO + {len(cxg_df)} CellxGene = "
              f"{len(combined)} total")

        combined["dup_key"] = (
            combined["first_author"].str.lower().fillna("") + "_" +
            combined["year"].astype(str)
        )
        dup_mask = combined.duplicated(subset=["dup_key"], keep=False)
        combined["possible_duplicate"] = dup_mask
        n_dups = dup_mask.sum()
        if n_dups:
            print(f"  {n_dups} rows flagged as possible duplicates -- review manually")
    else:
        combined = geo_df.copy()
        combined["possible_duplicate"] = False

    combined = flag_suitability(combined)
    combined = combined.sort_values(
        ["suitable", "control_donors", "covid_donors"],
        ascending=[False, False, False]
    ).reset_index(drop=True)

    combined_path = os.path.join(out_dir, "covid_pbmc_combined.csv")
    combined.to_csv(combined_path, index=False)
    print(f"Combined inventory saved: {combined_path}")

    print_combined_summary(combined)

    print("\nNext steps:")
    print("  1. Review results/covid_pbmc_combined.csv")
    print("  2. Resolve possible_duplicate=True rows")
    print("  3. Download h5ad files from CellxGene or GEO")
    print(f"  4. Recommended: use suitable datasets only "
          f"(>={MIN_CASES} COVID, >={MIN_CONTROLS} controls)")


if __name__ == "__main__":
    main()
