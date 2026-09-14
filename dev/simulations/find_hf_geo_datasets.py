"""
Query NCBI GEO for heart failure scRNA-seq datasets.

No special dependencies -- uses urllib + pandas only.

Usage:
    conda activate kirk
    python simulations/find_hf_geo_datasets.py

Output:
    simulations/results/hf_geo_datasets.csv
    simulations/results/hf_geo_suitable.csv
"""

import sys
import os
import json
import time
import urllib.request
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # repo root: dev/<this dir>/<file>

import pandas as pd

# -- Search terms --
# Search the GEO database directly (db=geo) for Series entries
# More comprehensive than gds which only indexes curated DataSets
DISEASE_TERMS = [
    "heart failure",
    "cardiomyopathy",
    "dilated cardiomyopathy",
    "hypertrophic cardiomyopathy",
    "ischemic cardiomyopathy",
    "cardiac",
    "myocardial infarction",
    "non-failing heart",
]

SC_TERMS = [
    "single cell",
    "single-cell",
    "scRNA",
    "scRNA-seq",
    "single nucleus",
    "single-nucleus",
    "snRNA",
    "snRNA-seq",
    "10x genomics",
    "10x chromium",
    "Drop-seq",
]

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RATE_LIMIT_SLEEP = 0.35


def fetch_json(url, timeout=15):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def fetch_text(url, timeout=15):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def build_queries():
    """Build search queries combining disease + single-cell terms."""
    queries = []
    for disease in DISEASE_TERMS:
        for sc in SC_TERMS:
            q = f'"{disease}" AND "{sc}" AND "Homo sapiens"[Organism]'
            queries.append(q)
    return queries


def search_geo(query, db="gds", retmax=500):
    """Search GEO via Entrez."""
    url = (f"{BASE_URL}esearch.fcgi"
           f"?db={db}&term={urllib.parse.quote(query)}"
           f"&retmax={retmax}&retmode=json")
    try:
        result = fetch_json(url).get("esearchresult", {})
        ids = result.get("idlist", [])
        count = result.get("count", 0)
        return ids, int(count)
    except Exception as e:
        print(f"  Search failed: {e}")
        return [], 0


def fetch_summaries_gds(ids, batch_size=50):
    """Fetch eSummary from gds database."""
    records = {}
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        id_str = ",".join(batch)
        url = f"{BASE_URL}esummary.fcgi?db=gds&id={id_str}&retmode=json"
        try:
            result = fetch_json(url, timeout=20).get("result", {})
            for uid, rec in result.items():
                if uid != "uids":
                    records[uid] = rec
        except Exception as e:
            print(f"  Fetch failed for batch {i}: {e}")
        time.sleep(RATE_LIMIT_SLEEP)
    return records


def search_geo_series_directly():
    """Search GEO Series via the geo database (not gds).

    Uses esearch on db=geo with type filter for Series, then
    converts UIDs to GSE accessions via esummary.
    """
    print("\n--- Searching GEO Series database ---", flush=True)

    all_ids = set()
    queries = build_queries()

    # Also add broad queries
    broad = [
        '"heart failure" AND "single cell RNA" AND "Homo sapiens"[Organism]',
        '"heart failure" AND "RNA sequencing" AND "single" AND "Homo sapiens"[Organism]',
        '"cardiomyopathy" AND "single cell" AND "Homo sapiens"[Organism]',
        '"non-failing" AND "single cell" AND heart AND "Homo sapiens"[Organism]',
        '"failing heart" AND "single cell" AND "Homo sapiens"[Organism]',
        '"heart" AND "single-cell transcriptom" AND "Homo sapiens"[Organism]',
        '"heart" AND "single-nucleus transcriptom" AND "Homo sapiens"[Organism]',
        '"cardiac" AND "scRNA-seq" AND "Homo sapiens"[Organism]',
        '"cardiac" AND "snRNA-seq" AND "Homo sapiens"[Organism]',
    ]
    queries.extend(broad)

    for i, query in enumerate(queries):
        short = query[:65] + "..." if len(query) > 65 else query
        # Search both gds and geo databases
        for db in ["gds"]:
            ids, count = search_geo(query, db=db, retmax=500)
            if ids:
                all_ids.update(ids)
                print(f"  [{i+1}/{len(queries)}] {db}: +{len(ids)} ({count} total) | {short}", flush=True)
        time.sleep(RATE_LIMIT_SLEEP)

    print(f"\n  Total unique GEO entries: {len(all_ids)}", flush=True)
    return list(all_ids)


def parse_record(uid, rec):
    acc = rec.get("accession", "")
    if not acc.startswith("GSE"):
        return None

    entry_type = rec.get("entrytype", "")
    if entry_type not in ("GSE", ""):
        return None

    summary = rec.get("summary", "")
    title = rec.get("title", "")

    return {
        "geo_accession": acc,
        "title": title,
        "summary": summary[:500] + "..." if len(summary) > 500 else summary,
        "n_samples": rec.get("n_samples", 0),
        "taxon": rec.get("taxon", ""),
        "platform_id": rec.get("gpl", ""),
        "experiment_type": rec.get("ptype", ""),
        "date_published": rec.get("pdat", ""),
        "geo_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",
    }


def classify_tissue(title, summary):
    text = (title + " " + summary).lower()
    if any(t in text for t in ["pbmc", "peripheral blood", "blood"]):
        return "blood/PBMC"
    if any(t in text for t in ["cardiac", "heart", "myocardium", "myocardial",
                                "ventricle", "atrium", "cardiomyocyte"]):
        return "heart tissue"
    if any(t in text for t in ["endotheli"]):
        return "endothelium"
    return "other/unclear"


def flag_case_control(title, summary):
    text = (title + " " + summary).lower()
    cc_terms = ["control", "healthy", "normal", "non-failing", "nonfailing",
                "vs.", "versus", "compared to", "case", "donor",
                "disease", "patient", "failing", "ischemic", "dilated",
                "hcm", "dcm", "icm"]
    return any(t in text for t in cc_terms)


def flag_snrna(title, summary):
    text = (title + " " + summary).lower()
    if any(t in text for t in ["snrna", "sn-rna", "single-nucleus",
                                "single nucleus", "snseq", "nuclei"]):
        return "snRNA-seq"
    return "scRNA-seq"


def classify_source(title, summary):
    """Classify as primary tissue, cultured cells/iPSC, or unclear."""
    text = (title + " " + summary).lower()
    culture_terms = ["ipsc", "ips ", "induced pluripotent", "organoid", "culture",
                     "cultured", "cell line", "in vitro", "differentiat",
                     "hipsc", "stem cell-derived", "hesc", "cell-derived"]
    tissue_terms = ["biopsy", "surgical", "explant", "transplant", "autopsy",
                    "tissue", "patient", "donor heart", "left ventricle",
                    "myocardial", "ventricular", "atrial", "septum",
                    "intraoperative", "endomyocardial"]
    is_culture = any(t in text for t in culture_terms)
    is_tissue = any(t in text for t in tissue_terms)
    if is_culture and is_tissue:
        return "both"
    if is_culture:
        return "iPSC/culture"
    if is_tissue:
        return "primary tissue"
    return "unclear"


def fetch_gsm_characteristics(geo_acc, timeout=30):
    """Fetch per-sample characteristics from GEO SOFT format.

    Returns
    -------
    dict: {gsm_id: {char_key: char_value, ...}, ...}
    """
    if not geo_acc.startswith("GSE"):
        return {}
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={geo_acc}&targ=gsm&view=brief&form=text")
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode('utf-8', errors='replace')
    except Exception as e:
        print(f"    Failed to fetch GSM data for {geo_acc}: {e}")
        return {}

    samples = {}
    current_gsm = None
    for line in text.split('\n'):
        if line.startswith('^SAMPLE'):
            current_gsm = line.split('=')[1].strip()
            samples[current_gsm] = {}
        elif line.startswith('!Sample_characteristics_ch1') and current_gsm:
            val = line.split('=', 1)[1].strip() if '=' in line else ''
            key_val = val.split(': ', 1) if ': ' in val else [val, '']
            if len(key_val) == 2:
                samples[current_gsm][key_val[0].strip()] = key_val[1].strip()
    return samples


def classify_condition(value):
    """Classify a sample characteristic value as case, control, or unknown."""
    v = value.lower().strip()
    control_terms = ["normal", "non-diseased", "nonfailing", "non-failing",
                     "healthy", "control", "donor", "nf", "non-disease",
                     "sham", "unaffected"]
    case_terms = ["heart failure", "hf", "dcm", "dilated", "icm", "ischemic",
                  "hcm", "hypertrophic", "cardiomyopathy", "failing",
                  "disease", "patient", "mi", "infarction", "arvc"]
    if any(t in v for t in control_terms):
        return "control"
    if any(t in v for t in case_terms):
        return "case"
    return "unknown"


def count_case_control(gsm_chars):
    """Count case/control samples from GSM characteristics.

    Looks for disease-related characteristic keys and classifies values.

    Returns
    -------
    n_case, n_control, n_unknown, condition_key, condition_values
    """
    if not gsm_chars:
        return 0, 0, 0, "", {}

    # Find the characteristic key most likely to be disease/condition
    # Priority tiers: exact disease keys > broader keys > value-based fallback
    from collections import Counter

    # Tier 1: keys that explicitly mention disease
    tier1_terms = ["disease state", "disease status", "disease", "diagnosis",
                   "pathology", "clinical diagnosis", "patient diagnosis"]
    # Tier 2: broader keys (but exclude misleading ones like "storage condition")
    tier2_terms = ["condition", "status", "group", "phenotype", "type"]
    tier2_exclude = ["storage", "culture", "media", "experimental", "sample type",
                     "cell type", "tissue type"]

    tier1_keys = []
    tier2_keys = []
    for gsm, chars in gsm_chars.items():
        for k in chars:
            kl = k.lower().strip()
            if any(t == kl or t in kl for t in tier1_terms):
                tier1_keys.append(k)
            elif any(t in kl for t in tier2_terms):
                if not any(ex in kl for ex in tier2_exclude):
                    tier2_keys.append(k)

    disease_keys = tier1_keys if tier1_keys else tier2_keys

    if not disease_keys:
        # Tier 3: try all keys and see if any values look like conditions
        for gsm, chars in gsm_chars.items():
            for k, v in chars.items():
                kl = k.lower()
                if any(ex in kl for ex in ["cell type", "tissue", "storage",
                                            "technology", "genotype"]):
                    continue
                cl = classify_condition(v)
                if cl != "unknown":
                    disease_keys.append(k)

    if not disease_keys:
        return 0, 0, len(gsm_chars), "", {}

    # Use the most common disease key
    key_counts = Counter(disease_keys)
    best_key = key_counts.most_common(1)[0][0]

    # Count conditions
    condition_values = {}
    n_case, n_control, n_unknown = 0, 0, 0
    for gsm, chars in gsm_chars.items():
        val = chars.get(best_key, "")
        if not val:
            n_unknown += 1
            continue
        condition_values[val] = condition_values.get(val, 0) + 1
        cl = classify_condition(val)
        if cl == "case":
            n_case += 1
        elif cl == "control":
            n_control += 1
        else:
            n_unknown += 1

    return n_case, n_control, n_unknown, best_key, condition_values


def enrich_with_gsm_data(df, max_datasets=None):
    """Fetch GSM characteristics for suitable datasets and add case/control counts."""
    print(f"\nFetching GSM-level metadata for {len(df)} datasets...", flush=True)
    print("  (this queries GEO per dataset, ~2s each due to rate limiting)\n", flush=True)

    n_case_list = []
    n_control_list = []
    n_unknown_list = []
    condition_key_list = []
    condition_detail_list = []

    to_process = df if max_datasets is None else df.head(max_datasets)

    for i, (idx, row) in enumerate(to_process.iterrows()):
        acc = row["geo_accession"]
        print(f"  [{i+1}/{len(to_process)}] {acc}...", end=" ", flush=True)

        gsm_chars = fetch_gsm_characteristics(acc)
        n_case, n_control, n_unknown, cond_key, cond_vals = count_case_control(gsm_chars)

        n_case_list.append(n_case)
        n_control_list.append(n_control)
        n_unknown_list.append(n_unknown)
        condition_key_list.append(cond_key)
        detail = "; ".join(f"{v}: {c}" for v, c in sorted(cond_vals.items(), key=lambda x: -x[1]))
        condition_detail_list.append(detail)

        print(f"case={n_case} ctrl={n_control} unk={n_unknown}  "
              f"[{cond_key}] {detail[:60]}", flush=True)

        time.sleep(RATE_LIMIT_SLEEP)

    to_process = to_process.copy()
    to_process["n_case_donors"] = n_case_list
    to_process["n_control_donors"] = n_control_list
    to_process["n_unknown_donors"] = n_unknown_list
    to_process["condition_key"] = condition_key_list
    to_process["condition_detail"] = condition_detail_list

    return to_process


def main():
    print("GEO Heart Failure scRNA-seq Dataset Search (comprehensive)")
    print("=" * 60)
    print(f"Disease terms: {len(DISEASE_TERMS)}")
    print(f"SC terms: {len(SC_TERMS)}")
    print()

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)

    # Comprehensive search
    all_ids = search_geo_series_directly()

    if not all_ids:
        print("No results found.")
        return

    # Fetch metadata
    print(f"\nFetching metadata for {len(all_ids)} entries...", flush=True)
    records = fetch_summaries_gds(all_ids)
    print(f"Retrieved {len(records)} records")

    # Parse
    rows = []
    for uid, rec in records.items():
        parsed = parse_record(uid, rec)
        if parsed:
            rows.append(parsed)

    if not rows:
        print("No GSE records found.")
        return

    df = pd.DataFrame(rows)

    # Classify
    df["tissue"] = df.apply(lambda r: classify_tissue(r["title"], r["summary"]), axis=1)
    df["likely_case_ctrl"] = df.apply(lambda r: flag_case_control(r["title"], r["summary"]), axis=1)
    df["seq_type"] = df.apply(lambda r: flag_snrna(r["title"], r["summary"]), axis=1)
    df["sample_source"] = df.apply(lambda r: classify_source(r["title"], r["summary"]), axis=1)

    df = df.sort_values("n_samples", ascending=False).reset_index(drop=True)

    # Save full results
    full_path = os.path.join(out_dir, "hf_geo_datasets.csv")
    df.to_csv(full_path, index=False)
    print(f"\nFull results saved: {full_path} ({len(df)} GSE entries)")

    # Filter suitable -- primary tissue, case-control, reasonable N
    suitable = df[
        (df["tissue"] == "heart tissue") &
        (df["likely_case_ctrl"] == True) &
        (df["n_samples"] >= 6) &
        (df["sample_source"].isin(["primary tissue", "both", "unclear"]))
    ].copy()

    suit_path = os.path.join(out_dir, "hf_geo_suitable.csv")
    suitable.to_csv(suit_path, index=False)
    print(f"Suitable datasets (pre-enrichment): {suit_path} ({len(suitable)} entries)")

    # Enrich primary tissue datasets with GSM-level case/control counts
    primary = suitable[suitable["sample_source"] == "primary tissue"].copy()
    if len(primary) > 0:
        enriched = enrich_with_gsm_data(primary)
        enriched_path = os.path.join(out_dir, "hf_geo_enriched.csv")
        enriched.to_csv(enriched_path, index=False)
        print(f"\nEnriched datasets saved: {enriched_path} ({len(enriched)} entries)")

        # Show top candidates
        good = enriched[(enriched["n_case_donors"] >= 3) & (enriched["n_control_donors"] >= 3)]
        print(f"\nDatasets with >= 3 case AND >= 3 control donors: {len(good)}")
        if len(good) > 0:
            for _, row in good.iterrows():
                print(f"  {row['geo_accession']:>12}  case={row['n_case_donors']:>3}  "
                      f"ctrl={row['n_control_donors']:>3}  "
                      f"{row['seq_type']:>10}  {str(row['title'])[:55]}")

    # Summary
    print("\n" + "=" * 60)
    print("ALL RESULTS SUMMARY")
    print("=" * 60)
    print(f"\nBy tissue type:")
    print(df["tissue"].value_counts().to_string())
    print(f"\nBy seq type:")
    print(df["seq_type"].value_counts().to_string())
    print(f"\nBy sample source:")
    print(df["sample_source"].value_counts().to_string())
    print(f"\nLikely case-control design: {df['likely_case_ctrl'].sum()}")

    print("\n" + "=" * 60)
    print(f"LIKELY SUITABLE (heart tissue, case-control, n>=6 samples)")
    print("=" * 60)
    if len(suitable) == 0:
        print("No datasets met all criteria -- review hf_geo_datasets.csv manually")
    else:
        display_cols = ["geo_accession", "seq_type", "sample_source", "n_samples",
                        "date_published", "title"]
        # Truncate titles and handle unicode for Windows console
        display = suitable[display_cols].copy()
        display["title"] = display["title"].str[:80].str.encode('ascii', 'replace').str.decode('ascii')
        print(display.to_string(index=False))

    print("\nNext steps:")
    print("  1. Review results/hf_geo_suitable.csv")
    print("  2. Check GEO pages for HF donor vs control donor counts")
    print("  3. Note: snRNA-seq is standard for heart tissue")
    print("  4. Look for datasets with raw counts + cell type annotations")


if __name__ == "__main__":
    main()
