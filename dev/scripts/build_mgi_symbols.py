#!/usr/bin/env python
"""
Build the bundled MGI mouse symbol table from MGI's marker report.

Download MRK_List2.rpt from https://www.informatics.jax.org/downloads/reports/
(current markers only; gzipped is fine), then, from the repo root:

    python dev/scripts/build_mgi_symbols.py MRK_List2.rpt.gz

Writes kosmic/reference/mgi/mgi_symbols.tsv.gz with one row per current
gene or pseudogene: 'Symbol', 'Locus group' (MGI's feature type reduced to
the four HGNC locus groups, so the Gene Names step counts both species the
same way) and 'Synonyms' (pipe-separated, as MGI lists them; they include
the gene's earlier symbols). Other marker types (QTLs, DNA segments,
transgenes, BAC ends) are left out: they are not in an expression matrix.

Also writes mgi_symbols.json beside it: the date the report was downloaded
(the file's modification date unless --downloaded is given), the build
date and the number of symbols. The Gene Names step records these in the
study's provenance, so the methods say which MGI release was used.
"""

import argparse
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd

OUT = Path("kosmic/reference/mgi/mgi_symbols.tsv.gz")
MARKER_TYPES = ("Gene", "Pseudogene")


def locus_group(feature_type: str) -> str:
    """MGI feature type -> HGNC locus group."""
    if feature_type == "protein coding gene":
        return "protein-coding gene"
    if "pseudogen" in feature_type:
        return "pseudogene"
    if "RNA gene" in feature_type:
        return "non-coding RNA"
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("report", help="MRK_List2.rpt or MRK_List2.rpt.gz")
    parser.add_argument("--out", default=str(OUT))
    parser.add_argument("--downloaded", default=None,
                        help="date the report was downloaded, YYYY-MM-DD "
                             "(default: the file's modification date)")
    args = parser.parse_args()
    downloaded = args.downloaded or date.fromtimestamp(
        Path(args.report).stat().st_mtime).isoformat()

    mrk = pd.read_csv(args.report, sep="\t", dtype=str).fillna("")
    keep = mrk["Marker Type"].isin(MARKER_TYPES) & (mrk["Status"] == "O")
    table = pd.DataFrame({
        "Symbol": mrk.loc[keep, "Marker Symbol"],
        "Locus group": mrk.loc[keep, "Feature Type"].map(locus_group),
        "Synonyms": mrk.loc[keep, "Marker Synonyms (pipe-separated)"],
    }).sort_values("Symbol")
    if table["Symbol"].duplicated().any():
        sys.exit("MGI report lists a symbol twice; check the input file.")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    # mtime 0: the same report always gives byte-identical output, so the
    # SHA-1 recorded in provenance changes only when the table does.
    table.to_csv(out, sep="\t", index=False,
                 compression={"method": "gzip", "mtime": 0})
    version = {
        "source": "MGI MRK_List2.rpt",
        "url": "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt",
        "downloaded": downloaded,
        "built": date.today().isoformat(),
        "n_symbols": int(len(table)),
    }
    out.with_name(out.name.replace(".tsv.gz", ".json")).write_text(
        json.dumps(version, indent=2) + "\n", encoding="utf-8")
    print(f"MGI report downloaded {downloaded}: {len(table):,} symbols "
          f"({table['Locus group'].value_counts().to_dict()}), "
          f"{(table['Synonyms'] != '').sum():,} with synonyms -> {out}")


if __name__ == "__main__":
    main()
