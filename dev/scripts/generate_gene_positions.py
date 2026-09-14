"""
Regenerate `kosmic/reference/gwas/gene_positions_GRCh37.tsv` from UCSC
RefFlat (hg19).

Run when you want to update the bundled gene-position reference. The
output is checked in so that runtime usage doesn't need network access.

Usage::

    python scripts/generate_gene_positions.py

The bundled file has columns: ``gene``, ``chrom``, ``start``, ``end``,
``strand``. One row per gene (collapsed across RefSeq isoforms by taking
the union of transcript bounds).
"""
from __future__ import annotations

import csv
import gzip
import io
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

UCSC_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz"
)
OUT_PATH = (
    Path(__file__).resolve().parents[2]
    / "kosmic" / "reference" / "gwas" / "gene_positions_GRCh37.tsv"
)


def main() -> None:
    print(f"Fetching {UCSC_URL} ...")
    with urllib.request.urlopen(UCSC_URL, timeout=60) as r:
        raw = r.read()
    print(f"  downloaded {len(raw) / 1e6:.1f} MB compressed")

    # Parse refFlat: geneName, name, chrom, strand, txStart, txEnd, ...
    records: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    with gzip.open(io.BytesIO(raw), "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            gene = parts[0].strip()
            chrom = parts[2].replace("chr", "")
            strand = parts[3]
            try:
                tx_s = int(parts[4])
                tx_e = int(parts[5])
            except ValueError:
                continue
            # Drop alt / random / hap contigs and mitochondria
            if "_" in chrom or chrom in ("Un", "M", "MT"):
                continue
            records[gene].append((chrom, strand, tx_s, tx_e))

    print(f"  unique gene names: {len(records)}")

    # Collapse per-gene; drop genes that map to >1 chromosome
    rows: list[tuple[str, str, int, int, str]] = []
    ambiguous = 0
    for gene, recs in records.items():
        chroms = {r[0] for r in recs}
        if len(chroms) > 1:
            ambiguous += 1
            continue
        chrom = next(iter(chroms))
        strand = recs[0][1]
        rows.append((
            gene, chrom,
            min(r[2] for r in recs),
            max(r[3] for r in recs),
            strand,
        ))
    print(f"  multi-chrom dropped: {ambiguous}")
    print(f"  final rows: {len(rows)}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "chrom", "start", "end", "strand"])
        for r in sorted(rows):
            w.writerow(r)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    sys.exit(main())
