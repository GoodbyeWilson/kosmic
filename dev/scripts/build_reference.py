#!/usr/bin/env python
"""
Build a lightweight reference JSON from an annotated h5ad atlas.

Usage:
    python scripts/build_reference.py input.h5ad output.json.gz \
        --cell-type-col cell_type --name "My Atlas" --species human
"""

import argparse
import sys
sys.path.insert(0, '.')  # run from the repo root

from kosmic.scrna.annotate.reference import build_reference_from_h5ad


def main():
    parser = argparse.ArgumentParser(description='Build lightweight reference from annotated h5ad')
    parser.add_argument('input', help='Path to annotated h5ad file')
    parser.add_argument('output', help='Output path (.json.gz)')
    parser.add_argument('--cell-type-col', default='cell_type', help='Column with cell type labels')
    parser.add_argument('--name', default=None, help='Reference name')
    parser.add_argument('--species', default='human', choices=['human', 'mouse'])
    parser.add_argument('--source', default='', help='Citation / DOI')
    parser.add_argument('--n-genes', type=int, default=2000, help='Number of HVGs to include')
    parser.add_argument('--min-cells', type=int, default=10, help='Min cells per type')
    args = parser.parse_args()

    print(f"Building reference from {args.input}...")
    meta = build_reference_from_h5ad(
        args.input, args.output,
        cell_type_col=args.cell_type_col,
        name=args.name,
        species=args.species,
        source=args.source,
        n_top_genes=args.n_genes,
        min_cells_per_type=args.min_cells,
    )

    print(f"Reference built:")
    print(f"  Name: {meta['name']}")
    print(f"  Species: {meta['species']}")
    print(f"  Cell types: {meta['n_cell_types']}")
    print(f"  Genes: {meta['n_genes']}")
    print(f"  Saved to: {args.output}")


if __name__ == '__main__':
    main()
