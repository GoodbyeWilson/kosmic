# ADR-013: An unfiltered 10x matrix is imported as its cells, split into samples by barcode suffix

- Status: Approved (2026-10-30). Proposed 2026-10-27.

## Problem

Some GEO deposits contain every barcode of the 10x whitelist rather than
the cells Cell Ranger called: an unfiltered `raw_feature_bc_matrix`, or a
`cellranger aggr` run on raw matrices. GSE234052 is the second kind: one
`matrix.mtx.gz` of 32,285 genes by 40,769,280 barcodes, which is six
samples times the v3 whitelist (6,794,880 barcodes). About 22,000 of the
barcodes are cells.

Every 10x reader in `kosmic/scrna/load/converters.py` read the matrix
with `mmread` and created one `obs` row per barcode, so Convert to h5ad
on this deposit ran out of memory and the application stopped. The
sample of each barcode is given only by its suffix (`-1` .. `-6`), and the
converters assigned one sample per file or folder.

A downloaded deposit should be importable however its authors arranged
it, and an import that removes barcodes changes the contract that raw
counts are preserved, so the rule has to be explicit and recorded.

## Decision

1. **Large matrices are read in chunks.** On the import path
   (`cells_only=True` in `load_10x_folder` and `load_prefixed_10x`), a
   matrix with more than 200,000 barcode columns is read by
   `kosmic/scrna/load/large_10x.py`. A first pass totals the UMIs and
   detected genes of every barcode; a second pass builds the matrix of the
   kept barcodes only, so memory is bounded by the cells.
2. **An unfiltered matrix keeps only barcodes that look like cells.** When
   more than half of the barcodes have fewer than 100 UMIs, the file is
   taken to hold every whitelist barcode, and barcodes with at least 500
   UMIs and at least 200 detected genes are kept. A large matrix that does
   not meet that test is a filtered one and keeps every barcode. A fixed
   threshold rather than a knee or EmptyDrops call keeps the rule simple
   to state in the methods; on GSE234052 the drop in UMIs between empty
   droplets and cells is sharp, and the kept counts per sample (1,755 to
   7,204) do not depend on the exact cut between 200 and 1,000 UMIs.
3. **The rule is recorded.** `uns['kosmic_import']` holds the barcodes in
   the file, the filter applied, and the numbers kept and dropped; the
   import records it in the study's 'load' provenance stage.
4. **Aggregated barcodes become samples.** When the barcodes of one file
   carry more than one numeric suffix, each suffix becomes a sample named
   `{file prefix}_{N}`. The user names the samples and assigns conditions
   afterwards, as for any other import. A file that already defines its
   samples keeps them when several files are merged.
5. **Reading without `cells_only` is unchanged.** The SoupX raw-matrix read
   in QC needs the empty droplets and does not filter.

## Consequences

- GSE234052 converts in about a minute with a peak of about 2 GB of memory,
  to 22,203 cells in six samples.
- The empty droplets of such a deposit are not imported. The deposit file
  in `raw_data/` is unchanged and remains available for ambient RNA
  estimation outside KOSMIC.
- The 200,000-barcode limit is above the cell count of any single 10x
  sample, so ordinary filtered matrices take the existing path. A
  filtered matrix above the limit (a large aggregated study) is read in
  chunks and keeps every barcode.
- An aggregated filtered matrix below the limit is now split into samples
  by suffix as well, where it was previously imported as one sample.
