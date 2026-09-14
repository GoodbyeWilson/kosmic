# Select Gene Sets

Choose which pathways or gene sets to test for differential expression.

## Built-in gene sets

KOSMIC includes curated metabolic pathway gene sets. 'Comprehensive Metabolic' covers 13 major metabolic pathways (~240 genes). 'Kirk Core' is a smaller focused set of 6 pathways.

## Custom gene sets

You can import your own gene sets from GMT, JSON, or CSV files. GMT is the standard format used by MSigDB and Enrichr. You can also search Enrichr's 225+ gene set libraries directly.

## Coverage table

The coverage table shows how many genes from each pathway are present in your dataset. Low coverage (<50%) may indicate species mismatch or that genes were filtered during QC. KOSMIC automatically converts between human and mouse gene naming conventions.
