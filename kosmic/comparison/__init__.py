"""Comparison of merged pseudobulk with meta-analysis (ADR-017).

The Methods Comparison workspace compares the genes called by one DESeq2
model fitted to every study's donors in the shared atlas (``_master``,
merged pseudobulk) with the genes called by REML, SumRank, Stouffer and
vote counting on the per-study results. This package holds the parts
that do not depend on Qt:

- ``design``: the merged model's design, built from the per-study DE run
  records, and the checks that it matches them.
- ``merged``: the atlas's saved merged pseudobulk run, or a new fit of it.
- ``ma``: REML, SumRank, Stouffer and vote counting, read from a saved run
  or pooled again with one set of settings.
- ``standard``: one long table for every method's result.
"""
