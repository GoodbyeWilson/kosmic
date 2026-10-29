# Gene Meta-Analysis (Consensus)

Pool per-study gene-level DE results with one or more methods. Single mode = one method; consensus mode = several methods, intersection of significance.

## Family checkboxes (Single / Consensus)

*TODO.*

### Vote counting

The **Vote counting** family passes a gene when at least *N* studies
call it significant (study FDR < 0.05) in the same direction and none
calls it significant in the opposite direction. **Vote threshold** (*N*
of *k*, all studies by default) appears when the family is ticked.
Switching to Consensus mode ticks every family except vote counting;
tick it to add the rule. Run with REML in Consensus mode, the consensus
genes are those with a pooled effect that at least *N* studies also
call on their own. It ignores effect sizes and gives no pooled
estimate, so use it beside DL or REML, not in place of them. Run on its
own, the page counts the genes that pass the rule rather than
significant genes. It has no p-value to draw a volcano from, so the
**Plots** tab shows the genes at each level of study agreement, split
up and down: the figure that Figure Export's
[Study Agreement](../figures/meta_vote_levels.md) page exports.

**Vote threshold** decides which genes the run calls significant: the
significant-gene count, the Venn diagram and the gene set the later
steps use. Every gene's level is kept whatever the threshold: the
results file has a `vote_level` column (for example `3/5` or
`conflict`), with `n_up`, `n_down` and `votes`.

In a run that includes vote counting, the **Results** tab shows, in
place of an FDR for vote counting, a column named after the threshold
(for example **≥ 5 of 5 agree**) that is ticked for the genes that
pass; the gene's
**Level** beside **k** (the studies in which the gene was tested) and a
**Studies agreeing** row with one box per level, 0/*k* to *k*/*k* and
*conflict*, each with its number of genes; untick levels to hide their
genes from the table. Figure Export's
[Study Agreement](../figures/meta_vote_levels.md) page draws the genes
at each level, split up and down, and can hatch the genes REML or
another method of the run calls significant. See
[Statistical approach](statistics.md#vote-counting).

Vote counting is new in the method list; it is enabled by default and
can be hidden in **View → Meta-Analysis Settings**.

## Top 50% and HKSJ modifiers

**Top 50% selection** changes how SumRank, DL and REML use each gene's
studies. It does not remove genes before pooling.

- With **SumRank**, each gene's score is the sum of its best ceil(k/2)
  study ranks (the best 3 of 5), tested against the null distribution of
  that partial sum. It asks whether a gene is highly ranked in at least
  half of the studies, which suits an effect present in a subset of
  cohorts. When the cohorts agree it is stricter than the full sum,
  because the other studies' evidence is not used: on five DCM
  cardiomyocyte studies it called 378 genes against 1,348 for plain
  SumRank.
- With **DL** and **REML**, each gene is pooled over its ceil(k/2) most
  precise studies (smallest standard error, at least 2).

**HKSJ correction** (Hartung-Knapp-Sidik-Jonkman) applies to DL and
REML. It replaces the normal test with a t-test on k - 1 degrees of
freedom, which is appropriate for one outcome pooled over a few
studies. Over thousands of genes it has little power when k is small:
with five studies a t of 20 gives p of about 4e-5, while
Benjamini-Hochberg over 11,000 genes needs about 4e-6 for the strongest
gene to pass. With fewer than 10 studies selected, a note under the box
gives the smallest reachable p-value, and a run that calls no gene says
so in the output panel.

The Methods page states each pooling method in words, including these
modifiers.

## Calibration (analytical / CC permutation)

*TODO.*

## Tabs: Plots / Pathway Forest / Results / Venn / Enrichment / Cross-Dataset Reproducibility / LOO Validation

*TODO.*

### Direction conflicts in the Results tab

The Results table has three columns that count, for each gene, the
studies in which it is significant (study FDR < 0.05) with a positive
log2 fold change (**Up**) and with a negative one (**Down**).
**Conflict** is marked when both counts are above zero: the gene is
significantly up in at least one study and significantly down in
another. The line above the table reports how many of the significant
genes have such a conflict, and the same count is written to the output
panel and the Methods page.

Tick **Show only genes with opposite-direction effects** to list only
these genes. They occur mainly with Fisher's method, SumRank and gwOP,
which do not use the sign of the effect; see
[Statistical approach](statistics.md#direction-of-effect-across-studies).

Studies whose DE results have no adjusted p-values are not counted; the
output panel names them when the run starts.

## Olink panel overlay

*TODO.*
