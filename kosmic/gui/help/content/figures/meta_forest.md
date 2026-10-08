# Meta-analysis forest plot

One gene (**Gene Forest Plot**) or pathway (**Pathway Forest Plot**)
across the studies of the latest meta-analysis.

- Each study's row shows its log2 fold change (diamond) and 95%
  confidence interval (grey line), from that study's own DE results.
  Diamonds are red when up in the disease group and blue when down; their
  size reflects the study's random-effects weight, 1 / (SE² + tau²).
- **Pooled** is the meta-analysis estimate, with its 95% confidence
  interval as the diamond's width. The same interval is shaded behind
  every row: red or blue when it excludes zero, grey when it crosses zero.
- The dashed line marks zero. **X limit** fixes a symmetric axis range;
  *Auto* fits it to the data.

The **Gene / pathway** list starts with the most significant. Exporting
also writes a text file next to the figure with the pooled estimate, FDR,
heterogeneity (I² and tau²) and each study's log2 fold change, SE, 95%
confidence interval and weight.
