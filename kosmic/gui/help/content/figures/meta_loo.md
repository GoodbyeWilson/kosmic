# Leave-one-out sensitivity

Whether one gene's or pathway's pooled result depends on any single
study. Each row leaves one study out and shows the pooled log2 fold
change and 95% confidence interval of the remaining studies. The solid
line and shaded band are the estimate and interval with all studies.

The studies are re-pooled with the method the meta-analysis recorded;
when it recorded a consensus of several methods, REML is used. At least
three studies with a result for the item are needed. **Level** chooses
gene- or pathway-level results.

Exporting also writes a text file next to the figure with the
all-studies result and, for each study left out, the re-pooled estimate,
confidence interval and I².
