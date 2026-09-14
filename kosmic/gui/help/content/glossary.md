# Glossary

*TODO: define key terms used throughout KOSMIC.*

### Pseudobulk
Aggregating per-cell counts up to one number per (sample, gene), so the unit of analysis is the donor / sample rather than the cell.

### Consensus
A gene called significant by every selected pooling method at FDR < 0.05.

### FDR
False discovery rate, controlled with the Benjamini-Hochberg procedure.

### log2FC
Log-2 fold change. Positive = up in disease, negative = down in disease.

### LOO Validation
Leave-one-out: rebuild the consensus on K-1 studies, test whether the held-out study replicates the call.

### CC permutation
Case/control label-permutation calibration of pooled p-values.

### HKSJ
Hartung-Knapp-Sidik-Jonkman: small-sample t-based variance correction for random-effects pooling.

### Top 50%
Pre-filter to the top half of genes by per-study p-value before pooling. Helps SumRank a lot, leaves DL essentially unchanged.

*(More terms to come.)*
