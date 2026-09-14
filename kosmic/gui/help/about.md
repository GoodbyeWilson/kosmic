## KOSMIC

*Kirk's Open-source Single-cell Meta-analysis Integration and Comparison*

We wanted to ask "is the expression of X gene reliably altered in publicly-available single-cell transcriptomic datasets". After lots of messing around in Jupyter Notebooks, Kaggle worksheets, migrated through Google Colab, we eventually thought our lives would be a lot easier if we just built a desktop application.

The key point in our software is to ask "what is the consensus" on whether the expression of a specific gene is altered in disease. This is a little different from how most single-cell datasets are approached. But it suits the work in the lab.

### Workspaces

KOSMIC is split into "workspaces" or tabs that provide an end-to-end meta-analysis pipeline.


  - **sc/snRNA-seq Analysis** - A single-dataset pipeline. Download from GEO, inspect what you got, run quality control, normalise, cluster, annotate cell types, then filter to the cells you actually care about.

  - **Differential Expression** - Once you've got annotated cells, you can perform pseudobulk differential expression analysis here. But will results from one study replicate in another?

  - **Meta-Analysis** - The whole point of KOSMIC. Pool DE results across multiple datasets to ask whether a gene is *reliably* altered, instead of altered in just one study. Random-effects models (DerSimonian-Laird, REML), classical p-value pooling (Fisher, Stouffer), and newer rank-based methods (SumRank, wOP, rOP). Has methods for assessing at the gene and pathway level.

  - **Figure Export** - We all want nice figures, so this tab creates them using Matplotlib.
