# KOSMIC

[![CI](https://github.com/GoodbyeWilson/kosmic/actions/workflows/ci.yml/badge.svg)](https://github.com/GoodbyeWilson/kosmic/actions/workflows/ci.yml)
[![Licence: GPL-3.0](https://img.shields.io/badge/licence-GPL--3.0-blue.svg)](LICENSE)

**Single-cell and single-nucleus RNA-seq analysis and cross-study
meta-analysis, without code.**

KOSMIC is a desktop application that takes a single-cell or
single-nucleus RNA-seq dataset — one you generated, or one of the
thousands deposited in public repositories — through quality control,
clustering, cell-type annotation, pseudobulk differential expression
and pathway analysis to publication-ready figures. When a question is
addressed by more than one dataset, KOSMIC pools their results in a
random-effects meta-analysis to find the genes and pathways that change
consistently across studies, and can build a shared cell-type atlas so
that the comparison is like for like. Every setting used is recorded
and rendered as a methods description.

![KOSMIC's scRNA Analysis workspace after clustering and annotation](kosmic/gui/help/content/scrna/img/cluster.png)

The name stands for *Kirk's Open-source Single-cell Meta-analysis
Integration and Comparison*.

## Key features

- **Complete single-dataset workflow** — import (`.h5ad`, Seurat, 10x,
  count matrices, or straight from GEO), QC, doublet removal,
  normalisation, clustering, annotation, ambient-RNA correction
- **Pseudobulk differential expression** with donors as the unit of
  replication (PyDESeq2), with covariate adjustment
- **Hypothesis-driven pathway analysis** — test the pathways you care
  about, across studies
- **Cross-study meta-analysis** — random-effects pooling of per-study
  results, with rank-based and p-value-combination methods alongside
  and permutation-calibrated significance
- **Consensus and validation** — the overlap between method families,
  leave-one-study-out reproducibility, comparison with proteomics panels
- **A shared atlas** — cluster and annotate several studies together and
  transfer the labels back
- **Publication figures**, and a **methods paragraph** generated from
  the recorded settings
- **No programming required**

## How KOSMIC works

A *project* is one biological question. Each dataset that addresses it
is added to the project as a *study* and processed on its own; the
per-study results are then combined.

```
Project
├── Study 1  ─  scRNA Analysis  ─  Differential Expression  ─┐
├── Study 2  ─  scRNA Analysis  ─  Differential Expression  ─┼─  Meta-Analysis  ─  Figures
└── Study 3  ─  scRNA Analysis  ─  Differential Expression  ─┘
```

The **Project** workspace manages the studies: adding them, importing
their data, and showing how far each has progressed. **scRNA Analysis**
takes one study from an expression matrix to annotated cell
populations, and can extract a cell type of interest as a study in its
own right. **Differential Expression** compares two conditions within a
cell type for one study. **Meta-Analysis** pools those results across
studies. **Figures** exports from any of them.

Studies are kept separate deliberately. A meta-analysis asks whether
independent studies agree; merging their cells into one experiment
would throw that question away. Where a common cell-type scheme is
needed, the studies can be combined into a shared atlas, annotated
once, and the labels transferred back to each study.

## Analysis modes

- **Hypothesis mode** — score a chosen pathway or gene programme in each
  study, pool the pathway-level effect across studies, and then look at
  which genes within it drive the signal. Multiple testing is over the
  pathways you chose, not the transcriptome.
- **Discovery mode** — genome-wide meta-analysis: which genes are
  consistently differentially expressed in this cell type across these
  datasets, followed by enrichment on the pooled result.
- **Methods comparison** — run every pooling method and evaluate them
  against each other by held-out replication.

The statistical reasoning behind the defaults is set out in the
documentation's [Statistical approach](kosmic/gui/help/content/meta/statistics.md).

## Getting started

**Installer.** One-click installers for Windows, macOS and Linux will be
published on the Releases page. Until the first release, install from
source.

**From source** — Python 3.11 or newer and Git are the prerequisites:

```sh
git clone https://github.com/GoodbyeWilson/kosmic.git
cd kosmic
python -m venv .venv
.venv\Scripts\activate          # Windows;  macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
python main.py
```

The step-by-step version, including installing Python, R for Seurat
files, and the optional extras, is in the
[installation guide](kosmic/gui/help/content/install.md).

## Documentation

The user guide — one page per screen, reachable with **F1** inside the
application — covers the [Project workspace](kosmic/gui/help/content/project/index.md),
[scRNA Analysis](kosmic/gui/help/content/scrna/index.md),
[Differential Expression](kosmic/gui/help/content/de/index.md),
[Meta-Analysis](kosmic/gui/help/content/meta/index.md) and
[Figures](kosmic/gui/help/content/figures/index.md), starting from
[what KOSMIC is](kosmic/gui/help/content/index.md). It will be published
at docs.scmetaanalysis.com.

Developers: [CONTRIBUTING.md](CONTRIBUTING.md) for setup and conventions,
[CODEBASE.md](CODEBASE.md) for the architecture, and a developer
reference generated from the source in the documentation site.

## Citation

A methods paper is in preparation. Until it is published, please cite
this repository.

## Licence

KOSMIC is released under the [GNU General Public License v3.0](LICENSE).
Its GUI toolkit (PyQt6) and clustering library (leidenalg) are
themselves GPL, so any distributed build is GPL whatever licence the
KOSMIC code carried. The licences and terms of the bundled reference
data are listed in [NOTICE.md](NOTICE.md); two gene-set files are
derived from KEGG, whose terms permit academic use only.

## Authors

- **Calum Wilson** — University of Strathclyde
- **Kirk Franks** — University of Strathclyde
