"""A study analysed per cell type has no whole-study DE file. Reopening it
must load the per-cell-type files, or the volcano (DE page and Figures)
has no results to show."""
from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from kosmic.gui.de_analysis.pages.gene_de_page import GeneDEPage


def _write(path, n_sig):
    pd.DataFrame({
        'names': ['A', 'B', 'C'],
        'logfoldchanges': [1.0, -1.0, 0.0],
        'abs_logfoldchange': [1.0, 1.0, 0.0],
        'pvals_adj': [0.01 if i < n_sig else 0.5 for i in range(3)],
        'disease_samples': [4, 4, 4],
        'control_samples': [5, 5, 5],
    }).to_csv(path, index=False)


def test_per_cell_type_files_are_loaded(tmp_path):
    _write(tmp_path / 'GSE1_Fibroblast_DE_deseq2.csv', n_sig=2)
    _write(tmp_path / 'GSE1_Ventricular_Cardiomyocyte_DE_deseq2.csv', n_sig=1)
    _write(tmp_path / 'GSE1_Fibroblast_DE_welch_cpm.csv', n_sig=0)   # other method
    page = SimpleNamespace(
        pval_filter=SimpleNamespace(value=lambda: 0.05),
        fc_filter=SimpleNamespace(value=lambda: 0.25),
        ws=SimpleNamespace(log_message=SimpleNamespace(emit=print)),
    )

    runs = GeneDEPage._load_cell_type_runs(page, tmp_path, 'GSE1', 'deseq2')

    assert [r.cell_type for r in runs] == ['Fibroblast', 'Ventricular_Cardiomyocyte']
    assert [r.n_significant for r in runs] == [2, 1]
    assert all(r.n_samples == 9 and r.n_genes_tested == 3 for r in runs)
    assert all(r.de_results is not None for r in runs)
