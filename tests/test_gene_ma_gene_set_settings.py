"""The Gene Meta-Analysis page's gene-set settings (current symbols only,
inconsistently merged genes left out) and permutation settings (null from the
pooled gene set, no CPM refilter on shuffles). All off by default."""
from __future__ import annotations

import os

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage  # noqa: E402


@pytest.fixture(scope="module")
def app():
    yield QApplication.instance() or QApplication([])


@pytest.fixture
def page(app, tmp_path):
    p = GeneMAPage()
    p._project_folder = str(tmp_path)
    p._output_selection = 'decontX_counts/Endothelial_Cell'
    return p


def _datasets():
    names = ["TTN", "NPPA", "ENSG00000000001", "C11orf58"]
    return [{"name": f"S{i}", "df": pd.DataFrame({
        "names": names, "pct_disease": [.5] * 4, "pct_control": [.5] * 4})}
        for i in (1, 2)]


def test_every_setting_is_off_by_default(page):
    page._if_check.setChecked(False)
    assert not page._symbols_check.isChecked()
    assert not page._merge_check.isChecked()
    assert not page._cc_universe_check.isChecked()
    assert not page._cc_refilter_check.isChecked()
    data = _datasets()
    assert page._apply_gene_set_rules(data) is data       # pooled as it is
    assert page._last_gene_set is None


def test_the_permutation_settings_show_only_for_a_permutation(page):
    page._cal_combo.setCurrentIndex(0)
    assert page._cc_universe_check.isHidden() and page._cc_refilter_check.isHidden()
    page._cal_combo.setCurrentIndex(1)
    assert not page._cc_universe_check.isHidden()


def test_current_symbols_only_drops_names_that_are_not_symbols(page):
    page._if_check.setChecked(False)
    page._symbols_check.setChecked(True)
    out = page._apply_gene_set_rules(_datasets())
    assert list(out[0]["df"]["names"]) == ["TTN", "NPPA", "C11orf58"]
    assert page._last_gene_set_record["rules"][0]["rule"] == "current symbols only"
    assert page._last_gene_set_record["rules"][0]["n_removed"] == 1


def test_the_permutation_gets_the_gene_set_and_the_refilter_choice(page):
    page._if_check.setChecked(False)
    page._symbols_check.setChecked(True)
    page._apply_gene_set_rules(_datasets())
    assert page._cc_gene_set_params('analytical', 'deseq2_vec') == {
        'cc_gene_universe': None, 'cc_refilter': True}
    assert page._cc_gene_set_params('cc_perm', 'deseq2_vec') == {
        'cc_gene_universe': None, 'cc_refilter': True}
    page._cc_universe_check.setChecked(True)
    page._cc_refilter_check.setChecked(True)
    assert page._cc_gene_set_params('cc_perm', 'deseq2_exact') == {
        'cc_gene_universe': ["C11orf58", "NPPA", "TTN"], 'cc_refilter': False}


def test_an_engine_that_does_not_refit_keeps_refiltering(page):
    page._cc_refilter_check.setChecked(True)
    assert page._cc_gene_set_params('cc_perm', 'deseq2_fast')['cc_refilter'] is True


def test_inconsistently_merged_genes_are_left_out_and_unchecked_studies_named(page, tmp_path):
    import anndata as ad
    import numpy as np
    from scipy.sparse import csr_matrix
    from kosmic.provenance import record_stage
    for accession, merged in (("S1", {"NPPA": ["NPPA", "NPPA-AS1"]}), ("S2", {})):
        processed = tmp_path / accession / "processed_data"
        processed.mkdir(parents=True)
        ad.AnnData(X=csr_matrix(np.ones((2, 2), dtype=np.float32)),
                   var=pd.DataFrame(index=["TTN", "NPPA"]),
                   obs=pd.DataFrame(index=["c0", "c1"])).write_h5ad(processed / f"{accession}.h5ad")
        record_stage(processed, accession, "gene_names", params={"merged": merged})
    data = _datasets()
    for d in data:
        d["path"] = str(tmp_path / d["name"] / "statistics" / "de.csv")
    data.append({**_datasets()[0], "name": "imported", "external": True})
    page._if_check.setChecked(False)
    page._symbols_check.setChecked(False)
    page._merge_check.setChecked(True)
    out = page._apply_gene_set_rules(data)
    assert "NPPA" not in set(out[0]["df"]["names"])
    assert page._last_gene_set_record["inconsistent_merge_unchecked"] == ["imported"]
