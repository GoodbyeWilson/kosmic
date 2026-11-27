"""Preranked GSEA: the DE page's 'Rank genes by' choice, and the meta-analysis
Enrichment step saving each library GSEA under the pooling run it ranked, which
Figure Export then reads with that run (kosmic.meta_analysis.saved)."""
import os
import time

import pandas as pd
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


DE = pd.DataFrame({'names': ['A', 'B', 'C'], 'logfoldchanges': [3.0, -0.5, 0.2],
                   'statistic': [1.0, -6.0, 4.0], 'pvals_adj': [0.2, 0.01, 0.03],
                   'pvals': [0.1, 0.001, 0.003]})


class _Ws:
    """The parts of the DE workspace the Enrichment page reads for GSEA."""

    def __init__(self, de):
        from PyQt6.QtCore import QObject, pyqtSignal

        class _Sig(QObject):
            log_message = pyqtSignal(str)
            status_message = pyqtSignal(str)
        self._sig = _Sig()
        self.log_message = self._sig.log_message
        self.status_message = self._sig.status_message
        self.de_results = de
        self.gene_de_genome_wide = de
        self.progress_bar = None
        self.gsea_library_name = None
        self.project_dir = None


def _de_page(de):
    from kosmic.gui.de_analysis.pages.enrichment_page import EnrichmentPage
    return EnrichmentPage(_Ws(de))


def _capture(monkeypatch):
    """Stop the GSEA worker from running; record the ranking it was given."""
    import kosmic.gui.de_analysis.pages.enrichment_page as ep
    import kosmic.gui.de_analysis.workers as workers
    seen = {}

    class FakeWorker:
        def __init__(self, rank_df, lib_id, *a, **k):
            seen['rank'] = rank_df.set_index('names')['logfoldchanges'].to_dict()
            seen['library'] = lib_id
    monkeypatch.setattr(workers, 'LibraryGseaWorker', FakeWorker)
    monkeypatch.setattr(ep, 'run_worker', lambda *a, **k: None)
    return seen


def test_the_de_page_runs_gsea_of_hallmark_ranked_by_the_statistic(app):
    page = _de_page(DE)
    assert page._library_combo.currentData() == 'MSigDB_Hallmark_2026.2_Hs'
    assert page._rank_by() == 'statistic'
    assert not hasattr(page, '_method_combo') and not hasattr(page, '_rank_combo')
    assert 'Wald statistic' in page._about.text()


def test_gsea_ranks_by_the_test_statistic(app, monkeypatch):
    seen = _capture(monkeypatch)
    page = _de_page(DE)
    page._run_gsea()
    assert seen['rank'] == {'A': 1.0, 'B': -6.0, 'C': 4.0}
    params = page._enrichment_params()
    assert params['method'] == 'gsea' and params['rank_by'] == 'statistic'
    assert (params['gsea_min_size'], params['gsea_max_size']) == (15, 500)


def test_a_result_without_the_statistic_is_not_ranked_by_something_else(app, monkeypatch):
    seen = _capture(monkeypatch)
    page = _de_page(DE.drop(columns='statistic'))
    page._run_gsea()
    assert 'rank' not in seen
    assert "no 'statistic' column" in page._status_label.text()
    assert page._run_btn.isEnabled()


def test_the_meta_step_ranks_by_the_pooled_z(app):
    from kosmic.gui.meta_analysis.pages.enrichment_step import MetaEnrichmentPage
    p = MetaEnrichmentPage()
    assert p._rank_by() == 'z_stat'
    assert 'pooled z-statistic' in p._about.text()


def test_the_meta_step_saves_gsea_under_the_run_it_ranked(app, tmp_path):
    from kosmic.gui.meta_analysis.pages.enrichment_step import MetaEnrichmentPage
    from kosmic.paths import meta_output_dir

    p = MetaEnrichmentPage()
    p.set_project_folder(tmp_path, selection='decontX_counts/Cardiomyocyte')
    p._gsea_library_id = 'MSigDB_Hallmark_2026.2_Hs'
    folder = meta_output_dir(tmp_path, 'decontX_counts/Cardiomyocyte')
    p.set_meta_results(DE.assign(z_stat=DE['statistic']), run_stem='consensus_reml_analytical')
    assert p._gsea_results_path() == folder / 'consensus_reml_analytical_gsea_MSigDB_Hallmark_2026.2_Hs.csv'
    p.set_meta_results(DE.assign(z_stat=DE['statistic']), run_stem='consensus_stouffer_analytical')
    assert p._gsea_results_path().name == 'consensus_stouffer_analytical_gsea_MSigDB_Hallmark_2026.2_Hs.csv'
    p.set_meta_results(DE.assign(z_stat=DE['statistic']))       # run not known: the old name
    assert p._gsea_results_path().name == 'gsea_library_results_MSigDB_Hallmark_2026.2_Hs.csv'


def test_the_run_context_names_the_saved_run(app):
    from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage
    page = GeneMAPage()
    assert page.get_run_context()['run_stem'] is None          # nothing run yet
    page._meta_df = pd.DataFrame({'names': ['A']})
    page._method_keys = ['reml']
    assert page.get_run_context()['run_stem'] == 'consensus_reml_analytical'


def test_figure_export_reads_each_runs_own_gsea(tmp_path):
    from kosmic.meta_analysis.saved import MetaRun, figure_inputs
    from kosmic.paths import meta_output_dir

    folder = meta_output_dir(tmp_path, 'decontX_counts/Cardiomyocyte')
    folder.mkdir(parents=True)
    runs = {}
    for k in ('reml', 'stouffer'):
        runs[k] = folder / f'consensus_{k}_analytical.csv'
        runs[k].write_text('names\nA\n')
        pd.DataFrame({'names': [f'SET_{k.upper()}'], 'nes': [1.0], 'pvals_adj': [0.01]}).to_csv(
            folder / f'consensus_{k}_analytical_gsea_MSigDB_Hallmark_2026.2_Hs.csv', index=False)
        time.sleep(0.02)
    # An older-style file, newer than both: it must not be taken for either run.
    pd.DataFrame({'names': ['SET_OLD'], 'nes': [1.0], 'pvals_adj': [0.01]}).to_csv(
        folder / 'gsea_library_results_MSigDB_Hallmark_2026.2_Hs.csv', index=False)
    for k in ('reml', 'stouffer'):
        fi = figure_inputs(MetaRun(tmp_path, 'decontX_counts/Cardiomyocyte', runs[k], None))
        assert list(fi['gsea_library_results']['names']) == [f'SET_{k.upper()}']
    # A run with no GSEA of its own falls back to the older-style file.
    other = folder / 'consensus_sumrank_analytical.csv'
    other.write_text('names\nA\n')
    fi = figure_inputs(MetaRun(tmp_path, 'decontX_counts/Cardiomyocyte', other, None))
    assert list(fi['gsea_library_results']['names']) == ['SET_OLD']
