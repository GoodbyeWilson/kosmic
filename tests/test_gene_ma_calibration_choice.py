"""The Gene Meta-Analysis page offers the exact and the vectorised permutation
(ADR-018, ADR-019).

The Calibration box names the engine; the page turns that into the
permutation's DE engine and, for the engines that re-fit each study, a resume
folder beside the selection's results. The tests cover the choice, what Welch
input does with it, the file name that keeps the engines' runs apart, and that
a run in progress turns the run button into a Stop button.
"""
from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from kosmic.gui.meta_analysis.pages.gene_ma_page import GeneMAPage  # noqa: E402
from kosmic.meta_analysis.saved import _GENE_RUN  # noqa: E402


@pytest.fixture(scope="module")
def app():
    inst = QApplication.instance() or QApplication([])
    yield inst


@pytest.fixture
def page(app, tmp_path):
    p = GeneMAPage()
    p._project_folder = str(tmp_path)
    p._output_selection = 'decontX_counts/Endothelial_Cell'
    return p


def test_the_box_offers_analytical_fast_and_exact(page):
    items = [page._cal_combo.itemText(i) for i in range(page._cal_combo.count())]
    assert items == ["Analytical", "CC permutation (fast)",
                     "CC permutation (exact, slow)"]


@pytest.mark.parametrize('index,key,engine', [
    (0, 'analytical', 'fast'), (1, 'cc_perm', 'fast'), (2, 'cc_perm', 'exact')])
def test_the_box_gives_a_calibration_and_an_engine(page, index, key, engine):
    page._cal_combo.setCurrentIndex(index)
    assert page._calibration_key() == key
    if key == 'cc_perm':
        assert page._cc_engine() == engine


def test_fast_uses_the_vectorised_engine_with_its_own_resume_folder(page, tmp_path):
    page._cal_combo.setCurrentIndex(1)
    engine, folder = page._cc_engine_settings('cc_perm', 'deseq2')
    assert engine == 'deseq2_vec'
    assert folder.replace('\\', '/').endswith(
        'meta_analysis/decontX_counts/Endothelial_Cell/cc_permutation_vec')
    assert folder.startswith(str(tmp_path))


def test_fast_falls_back_to_the_approximation_when_the_optimiser_is_unusable(
        page, monkeypatch):
    from kosmic.meta_analysis import vec_lbfgsb
    monkeypatch.setattr(vec_lbfgsb, '_SETULB', None)
    monkeypatch.setattr(vec_lbfgsb, '_AVAILABLE', None)
    messages = []
    page.log_message.connect(messages.append)
    page._cal_combo.setCurrentIndex(1)
    assert page._cc_engine_settings('cc_perm', 'deseq2') == ('deseq2_fast', None)
    assert any('approximation' in m for m in messages)


def test_exact_uses_the_exact_engine_with_a_folder_beside_the_results(page, tmp_path):
    page._cal_combo.setCurrentIndex(2)
    engine, folder = page._cc_engine_settings('cc_perm', 'deseq2')
    assert engine == 'deseq2_exact'
    assert folder.replace('\\', '/').endswith(
        'meta_analysis/decontX_counts/Endothelial_Cell/cc_permutation')
    assert folder.startswith(str(tmp_path))


def test_welch_input_keeps_the_welch_engine_under_fast(page):
    page._cal_combo.setCurrentIndex(1)
    assert page._cc_engine_settings('cc_perm', 'welch_cpm') == ('welch_cpm', None)


def test_welch_input_keeps_the_welch_engine_under_exact(page):
    page._cal_combo.setCurrentIndex(2)
    messages = []
    page.log_message.connect(messages.append)
    assert page._cc_engine_settings('cc_perm', 'welch_cpm') == ('welch_cpm', None)
    assert any('Welch' in m for m in messages)


def test_analytical_calibration_never_selects_the_exact_engine(page):
    page._cal_combo.setCurrentIndex(0)
    assert page._cc_engine_settings('analytical', 'deseq2') == ('deseq2_fast', None)


def test_saved_runs_of_every_engine_are_listed(page):
    for stem, cal in (('consensus_reml_analytical.csv', 'analytical'),
                      ('consensus_reml_cc.csv', 'cc'),
                      ('consensus_reml_ccexact.csv', 'ccexact'),
                      ('consensus_reml_ccvec.csv', 'ccvec'),
                      ('consensus_reml_sumrank_ccvec.csv', 'ccvec'),
                      ('consensus_reml_sumrank_ccexact.csv', 'ccexact')):
        m = _GENE_RUN.match(stem)
        assert m and m.group(2) == cal


def test_an_exact_run_turns_the_run_button_into_stop(page):
    page._run_btn.setText("Run Meta-Analysis")
    page._after_worker_start({'calibration': 'cc_perm', 'de_method': 'deseq2_exact'})
    assert page._exact_run_active is True
    assert page._run_btn.isEnabled() and 'Stop' in page._run_btn.text()
    page._reset_run_button()
    assert page._exact_run_active is False
    assert page._run_btn.text() == "Run Meta-Analysis"


def test_a_vectorised_run_also_turns_the_run_button_into_stop(page):
    page._run_btn.setText("Run Meta-Analysis")
    page._after_worker_start({'calibration': 'cc_perm', 'de_method': 'deseq2_vec'})
    assert page._exact_run_active is True and 'Stop' in page._run_btn.text()
    page._reset_run_button()


def test_the_saved_file_name_tags_the_engine_that_made_the_run(page):
    page._cal_combo.setCurrentIndex(0)
    assert page._cal_tag() == 'analytical'
    page._cal_combo.setCurrentIndex(1)
    for engine, tag in (('deseq2_vec', 'ccvec'), ('deseq2_exact', 'ccexact'),
                        ('deseq2_fast', 'cc'), ('welch_cpm', 'cc')):
        page._last_consensus_params = {'de_method': engine}
        assert page._cal_tag() == tag


def test_a_fast_run_leaves_the_run_button_alone(page):
    page._run_btn.setText("Run Meta-Analysis")
    page._after_worker_start({'calibration': 'cc_perm', 'de_method': 'deseq2_fast'})
    assert page._exact_run_active is False
    assert page._run_btn.text() == "Run Meta-Analysis"
