"""Running a comparison end to end, and the Methods Comparison page (ADR-017)."""
import json
import os

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from kosmic.comparison.agreement import Rule  # noqa: E402
from kosmic.comparison.run import available_selections, check, run_comparison  # noqa: E402

CT, SLUG = 'Ventricular Cardiomyocyte', 'Ventricular_Cardiomyocyte'
STUDIES = ('S1', 'S2', 'S3')
ALL = ['pb_merged', 'ma_reml', 'ma_sumrank', 'ma_stouffer', 'ma_vote']


def _obs(study):
    rows = []
    for role in ('disease', 'control'):
        for i in range(3):
            rows += [(f'{study}_{role}{i}', role, ('female', 'male')[i % 2], CT)] * 20
    return pd.DataFrame(rows, columns=['sample', '_role', 'sex_inferred', 'cell_type'])


def _de_table(seed, n=60):
    rng = np.random.default_rng(seed)
    lfc = rng.normal(0, 0.05, n)
    lfc[:5] += 2.0
    p = np.clip(np.exp(-np.abs(lfc) * 8), 1e-12, 1)
    return pd.DataFrame({'names': [f'G{i}' for i in range(n)], 'logfoldchanges': lfc,
                         'se': 0.2, 'pvals': p, 'pvals_adj': np.minimum(p * 10, 1),
                         'pct_disease': 0.6, 'pct_control': 0.6})


@pytest.fixture
def project(tmp_path):
    rng = np.random.default_rng(0)
    atlas_obs = []
    for k, study in enumerate(STUDIES):
        obs = _obs(study)
        proc = tmp_path / study / 'processed_data'
        proc.mkdir(parents=True)
        a = anndata.AnnData(X=np.zeros((len(obs), 2), np.float32),
                            obs=obs.set_axis([f'c{i}' for i in range(len(obs))]))
        a.uns['sex_column'] = 'sex_inferred'
        a.write_h5ad(proc / f'{study}.h5ad')
        run = tmp_path / study / 'results/de_analysis/statistics/original_counts' / SLUG
        run.mkdir(parents=True)
        path = run / f'{study}_{SLUG}_DE_deseq2.csv'
        _de_table(k).to_csv(path, index=False)
        path.with_suffix('.json').write_text(json.dumps({
            'method': 'deseq2', 'cell_type': CT, 'cell_type_column': 'cell_type_atlas',
            'sample_column': 'sample', 'condition_column': '_role', 'counts_source': 'raw',
            'min_cells': 5, 'min_counts': 0, 'covariates': ['sex_inferred'],
            'dataset_file': f'{study}.h5ad'}), encoding='utf-8')
        atlas_obs.append(obs.rename(columns={'sex_inferred': 'sex'}).assign(study=study))
    obs = pd.concat(atlas_obs, ignore_index=True)
    obs.index = [f'm{i}' for i in range(len(obs))]
    X = rng.poisson(8, size=(len(obs), 60)).astype(np.float32)
    X[obs['_role'].values == 'disease', :5] *= 4
    master = tmp_path / '_master/processed_data'
    master.mkdir(parents=True)
    anndata.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f'G{i}' for i in range(60)])
                    ).write_h5ad(master / 'master.h5ad')
    return tmp_path


def test_selections_are_cell_types_run_in_two_or_more_studies(project):
    (sel,) = available_selections(project)
    assert sel.cell_type == CT and sel.counts == 'original_counts'
    assert sel.path == f'original_counts/{SLUG}' and sel.studies == list(STUDIES)
    assert sel.names[0] == f'S1_{SLUG}'
    assert sel.restricted(['S1', 'S3']).studies == ['S1', 'S3']


def test_a_first_run_fits_merged_pseudobulk_and_saves_the_comparison(project):
    (sel,) = available_selections(project)
    c = check(project, sel)
    assert c.design.formula == '~sex + study + condition' and c.saved is None
    b = run_comparison(project, sel, ALL, Rule())
    assert b.params['sources']['pb_merged']['source'] == 'fitted'
    assert {b.params['sources'][m]['source'] for m in ALL[1:]} == {'pooled'}
    assert set(b.methods) == set(ALL)
    assert (project / 'methods_comparison' / sel.path / 'results.csv').is_file()
    # The genes changed in every study are called by REML.
    assert {f'G{i}' for i in range(5)} <= b.sets['ma_reml']
    # A second run reads the merged run it saved.
    again = run_comparison(project, sel, ['pb_merged', 'ma_reml'], Rule())
    assert again.params['sources']['pb_merged']['source'] == 'saved'


def _break_saved_design(project, sel):
    c = check(project, sel)
    record = json.loads(c.saved.path.with_suffix('.json').read_text())
    record['covariates'] = ['sex']          # study left out
    c.saved.path.with_suffix('.json').write_text(json.dumps(record))
    return check(project, sel)


def test_a_differing_saved_run_needs_a_choice(project):
    (sel,) = available_selections(project)
    run_comparison(project, sel, ['pb_merged', 'ma_reml'], Rule())
    c = _break_saved_design(project, sel)
    assert c.needs_choice
    with pytest.raises(ValueError, match='choose'):
        run_comparison(project, sel, ['pb_merged', 'ma_reml'], Rule(), checked=c)
    anyway = run_comparison(project, sel, ['pb_merged', 'ma_reml'], Rule(),
                            merged_action='saved', checked=c)
    assert any('study' in d for d in anyway.params['design_differences'])
    resolved = run_comparison(project, sel, ['pb_merged', 'ma_reml'], Rule(),
                              merged_action='refit', checked=c)
    assert resolved.params['sources']['pb_merged']['source'] == 'fitted'
    assert check(project, sel).saved.matches


def test_one_method_is_not_a_comparison(project):
    (sel,) = available_selections(project)
    with pytest.raises(ValueError, match='two methods'):
        run_comparison(project, sel, ['ma_reml'], Rule())


# -- the page -------------------------------------------------------------------

@pytest.fixture(scope='module')
def app():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def workspace(app):
    from PyQt6.QtCore import QEvent
    from kosmic.gui.comparison.workspace import ComparisonWorkspace
    ws = ComparisonWorkspace()
    yield ws
    for w in list(ws.compare_page._workers):
        w.wait()
    ws.close()
    ws.deleteLater()
    app.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def _wait(app, page):
    for w in list(page._workers):
        w.wait()
    app.processEvents()


def test_the_page_lists_the_selection_and_checks_its_design(app, workspace, project):
    workspace.set_project_directory_external(project)
    page = workspace.compare_page
    assert page.selection_combo.count() == 1 and page.study_list.count() == 3
    _wait(app, page)
    assert page.checked is not None
    assert '~sex + study + condition' in page.design_label.text()
    assert page.run_btn.isEnabled()
    page.method_checks['ma_sumrank'].setChecked(False)
    for m in ('ma_reml', 'ma_stouffer', 'ma_vote'):
        page.method_checks[m].setChecked(False)
    assert not page.run_btn.isEnabled()          # one method left


def test_running_shows_the_summary_and_a_rule_change_recuts_it(app, workspace, project):
    workspace.set_project_directory_external(project)
    page = workspace.compare_page
    _wait(app, page)
    page._on_run()
    _wait(app, page)
    assert workspace.bundle is not None and page.summary_table.rowCount() == 5
    assert page.pairwise_table.rowCount() == 10
    page.rule_combo.setCurrentIndex(page.rule_combo.findData('top_n'))
    page.top_spin.setValue(3)
    assert workspace.bundle.rule.rule == 'top_n'
    assert all(len(s) == 3 for s in workspace.bundle.sets.values())
    # The saved comparison follows the rule shown, and is shown again on return.
    from kosmic.comparison.bundle import load
    assert load(project / 'methods_comparison' / f'original_counts/{SLUG}').rule.top_n == 3


def test_a_differing_saved_run_asks_and_cancel_runs_nothing(app, workspace, project,
                                                           monkeypatch):
    (sel,) = available_selections(project)
    run_comparison(project, sel, ['pb_merged', 'ma_reml'], Rule())
    _break_saved_design(project, sel)
    workspace.set_project_directory_external(project)
    page = workspace.compare_page
    _wait(app, page)
    asked = {}

    def choose(parent, title, text, options, escape=None):
        asked['options'] = options
        asked['text'] = text
        return asked.get('answer', 'Cancel')

    started = []
    monkeypatch.setattr('kosmic.gui.shared.dialogs.choose', choose)
    monkeypatch.setattr('kosmic.gui.comparison.pages.compare_page._comparison_worker',
                        lambda *a: started.append(a) or _NoWorker())
    page._on_run()
    assert asked['options'] == ['Resolve', 'Run anyway', 'Cancel'] and not started
    assert 'does not include study' in asked['text']
    asked['answer'] = 'Run anyway'
    page._on_run()
    assert started[-1][4] == 'saved'
    asked['answer'] = 'Resolve'
    page._on_run()
    assert started[-1][4] == 'refit'


class _NoWorker:
    """Stands in for a worker in the dialog test: nothing runs."""

    def __init__(self):
        from PyQt6.QtCore import QObject, pyqtSignal

        class _S(QObject):
            sig = pyqtSignal(object)
        self._s = _S()
        self.progress = self.progress_pct = self.finished_ok = self.failed = \
            self.finished = self._s.sig

    def start(self):
        pass

    def isRunning(self):
        return False

    def wait(self):
        return True
