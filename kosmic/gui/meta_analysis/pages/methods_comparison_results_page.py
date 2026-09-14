# Consensus Evaluation Page (Methods Comparison, step 2 of 2).
#
# Pareto-frontier combination explorer. Given the per-method results
# from Individual Methods, builds and compares arbitrary method
# combinations (intersection / union / k-of-N voting for k in {3, 5, 7})
# and highlights the Pareto frontier across the selected x/y metric
# pair on the scatter plot.
#
# LOO direction replication and AUC are filled in for single methods
# from the Cross-Study Validation cache; LOO for combinations is blank.

from __future__ import annotations

from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import pyqtgraph as pg
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QAbstractItemView, QComboBox, QFrame, QGridLayout,
    QHBoxLayout, QHeaderView, QLabel, QTabWidget, QVBoxLayout, QWidget,
)

from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets import Column, ResultsTable, SecondaryLabel, SimplePage
from kosmic.gui.meta_analysis.pages.gene_ma_page import _METHODS as METHODS_REGISTRY
from kosmic.meta_analysis.consensus_loo import _replication_test
from kosmic import DEFAULT_FDR


# Combination computation (pure, no GUI)

def _gene_set_at_fdr(df: pd.DataFrame, fdr: float = DEFAULT_FDR) -> Set[str]:
    if df is None or len(df) == 0 or 'fdr' not in df.columns:
        return set()
    return set(df.loc[df['fdr'] < fdr, 'names'].astype(str).tolist())


def _combine_gene_sets(method_sets: List[Set[str]], rule: str,
                       k: int = 1) -> Set[str]:
    if not method_sets:
        return set()
    if rule == 'intersection':
        out = set(method_sets[0])
        for s in method_sets[1:]:
            out &= s
        return out
    if rule == 'union':
        out = set()
        for s in method_sets:
            out |= s
        return out
    if rule == 'voting':
        # Gene must be in >= k of the method sets
        from collections import Counter
        c: Counter = Counter()
        for s in method_sets:
            c.update(s)
        return {g for g, n in c.items() if n >= k}
    raise ValueError(f"unknown combination rule: {rule}")


def _compute_combo_metrics(method_dfs: Dict[str, pd.DataFrame],
                           selected: List[str],
                           rule: str,
                           k: int,
                           fdr: float) -> dict:
    """
    Compute full-K cheap metrics for a combination from cached per-
    method DataFrames. Returns: n_called, n_meta, n_obvious, pct_meta,
    mean_abs_lfc, mean_concord, pct_unanimous.
    """
    sets = [_gene_set_at_fdr(method_dfs.get(m), fdr) for m in selected]
    combined = _combine_gene_sets(sets, rule, k)
    n_called = len(combined)
    if n_called == 0:
        return {
            'n_called': 0, 'n_meta': 0, 'n_obvious': 0, 'pct_meta': 0.0,
            'mean_abs_lfc': 0.0, 'mean_concord': 0.0,
            'pct_unanimous': 0.0,
            'gene_set': set(),
        }
    # Pull per-gene metadata from the first selected method's DF
    # (n_contributing_sig is gene-property, same across methods)
    base = next((method_dfs[m] for m in selected
                 if method_dfs.get(m) is not None
                 and len(method_dfs[m]) > 0), None)
    if base is None:
        return {
            'n_called': n_called, 'n_meta': 0, 'n_obvious': 0,
            'pct_meta': 0.0, 'mean_abs_lfc': 0.0, 'mean_concord': 0.0,
            'pct_unanimous': 0.0, 'gene_set': combined,
        }
    sig = base[base['names'].isin(combined)].copy()
    n_meta = int((sig['n_contributing_sig'] == 0).sum()) \
        if 'n_contributing_sig' in sig.columns else 0
    n_obv = n_called - n_meta
    pct_meta = 100.0 * n_meta / n_called if n_called > 0 else 0.0
    mean_abs_lfc = float(np.abs(sig['logfoldchanges']).mean()) \
        if 'logfoldchanges' in sig.columns and len(sig) > 0 else 0.0
    if 'frac_agreeing' in sig.columns and len(sig) > 0:
        frac = sig['frac_agreeing'].to_numpy(dtype=float)
        finite = np.isfinite(frac)
        mean_concord = (100.0 * float(frac[finite].mean())
                        if finite.sum() > 0 else 0.0)
        pct_unanimous = (100.0 * float(np.isclose(frac[finite], 1.0).sum())
                         / n_called)
    else:
        mean_concord = 0.0
        pct_unanimous = 0.0
    return {
        'n_called': n_called,
        'n_meta': n_meta,
        'n_obvious': n_obv,
        'pct_meta': pct_meta,
        'mean_abs_lfc': mean_abs_lfc,
        'mean_concord': mean_concord,
        'pct_unanimous': pct_unanimous,
        'gene_set': combined,
    }


# Pareto frontier

def _pareto_frontier(points: List[Tuple[float, float]],
                     maximize_x: bool, maximize_y: bool) -> List[int]:
    """
    Return indices of points on the Pareto frontier.

    A point (x, y) is on the frontier if no other point dominates it
    on both axes. Domination direction is set by 'maximize_x' /
    'maximize_y': True means "higher is better" on that axis.
    """
    n = len(points)
    if n == 0:
        return []
    pts = np.asarray(points, dtype=float)
    sx = 1.0 if maximize_x else -1.0
    sy = 1.0 if maximize_y else -1.0
    pts_oriented = pts * np.array([sx, sy])

    on_frontier = []
    for i in range(n):
        if not np.isfinite(pts_oriented[i]).all():
            continue
        dominated = False
        for j in range(n):
            if i == j:
                continue
            if not np.isfinite(pts_oriented[j]).all():
                continue
            # j dominates i if j >= i on both axes and strictly > on one
            if (pts_oriented[j, 0] >= pts_oriented[i, 0]
                    and pts_oriented[j, 1] >= pts_oriented[i, 1]
                    and (pts_oriented[j, 0] > pts_oriented[i, 0]
                         or pts_oriented[j, 1] > pts_oriented[i, 1])):
                dominated = True
                break
        if not dominated:
            on_frontier.append(i)
    return on_frontier


def _combination_direction_loo(per_method_loo: dict,
                                methods: List[str],
                                rule: str,
                                k: int,
                                replication_threshold: float = DEFAULT_FDR,
                                ) -> dict:
    """
    Direction-replication LOO for a method combination, from cached
    per-fold stashes written by 'run_consensus_loo(stash_fold_details=True)'.

    Returns a dict with:
        overall_dir_pct  -- mean % of combination genes whose logFC sign
                           matches the held-out study, averaged over folds
        meta_dir_pct     -- same but restricted to meta-discovered genes
                           (n_contributing_sig == 0 in the K-1 fold)
        wrong_dir_pct    -- mean % with opposite sign (100 - overall_dir,
                           after excluding missing)
        n_folds          -- number of folds where the combination had
                           >= 1 testable gene (non-empty, non-all-missing)
    Any field is 'float('nan')' if no fold had testable genes.
    """
    nan = float('nan')
    missing = {'overall_dir_pct': nan, 'meta_dir_pct': nan,
               'wrong_dir_pct': nan, 'n_folds': 0}

    if not per_method_loo or not methods:
        return missing
    # Validate every requested method has a stash
    for m in methods:
        res = per_method_loo.get(m)
        if res is None or 'folds' not in res:
            return missing
        if not any('_stash' in f for f in res['folds']):
            return missing

    # Align folds by fold_index across methods (they all iterate the
    # same datasets in the same order, but match by index to be safe).
    n_folds_expected = min(len(per_method_loo[m]['folds']) for m in methods)

    overall_vals, meta_vals, wrong_vals = [], [], []
    for j in range(n_folds_expected):
        # Collect each method's called-gene set for this fold
        sets = []
        for m in methods:
            f = per_method_loo[m]['folds'][j]
            stash = f.get('_stash')
            if stash is None:
                sets = None
                break
            cg = stash['called_genes_by_method'].get(m)
            if cg is None:
                sets = None
                break
            sets.append(set(cg))
        if not sets:
            continue

        # Apply the combination rule
        if rule == 'intersection':
            combo = set(sets[0])
            for s in sets[1:]:
                combo &= s
        elif rule == 'union':
            combo = set()
            for s in sets:
                combo |= s
        elif rule == 'voting':
            from collections import Counter
            c = Counter()
            for s in sets:
                c.update(s)
            combo = {g for g, n in c.items() if n >= k}
        else:
            continue

        if not combo:
            continue

        # Use the first method's fold cache for the reference
        # logfoldchanges + held-out + n_sig_per_gene (they're identical
        # across methods for the same fold).
        ref_fold = per_method_loo[methods[0]]['folds'][j]
        ref_stash = ref_fold['_stash']
        all_genes = ref_stash['all_genes_with_lfc']
        held_out = ref_stash['held_out_min']
        n_sig_per_gene = ref_stash['n_sig_per_gene']

        combo_df = all_genes[all_genes['names'].isin(combo)].copy()
        if len(combo_df) == 0:
            continue

        rep = _replication_test(
            combo_df, held_out, mode='direction_only',
            fdr_threshold=replication_threshold)

        n_testable = int((rep['reason'] != 'missing').sum())
        if n_testable == 0:
            continue
        n_rep = int((rep['reason'] == 'replicated').sum())
        n_wrong = int((rep['reason'] == 'wrong_direction').sum())
        overall = 100.0 * n_rep / n_testable
        wrong = 100.0 * n_wrong / n_testable
        overall_vals.append(overall)
        wrong_vals.append(wrong)

        # Meta-discovered subset: genes with n_contributing_sig == 0
        meta_genes = [g for g in combo if n_sig_per_gene.get(g, 0) == 0]
        if meta_genes:
            meta_df = all_genes[all_genes['names'].isin(meta_genes)].copy()
            if len(meta_df) > 0:
                meta_rep = _replication_test(
                    meta_df, held_out, mode='direction_only',
                    fdr_threshold=replication_threshold)
                m_test = int((meta_rep['reason'] != 'missing').sum())
                if m_test > 0:
                    m_rep = int((meta_rep['reason'] == 'replicated').sum())
                    meta_vals.append(100.0 * m_rep / m_test)

    import numpy as _np
    def _mean(vs):
        return float(_np.mean(vs)) if vs else nan
    return {
        'overall_dir_pct': _mean(overall_vals),
        'meta_dir_pct':    _mean(meta_vals),
        'wrong_dir_pct':   _mean(wrong_vals),
        'n_folds':         len(overall_vals),
    }


def _robust_frontier_scores(rows: List[dict]) -> List[int]:
    """
    For each row, count how many of the C(n_metrics, 2) axis pairs it
    is non-dominated in. High = robustly good across many trade-offs.
    """
    n = len(rows)
    if n == 0:
        return []
    # Pre-extract metric arrays (NaN -> NaN)
    metric_keys = [k for k, _, _ in _METRICS]
    metric_max = [mx for _, _, mx in _METRICS]
    M = np.full((n, len(metric_keys)), np.nan, dtype=float)
    for i, r in enumerate(rows):
        for j, k in enumerate(metric_keys):
            v = r.get(k, float('nan'))
            try:
                M[i, j] = float(v)
            except (TypeError, ValueError):
                M[i, j] = float('nan')
    scores = np.zeros(n, dtype=int)
    for a in range(len(metric_keys)):
        for b in range(a + 1, len(metric_keys)):
            pts = list(zip(M[:, a].tolist(), M[:, b].tolist()))
            frontier = _pareto_frontier(pts, metric_max[a], metric_max[b])
            for i in frontier:
                scores[i] += 1
    return scores.tolist()


# Page

# Metric definitions: (column key, display label, maximize?)
_METRICS = [
    ('n_called',     'n called',         True),
    ('n_meta',       'n meta-discovered', True),
    ('pct_meta',     '% meta',           True),
    ('mean_abs_lfc', 'mean |logFC|',     True),
    ('mean_concord', 'mean concord %',   True),
    ('pct_unanimous', '% unanimous',     True),
    ('meta_dir_pct', 'meta dir % (LOO)', True),
    ('overall_dir_pct', 'overall dir % (LOO)', True),
    ('auc_z',        'AUC z (LOO)',      True),
    ('auc_u',        'AUC UCell (LOO)',  True),
    ('wrong_dir_pct', 'wrong-dir % (LOO)', False),
]
_METRIC_BY_KEY = {k: (lbl, mx) for k, lbl, mx in _METRICS}

# Fixed small-multiples panels: the 4 trade-offs that actually matter
# for method selection. (x_key, y_key, panel title)
_PANELS = [
    ('n_called', 'overall_dir_pct', 'Volume vs replication'),
    ('n_called', 'wrong_dir_pct',   'Volume vs wrongness'),
    ('n_meta',   'mean_concord',    'Discovery vs concord'),
    ('n_called', 'auc_z',           'Yield vs predictive AUC'),
]

# Unambiguous short codes for use in combination names. Falls back to
# the first word of the label from METHODS_REGISTRY if a key is missing
# here (but every key in the registry should have an entry).
_METHOD_SHORT = {
    'dl': 'DL',
    'reml': 'REML',
    'stouffer': 'Stouffer',
    'fisher': 'Fisher',
    'sumrank': 'SumRank',
    'gwop': 'gwOP',
}


class MethodsComparisonResultsPage(SimplePage):
    """Methods Comparison -- Pareto-frontier combination explorer."""

    help_id = "meta/methods_comparison_results"
    SCROLLABLE = True
    CONTENT_MARGINS = (12, 12, 12, 12)
    CONTENT_SPACING = 8

    log_message = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._method_results = None  # from Individual Methods Discovery
        self._loo_results = None     # from Individual Methods Cross-Study Validation
        self._project_folder = None
        self._rows: List[dict] = []  # rows in the comparison table
        self.progress_bar = None
        self._setup_ui()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _setup_ui(self):
        # SimplePage(SCROLLABLE=True) provides the QScrollArea + body
        # widget that caps minimumSizeHint regardless of inner table
        # width. Without that cap, a wide internal layout would propagate
        # up through the workspace QStackedLayout and force the main
        # window's sidebar splitter to collapse the workflow sidebar
        # when this page is active.
        outer = self.body_layout

        title = QLabel("Consensus Evaluation -- Trade-off Frontier")
        title.setProperty("role", "page_title")
        outer.addWidget(title)

        subtitle = SecondaryLabel(
            "Every single method plus every pairwise intersection and "
            "union, plus k-of-N voting combinations, scored on the same "
            "metrics. A method or combination is non-dominated (on the "
            "trade-off frontier) if no other row beats it on BOTH "
            "selected axes. Everything off the frontier is strictly "
            "worse than something on it.")
        subtitle.setWordWrap(True)
        outer.addWidget(subtitle)

        self._source_label = SecondaryLabel(
            "Run Individual Methods first. Cross-Study Validation also "
            "recommended (enables LOO-based metrics: direction "
            "replication and AUC).")
        self._source_label.setWordWrap(True)
        outer.addWidget(self._source_label)

        # FDR threshold: single top-level control. Applies to every row.
        ctrl_row = QHBoxLayout()
        ctrl_row.addWidget(QLabel("FDR threshold:"))
        self._fdr_combo = QComboBox()
        self._fdr_combo.addItems(["0.01", "0.05", "0.1", "0.25"])
        self._fdr_combo.setCurrentText("0.05")
        self._fdr_combo.currentIndexChanged.connect(self._rebuild_table)
        ctrl_row.addWidget(self._fdr_combo)
        ctrl_row.addStretch()
        outer.addLayout(ctrl_row)

        # Tabs: Comparison table / Trade-off frontier scatter
        self._tabs = QTabWidget()
        outer.addWidget(self._tabs, 1)

        # ---- Tab 1: Comparison table ------------------------------------
        table_tab = QWidget()
        tb_lay = QVBoxLayout(table_tab)
        tb_lay.setContentsMargins(6, 6, 6, 6)

        # Build schema: fixed cols + per-metric cols.  Metric fmt is
        # derived from the key (n_* → 'd', *_pct / mean_concord /
        # pct_unanimous → '.1%', auc_* → '.3f', others → '.3f').
        def _metric_fmt(key):
            if key in ('n_called', 'n_meta'):
                return 'd'
            if key.endswith('_pct') or key in ('mean_concord', 'pct_unanimous'):
                return '.1%'
            return '.3f'

        schema = [
            Column('Name',         'name',         's'),
            Column('Type',         'type',         's'),
            Column('Robust score', 'robust_score', 'd'),
        ]
        for key, lbl, _ in _METRICS:
            schema.append(Column(lbl, key, _metric_fmt(key), na_text='-'))

        self._table = ResultsTable()
        self._table.set_schema(schema)
        self._table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Interactive)
        self._table.horizontalHeader().setDefaultSectionSize(110)
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setSortingEnabled(True)
        # (Copy is inherited from CopyableTableWidget via ResultsTable.)
        tb_lay.addWidget(self._table)

        self._tabs.addTab(table_tab, "Comparison table")

        # ---- Tab 2: Trade-off frontier (fixed small-multiples) ----------
        # Four decision-relevant panels instead of a user-chosen axis
        # picker (which would force the user to pick among 55 pairings).
        scatter_tab = QWidget()
        pb_lay = QVBoxLayout(scatter_tab)
        pb_lay.setContentsMargins(6, 6, 6, 6)

        hint = SecondaryLabel(
            "Four fixed trade-offs. Blue points + line = non-dominated "
            "on that panel's axes; grey = dominated. See the Comparison "
            "table's 'Robust score' column for an overall ranking "
            "across all axis pairings.")
        hint.setWordWrap(True)
        pb_lay.addWidget(hint)

        grid = QGridLayout()
        grid.setSpacing(6)
        from kosmic.gui.shared.plots import InteractivePlot
        self._panel_plots: List[Tuple[pg.PlotWidget, QLabel]] = []
        for idx, (xk, yk, title) in enumerate(_PANELS):
            row, col = divmod(idx, 2)
            frame = QFrame()
            frame.setFrameShape(QFrame.Shape.StyledPanel)
            f_lay = QVBoxLayout(frame)
            f_lay.setContentsMargins(4, 4, 4, 4)
            lbl = QLabel(f"<b>{title}</b>")
            f_lay.addWidget(lbl)
            plot = InteractivePlot(
                title='', bottom_label=xk, left_label=yk,
                unavailable_message='Run Individual Methods to populate the trade-off frontier.',
            )
            f_lay.addWidget(plot, 1)
            cap = SecondaryLabel("")
            cap.setWordWrap(True)
            f_lay.addWidget(cap)
            grid.addWidget(frame, row, col)
            self._panel_plots.append((plot, cap))
        pb_lay.addLayout(grid, 1)

        self._tabs.addTab(scatter_tab, "Trade-off frontier")

    # ------------------------------------------------------------------
    # External activation
    # ------------------------------------------------------------------

    def set_method_results(self, method_results, project_folder=None,
                           loo_results=None):
        """
        Receive cached per-method results from Individual Methods,
        plus optional per-method LOO results from Cross-Study Validation.
        """
        self._method_results = method_results
        if project_folder:
            self._project_folder = project_folder
        if loo_results is not None:
            self._loo_results = loo_results
        else:
            # Fallback: try to pull from Individual Methods page if
            # workspace didn't pass it explicitly
            try:
                ind_page = self._get_individual_methods_page()
                if ind_page is not None and hasattr(ind_page, 'get_loo_results'):
                    self._loo_results = ind_page.get_loo_results()
            except Exception:
                pass
        self._rebuild_table()

    def _get_individual_methods_page(self):
        w = self.parent()
        while w is not None:
            if hasattr(w, '_methods_comparison_page'):
                return w._methods_comparison_page
            w = w.parent() if hasattr(w, 'parent') else None
        return None

    # ------------------------------------------------------------------
    # Table build (auto-enumerate full menu of methods + combinations)
    # ------------------------------------------------------------------

    def _current_fdr(self) -> float:
        try:
            return float(self._fdr_combo.currentText())
        except ValueError:
            return DEFAULT_FDR

    def _rebuild_table(self):
        """
        Rebuild the comparison table from scratch using the cached
        per-method results + LOO results. Auto-enumerates the full menu
        of combinations (singles, pairwise intersections, pairwise
        unions, voting k-of-N for k in {3, 5, 7}).
        """
        if self._method_results is None:
            self._source_label.setText(
                "Run Individual Methods first. Cross-Study Validation "
                "recommended for LOO metrics.")
            self._table.setRowCount(0)
            self._rows = []
            self._refresh_pareto()
            return

        keys = self._method_results.get('method_keys', [])
        method_dfs = self._method_results.get('method_dfs', {})
        key_to_label = dict(METHODS_REGISTRY)
        fdr = self._current_fdr()

        # Pull LOO data from cache (or fall back to peeking at the
        # Individual Methods page directly).
        loo_per_method = {}
        if self._loo_results is not None and 'per_method' in self._loo_results:
            loo_per_method = self._loo_results['per_method']
        else:
            ind_page = self._get_individual_methods_page()
            if ind_page is not None and hasattr(ind_page, 'get_loo_results'):
                loo_data = ind_page.get_loo_results()
                if loo_data is not None and 'per_method' in loo_data:
                    loo_per_method = loo_data['per_method']
                    self._loo_results = loo_data

        self._rows = []

        # ---- Singles ---------------------------------------------------
        valid_keys: List[str] = []
        for k in keys:
            df = method_dfs.get(k)
            if df is None or len(df) == 0:
                continue
            valid_keys.append(k)
            n_called = int((df['fdr'] < fdr).sum())
            sig = df[df['fdr'] < fdr]
            n_meta = (int((sig['n_contributing_sig'] == 0).sum())
                      if 'n_contributing_sig' in sig.columns else 0)
            pct_meta = 100.0 * n_meta / n_called if n_called > 0 else 0.0
            mean_abs_lfc = (float(np.abs(sig['logfoldchanges']).mean())
                            if 'logfoldchanges' in sig.columns
                            and len(sig) > 0 else 0.0)
            if 'frac_agreeing' in sig.columns and n_called > 0:
                frac = sig['frac_agreeing'].to_numpy(dtype=float)
                finite = np.isfinite(frac)
                mean_concord = (100.0 * float(frac[finite].mean())
                                if finite.sum() > 0 else 0.0)
                pct_unanimous = (100.0 * float(np.isclose(frac[finite], 1.0).sum())
                                 / n_called)
            else:
                mean_concord = 0.0
                pct_unanimous = 0.0

            # LOO-derived metrics
            res = loo_per_method.get(k, {})
            partitions = res.get('partition_by_n_contributing', {})
            meta_dir = float(partitions.get(0, {}).get(
                'mean_replication_pct_direction_only', float('nan')))
            overall_dir = float(res.get(
                'avg_replication_pct_direction', float('nan')))
            tn = res.get('top_n_per_method', {}).get(k, {})
            auc_z = float(tn.get('mean_auc_sum_of_z', float('nan')))
            auc_u = float(tn.get('mean_auc_ucell', float('nan')))
            # wrong_dir_pct on the FULL consensus set (FDR<threshold),
            # so it's apples-to-apples with the combination rows
            # (which compute wrong-dir on their full called set).
            # In direction_only mode, every testable gene is either
            # replicated or wrong_direction, so:
            #     wrong_dir_pct = 100 - overall_dir_pct
            # The earlier value 'tn['mean_wrong_dir_pct']' was computed
            # on the top-1000 slice, not the full consensus, which made
            # it incomparable to combinations.
            if overall_dir == overall_dir:   # not NaN
                wd = 100.0 - overall_dir
            else:
                wd = float('nan')

            self._rows.append({
                'name': key_to_label.get(k, k),
                'type': 'Single method',
                'methods': [k],
                'rule': 'single',
                'k': 0,
                'n_called':       float(n_called),
                'n_meta':         float(n_meta),
                'pct_meta':       float(pct_meta),
                'mean_abs_lfc':   float(mean_abs_lfc),
                'mean_concord':   float(mean_concord),
                'pct_unanimous':  float(pct_unanimous),
                'meta_dir_pct':   meta_dir,
                'overall_dir_pct': overall_dir,
                'auc_z':          auc_z,
                'auc_u':          auc_u,
                'wrong_dir_pct':  wd,
            })

        # ---- Pairwise intersections + unions ---------------------------
        def _short(key: str) -> str:
            # Use the unambiguous short code so DL / DL-maj / DL-str are
            # distinguishable (same for gwOP / gwOP-str). Fall back to
            # first word of the registry label if a key isn't mapped.
            if key in _METHOD_SHORT:
                return _METHOD_SHORT[key]
            return key_to_label.get(key, key).split(' ')[0]

        # Only compute per-combination LOO when every method has a
        # fold-stash. (True iff per-method LOO was run with
        # stash_fold_details=True -- i.e. after the v2 worker change.)
        stash_ready = (
            bool(loo_per_method)
            and all(
                any('_stash' in f for f in loo_per_method[k].get('folds', []))
                for k in loo_per_method))

        for i in range(len(valid_keys)):
            for j in range(i + 1, len(valid_keys)):
                a, b = valid_keys[i], valid_keys[j]
                for rule, sym, rtype in (
                        ('intersection', ' ∩ ', 'Intersection'),
                        ('union',        ' ∪ ', 'Union')):
                    name = f"{_short(a)}{sym}{_short(b)}"
                    m = _compute_combo_metrics(
                        method_dfs, [a, b], rule, 0, fdr)
                    # Direction-replication LOO from cached stashes.
                    # AUC is not yet computed for combinations (follow-up).
                    if stash_ready and a in loo_per_method and b in loo_per_method:
                        dloo = _combination_direction_loo(
                            loo_per_method, [a, b], rule, 0,
                            replication_threshold=fdr)
                    else:
                        dloo = {'overall_dir_pct': float('nan'),
                                'meta_dir_pct': float('nan'),
                                'wrong_dir_pct': float('nan'),
                                'n_folds': 0}
                    self._rows.append({
                        'name': name,
                        'type': rtype,
                        'methods': [a, b],
                        'rule': rule,
                        'k': 0,
                        'n_called':       float(m['n_called']),
                        'n_meta':         float(m['n_meta']),
                        'pct_meta':       float(m['pct_meta']),
                        'mean_abs_lfc':   float(m['mean_abs_lfc']),
                        'mean_concord':   float(m['mean_concord']),
                        'pct_unanimous':  float(m['pct_unanimous']),
                        'meta_dir_pct':    dloo['meta_dir_pct'],
                        'overall_dir_pct': dloo['overall_dir_pct'],
                        'auc_z':          float('nan'),
                        'auc_u':          float('nan'),
                        'wrong_dir_pct':  dloo['wrong_dir_pct'],
                    })

        # ---- Voting k-of-N over the full set ---------------------------
        n_m = len(valid_keys)
        for k_vote in (3, 5, 7):
            if k_vote > n_m:
                continue
            m = _compute_combo_metrics(
                method_dfs, valid_keys, 'voting', k_vote, fdr)
            # Voting uses all valid_keys; LOO direction replication
            # requires every key to have stashed fold data.
            if (stash_ready
                    and all(k in loo_per_method for k in valid_keys)):
                dloo = _combination_direction_loo(
                    loo_per_method, list(valid_keys), 'voting', k_vote,
                    replication_threshold=fdr)
            else:
                dloo = {'overall_dir_pct': float('nan'),
                        'meta_dir_pct': float('nan'),
                        'wrong_dir_pct': float('nan'),
                        'n_folds': 0}
            self._rows.append({
                'name': f"\u2265{k_vote} of {n_m} methods",
                'type': f"Voting (k={k_vote})",
                'methods': list(valid_keys),
                'rule': 'voting',
                'k': k_vote,
                'n_called':       float(m['n_called']),
                'n_meta':         float(m['n_meta']),
                'pct_meta':       float(m['pct_meta']),
                'mean_abs_lfc':   float(m['mean_abs_lfc']),
                'mean_concord':   float(m['mean_concord']),
                'pct_unanimous':  float(m['pct_unanimous']),
                'meta_dir_pct':    dloo['meta_dir_pct'],
                'overall_dir_pct': dloo['overall_dir_pct'],
                'auc_z':          float('nan'),
                'auc_u':          float('nan'),
                'wrong_dir_pct':  dloo['wrong_dir_pct'],
            })

        # Robust-frontier score: how many of the C(n_metrics, 2) = 55
        # axis pairs is each row non-dominated in? Summary over all
        # possible trade-offs, so user doesn't have to pick axes.
        scores = _robust_frontier_scores(self._rows)
        for r, s in zip(self._rows, scores):
            r['robust_score'] = int(s)

        n_singles = len(valid_keys)
        n_pairs = n_singles * (n_singles - 1) // 2
        n_vote = sum(1 for k_v in (3, 5, 7) if k_v <= n_singles)
        n_metrics = len(_METRICS)
        n_axis_pairs = n_metrics * (n_metrics - 1) // 2
        if not loo_per_method:
            loo_note = (" (no LOO data; run Cross-Study Validation on "
                        "Individual Methods for LOO metrics)")
        elif stash_ready:
            loo_note = (f". LOO data loaded for {len(loo_per_method)} "
                        "methods, including per-fold stashes -- "
                        "combinations show direction-replication LOO "
                        "(AUC for combinations is a follow-up)")
        else:
            loo_note = (f". LOO data loaded for {len(loo_per_method)} "
                        "methods, but from an older run without "
                        "fold stashes -- rerun Cross-Study Validation "
                        "to populate combination LOO metrics")
        self._source_label.setText(
            f"Enumerated {len(self._rows)} rows: {n_singles} singles + "
            f"{n_pairs} pairwise intersections + {n_pairs} pairwise unions + "
            f"{n_vote} voting rules. Robust score = # of {n_axis_pairs} "
            f"axis pairs where row is non-dominated" + loo_note)

        self._refresh_table()
        self._refresh_pareto()

    def _refresh_table(self):
        # Upstream stores percent values as 0-100; widget '.1%' expects
        # 0-1, so divide pct cols once before handing off.
        pct_keys = {k for k, _, _ in _METRICS
                    if k.endswith('_pct') or k in ('mean_concord',
                                                   'pct_unanimous')}
        df_rows = []
        for r in self._rows:
            row = {'name': r['name'], 'type': r['type'],
                   'robust_score': int(r.get('robust_score', 0))}
            for k, _, _ in _METRICS:
                v = r.get(k, float('nan'))
                if v == v and k in pct_keys:
                    v = v / 100.0
                row[k] = v
            df_rows.append(row)

        secondary = QColor(get_color('fg_secondary'))

        def _decorate(item, row, col_i):
            if row['type'] == 'Single method':
                item.setForeground(secondary)

        self._table.set_data(pd.DataFrame(df_rows), decorate=_decorate)
        # Default sort: robust score descending.
        self._table.sortItems(2, Qt.SortOrder.DescendingOrder)

    # ------------------------------------------------------------------
    # Pareto scatter
    # ------------------------------------------------------------------

    def _refresh_pareto(self):
        """Redraw all four small-multiple panels."""
        for plot, cap in self._panel_plots:
            plot.clear_plot_items()
            cap.setText("")
        if not self._rows:
            for plot, _cap in self._panel_plots:
                plot.show_unavailable_message()
            return
        for plot, _cap in self._panel_plots:
            plot.hide_unavailable_message()

        labels = [r['name'] for r in self._rows]

        for (xk, yk, title), (plot, cap) in zip(_PANELS, self._panel_plots):
            x_lbl, x_max = _METRIC_BY_KEY.get(xk, (xk, True))
            y_lbl, y_max = _METRIC_BY_KEY.get(yk, (yk, True))
            pts = []
            for r in self._rows:
                try:
                    x = float(r.get(xk, float('nan')))
                    y = float(r.get(yk, float('nan')))
                except (TypeError, ValueError):
                    x = y = float('nan')
                pts.append((x, y))

            frontier_idx = set(_pareto_frontier(pts, x_max, y_max))

            # Non-frontier points first (grey)
            nf_x = [p[0] for i, p in enumerate(pts)
                    if i not in frontier_idx
                    and np.isfinite(p[0]) and np.isfinite(p[1])]
            nf_y = [p[1] for i, p in enumerate(pts)
                    if i not in frontier_idx
                    and np.isfinite(p[0]) and np.isfinite(p[1])]
            if nf_x:
                plot.plot(nf_x, nf_y, pen=None, symbol='o', symbolSize=7,
                          symbolBrush=(150, 150, 150, 160), symbolPen=None)

            # Frontier points (blue line + markers)
            f_sel = [i for i in frontier_idx
                     if np.isfinite(pts[i][0]) and np.isfinite(pts[i][1])]
            if f_sel:
                f_x = [pts[i][0] for i in f_sel]
                f_y = [pts[i][1] for i in f_sel]
                order = np.argsort(f_x)
                plot.plot([f_x[i] for i in order], [f_y[i] for i in order],
                          pen=pg.mkPen('#1f77b4', width=2),
                          symbol='o', symbolSize=10,
                          symbolBrush=(31, 119, 180), symbolPen=None)
                # Label only frontier points (keeps scatter readable)
                for i in f_sel:
                    txt = pg.TextItem(labels[i], color='#1f77b4',
                                      anchor=(0, 1))
                    txt.setPos(pts[i][0], pts[i][1])
                    plot.addItem(txt)

            plot.setLabel('bottom', x_lbl
                          + (' (lower better)' if not x_max
                             else ' (higher better)'))
            plot.setLabel('left', y_lbl
                          + (' (lower better)' if not y_max
                             else ' (higher better)'))

            n_total = sum(1 for p in pts
                          if np.isfinite(p[0]) and np.isfinite(p[1]))
            cap.setText(
                f"Non-dominated: {len(f_sel)} of {n_total}. "
                f"Frontier: {', '.join(labels[i] for i in f_sel)}")
