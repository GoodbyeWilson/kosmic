# LOO Validation tab: held-out replication + method benchmark.
#
# Parent <-> tab contract:
# * 'characterise_clicked' / 'benchmark_clicked(level)' signals out
# * 'populate' / 'populate_benchmark' in
# * 'n_methods_getter' / 'mode_getter' closures injected at construction
from __future__ import annotations

from typing import Callable, Optional

import pandas as pd
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import (
    QAbstractItemView, QHBoxLayout, QLabel, QMenu,
    QPushButton, QSplitter, QToolButton, QVBoxLayout, QWidget,
)

from kosmic.gui.shared import borderless
from kosmic.gui.shared.widgets import (
    Column, ResultsTable, SecondaryLabel, SettingsGroup,
)


_DEFAULT_STATUS = (
    "Run a normal consensus first; then click Run LOO. Each fold "
    "rebuilds the consensus on K-1 studies and tests whether the "
    "held-out study independently calls those genes. With CC perm "
    "calibration this can take K x your normal run time.")

_DEFAULT_GENE_PROMPT = (
    "Click a fold above to see its per-gene replication.")

_DEFAULT_BENCH_STATUS = (
    "Click 'Benchmark Methods' above to compare the consensus "
    "against each single constituent method on the same K folds.")

# Order of base methods used in the States column. Must match
# _BASE_METHOD_KEYS in kosmic.meta_analysis.consensus_loo.
_BENCHMARK_BASE_METHODS = ('dl', 'sumrank', 'gwop')

_BENCHMARK_LEVELS = [
    ('current',
     "Chosen methods only (each alone + their consensus, ~5 min)",
     "Tests just the methods you've ticked above, plus the "
     "consensus of all of them. Fastest option."),
    ('quick',
     "Quick combinatorial (~10 configs, ~10 min)",
     "All method singles (both Top50 states each) + all "
     "multi-method subsets at your current Top50 settings. "
     "Good first pass for 3 methods."),
    ('full',
     "Full combinatorial (26 configs, ~26 min)",
     "Each method has 3 states: Off / On (all studies) / "
     "On (Top 50%). 3^3 - 1 = 26 non-empty configs for the "
     "3 methods. The principled answer for 'which "
     "configuration is best'."),
    ('exhaustive',
     "Exhaustive (63 configs, ~1 hour)",
     "Treats each method+Top50 variant as independent; "
     "allows combinations like DL_AllStudies + DL_Top50 in "
     "the same consensus. 2^6 - 1 = 63. Most thorough."),
]


class LOOTab(QWidget):
    """LOO Validation tab UI + per-fold drill-down."""

    benchmark_clicked = pyqtSignal(str)   # level
    characterise_clicked = pyqtSignal()

    def __init__(
        self,
        n_methods_getter: Callable[[], int],
        mode_getter: Callable[[], str],
        parent: Optional[QWidget] = None,
    ):
        super().__init__(parent)
        self._n_methods_getter = n_methods_getter
        self._mode_getter = mode_getter
        self._result: Optional[dict] = None
        self._build_ui()

    # -- UI construction -----------------------------------------------

    def _build_ui(self) -> None:
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)

        _selectable = (Qt.TextInteractionFlag.TextSelectableByMouse
                       | Qt.TextInteractionFlag.TextSelectableByKeyboard)

        # Top control row.
        top = QHBoxLayout()

        self.characterise_btn = QPushButton("Characterise Replicators")
        self.characterise_btn.setToolTip(
            "After running an LOO, click to pool per-gene features "
            "across folds and compare |log2FC|, heterogeneity (I2), "
            "n_studies, and SE between replicated genes and wrong-"
            "direction failures. Opens a dialog with the comparison.")
        self.characterise_btn.clicked.connect(self.characterise_clicked.emit)
        top.addWidget(self.characterise_btn)

        self.benchmark_btn = QToolButton()
        self.benchmark_btn.setText("Benchmark Methods ▼")
        self.benchmark_btn.setPopupMode(
            QToolButton.ToolButtonPopupMode.InstantPopup)
        self.benchmark_btn.setToolTip(
            "Compare LOO replication across method configurations. "
            "Choose a sweep depth from the dropdown. All sweeps use "
            "analytical calibration (CC perm x N configs would take "
            "many hours).")
        bench_menu = QMenu(self.benchmark_btn)
        for level, label, tip in _BENCHMARK_LEVELS:
            act = QAction(label, self.benchmark_btn)
            act.setToolTip(tip)
            act.triggered.connect(
                lambda _checked=False, lv=level:
                    self.benchmark_clicked.emit(lv))
            bench_menu.addAction(act)
        self.benchmark_btn.setMenu(bench_menu)
        top.addWidget(self.benchmark_btn)

        lay.addLayout(top)

        # Status + headline.
        self.status_label = SecondaryLabel(_DEFAULT_STATUS)
        self.status_label.setWordWrap(True)
        self.status_label.setTextInteractionFlags(_selectable)
        lay.addWidget(self.status_label)

        self.headline_label = QLabel("")
        self.headline_label.setWordWrap(True)
        self.headline_label.setProperty("role", "headline")
        self.headline_label.setTextInteractionFlags(_selectable)
        lay.addWidget(self.headline_label)

        # Splitter: per-fold summary on top, drill-down per-gene below.
        splitter = QSplitter(Qt.Orientation.Vertical)

        fold_box = QWidget()
        fold_lay = borderless(QVBoxLayout, fold_box)
        fold_lay.addWidget(QLabel("Per-fold summary (click a row to drill in):"))
        self.fold_table = ResultsTable()
        self.fold_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.fold_table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self.fold_table.setSortingEnabled(True)
        self.fold_table.setAlternatingRowColors(True)
        self.fold_table.itemSelectionChanged.connect(self._on_fold_selected)
        fold_lay.addWidget(self.fold_table)
        splitter.addWidget(fold_box)

        gene_box = QWidget()
        gene_lay = borderless(QVBoxLayout, gene_box)
        self.gene_label = QLabel(_DEFAULT_GENE_PROMPT)
        self.gene_label.setTextInteractionFlags(_selectable)
        gene_lay.addWidget(self.gene_label)
        self.gene_table = ResultsTable()
        self.gene_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.gene_table.setSortingEnabled(True)
        self.gene_table.setAlternatingRowColors(True)
        gene_lay.addWidget(self.gene_table)
        splitter.addWidget(gene_box)

        splitter.setSizes([300, 400])
        lay.addWidget(splitter, stretch=1)

        # Method benchmark section.
        self.bench_box = SettingsGroup(
            "Method Benchmark (held-out replication per method)")
        self.bench_box.setToolTip(
            "Click 'Benchmark Methods' above to run LOO separately "
            "for each single-method configuration and compare against "
            "the full consensus. Tests whether the consensus actually "
            "improves reproducibility over its parts.")
        bench_lay = borderless(QVBoxLayout)
        self.bench_box.add_layout(bench_lay)
        self.benchmark_status = SecondaryLabel(_DEFAULT_BENCH_STATUS)
        self.benchmark_status.setWordWrap(True)
        self.benchmark_status.setTextInteractionFlags(_selectable)
        bench_lay.addWidget(self.benchmark_status)

        self.benchmark_table = ResultsTable()
        self.benchmark_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.benchmark_table.setSortingEnabled(True)
        self.benchmark_table.setAlternatingRowColors(True)
        self.benchmark_table.setMaximumHeight(200)
        bench_lay.addWidget(self.benchmark_table)
        lay.addWidget(self.bench_box)

    # -- public API ----------------------------------------------------

    @property
    def result(self) -> Optional[dict]:
        """
        Most recent LOO result, or 'None' if never run.

        The replicator-characterisation dialog (kept on the parent
        because it needs cross-workspace state) reads this.
        """
        return self._result

    def set_status(self, text: str) -> None:
        self.status_label.setText(text)

    def set_headline(self, text: str) -> None:
        self.headline_label.setText(text)

    def set_benchmark_status(self, text: str) -> None:
        self.benchmark_status.setText(text)

    def clear(self) -> None:
        """Reset everything to the empty state (called on dataset reload)."""
        self._result = None
        self.fold_table.setRowCount(0)
        self.fold_table.setColumnCount(0)
        self.gene_table.setRowCount(0)
        self.gene_table.setColumnCount(0)
        self.benchmark_table.setRowCount(0)
        self.benchmark_table.setColumnCount(0)
        self.set_status(_DEFAULT_STATUS)
        self.set_headline("")
        self.gene_label.setText(_DEFAULT_GENE_PROMPT)
        self.set_benchmark_status(_DEFAULT_BENCH_STATUS)

    def populate(self, result: dict) -> None:
        """Render the LOO result into the headline + per-fold table."""
        from kosmic.meta_analysis.consensus_loo import REPLICATION_LABELS

        self._result = result

        mode_label = REPLICATION_LABELS.get(
            result['replication_mode'], result['replication_mode'])
        is_hypothesis = self._mode_getter() == 'hypothesis'

        # -- Per-method top-N fair-comparison summary (discovery only).
        tn_summary = result.get('top_n_per_method', {})
        tn_size = result.get('top_n', 1000)
        tn_lines = ""
        if tn_summary and not is_hypothesis:
            tn_parts = []
            for mkey, stats in tn_summary.items():
                auc_z = stats.get('mean_auc_sum_of_z', float('nan'))
                auc_u = stats.get('mean_auc_ucell', float('nan'))
                auc_suffix = ""
                # Only show AUC if it was computed (pseudobulk available).
                # 'x == x' filters NaN without importing numpy.
                if auc_z == auc_z or auc_u == auc_u:
                    parts = []
                    if auc_z == auc_z:
                        parts.append(f"AUCz {auc_z:.2f}")
                    if auc_u == auc_u:
                        parts.append(f"AUCu {auc_u:.2f}")
                    auc_suffix = f", {' / '.join(parts)}"
                tn_parts.append(
                    f"{mkey}: {stats['mean_replication_pct']:.1f}% "
                    f"(wd {stats['mean_wrong_dir_pct']:.1f}%{auc_suffix})"
                )
            tn_lines = (
                f"\nPer-method top-{tn_size} (fair comparison, fixed gene-set "
                f"size):  " + "  |  ".join(tn_parts)
            )

        # -- Consensus-level AUC line.
        auc_line = ""
        avg_auc_z = result.get('avg_auc_sum_of_z', float('nan'))
        avg_auc_u = result.get('avg_auc_ucell', float('nan'))
        if avg_auc_z == avg_auc_z or avg_auc_u == avg_auc_u:
            parts = []
            if avg_auc_z == avg_auc_z:
                parts.append(
                    f"AUC (sum-of-signed-z): {avg_auc_z:.3f} "
                    f"(range {result['min_auc_sum_of_z']:.3f}-"
                    f"{result['max_auc_sum_of_z']:.3f})"
                )
            if avg_auc_u == avg_auc_u:
                parts.append(
                    f"AUC (UCell-rank): {avg_auc_u:.3f} "
                    f"(range {result['min_auc_ucell']:.3f}-"
                    f"{result['max_auc_ucell']:.3f})"
                )
            auc_line = (
                "\nHeld-out case/control prediction on the consensus gene "
                "list:  " + "  |  ".join(parts))

        # -- Per-method SOLO meta-discovery counts.
        solo_sum = result.get('solo_per_method', {})
        solo_line = ""
        if solo_sum and not is_hypothesis:
            parts = []
            for mkey, s in solo_sum.items():
                n_meta = s.get('mean_n_meta_discovered', 0.0)
                n_total = s.get('mean_n_total', 0.0)
                meta_dir = s.get('mean_meta_direction_pct', 0.0)
                if n_meta >= 0.5:
                    parts.append(
                        f"{mkey}: {n_meta:.0f}/{n_total:.0f} meta "
                        f"(dir {meta_dir:.1f}%)"
                    )
                else:
                    parts.append(f"{mkey}: 0/{n_total:.0f} meta")
            solo_line = (
                "\nPer-method SOLO meta-discoveries (FDR<0.05 standalone, "
                "not intersection):  " + "  |  ".join(parts))

        # -- Meta-discovered vs obvious partition.
        part_obvious = result.get('partition_obvious', {})
        part_meta = result.get('partition_meta_discovered', {})
        part_line = ""
        if part_obvious or part_meta:
            n_obv = part_obvious.get('mean_n', 0.0)
            n_meta = part_meta.get('mean_n', 0.0)
            dir_obv = part_obvious.get('mean_replication_pct_direction_only', 0.0)
            dir_meta = part_meta.get('mean_replication_pct_direction_only', 0.0)
            strict_obv = part_obvious.get('mean_replication_pct_strict', 0.0)
            strict_meta = part_meta.get('mean_replication_pct_strict', 0.0)
            meta_label = ("Pooling-only" if is_hypothesis
                          else "Meta-discovered")
            part_line = (
                f"\n{meta_label} ({n_meta:.0f} genes, FDR-sig in 0 "
                f"contributing studies): direction {dir_meta:.1f}%, "
                f"strict {strict_meta:.1f}%  |  "
                f"Individually significant ({n_obv:.0f} genes, FDR-sig "
                f"in >=1 study): direction {dir_obv:.1f}%, "
                f"strict {strict_obv:.1f}%"
            )

        # -- Granular per-n-study breakdown.
        part_by_n = result.get('partition_by_n_contributing', {})
        by_n_line = ""
        if part_by_n:
            by_n_parts = []
            for n_sig in sorted(part_by_n.keys()):
                stats = part_by_n[n_sig]
                if stats['mean_n'] < 0.5:
                    continue
                by_n_parts.append(
                    f"[{n_sig} studies sig: {stats['mean_n']:.0f} genes, "
                    f"dir {stats['mean_replication_pct_direction_only']:.1f}%, "
                    f"strict {stats['mean_replication_pct_strict']:.1f}%]"
                )
            if by_n_parts:
                by_n_line = ("\nBy n contributing studies FDR-sig:  "
                             + "  ".join(by_n_parts))

        self.set_headline(
            f"Average replication across {result['n_studies']} folds "
            f"(rates exclude missing genes from the denominator):  "
            f"Strict {result['avg_replication_pct_strict']:.1f}% "
            f"(range {result['min_replication_pct_strict']:.1f}-"
            f"{result['max_replication_pct_strict']:.1f}%)  |  "
            f"Loose {result['avg_replication_pct_loose']:.1f}% "
            f"(range {result['min_replication_pct_loose']:.1f}-"
            f"{result['max_replication_pct_loose']:.1f}%)  |  "
            f"Direction-only {result['avg_replication_pct_direction']:.1f}% "
            f"(range {result['min_replication_pct_direction']:.1f}-"
            f"{result['max_replication_pct_direction']:.1f}%)\n"
            f"Drill-down (per-gene) reasons reflect the selected primary "
            f"mode: {mode_label}"
            + auc_line + solo_line + tn_lines + part_line + by_n_line)

        # -- Per-fold table.
        n_methods = self._n_methods_getter()
        if n_methods <= 1:
            sig_col_name = 'Significant'
            sig_col_tip = ("Genes called significant (FDR < 0.05) in "
                           "the K-1 fold. This is the gene set the "
                           "held-out study is asked to replicate.")
        else:
            sig_col_name = 'Consensus'
            sig_col_tip = ("Genes called significant in the K-1 fold "
                           "consensus (FDR < 0.05 in all methods). "
                           "This is the gene set the held-out study "
                           "is asked to replicate.")
        cols = [
            Column('Held out',     'study_held_out',            's',
                   "Study held out from the K-1 fold. The held-out "
                   "study's standalone DE result is the test set."),
            Column('K-1',          'n_studies_in_fold',         'd',
                   "Number of studies used to build this fold's analysis."),
            Column(sig_col_name,   'consensus_size',            'd', sig_col_tip),
            Column('Missing',      'n_missing',                 'd',
                   "Consensus genes not measured in the held-out study (e.g. "
                   "filtered out by per-study expression cutoffs). EXCLUDED "
                   "from rate denominators -- a gene that wasn't measured "
                   "can't have its direction or p-value tested."),
            Column('Strict %',     'replication_pct_strict',    '.1f',
                   "n_replicated / n_testable, where: replicated = held-out "
                   "FDR < 0.05 AND same direction; n_testable = consensus "
                   "size minus missing."),
            Column('Loose %',      'replication_pct_loose',     '.1f',
                   "n_replicated / n_testable, where: replicated = held-out "
                   "nominal p < 0.05 AND same direction; n_testable = "
                   "consensus size minus missing."),
            Column('Direction %',  'replication_pct_direction', '.1f',
                   "n_replicated / n_testable, where: replicated = same sign "
                   "of log2FC, no p-value test; n_testable = consensus size "
                   "minus missing. This is the upper bound on replication "
                   "when individual studies are underpowered."),
        ]

        # Append AUC columns if any fold computed them.
        any_auc_z = any(
            (f.get('auc_consensus') or {}).get('auc_sum_of_z', float('nan'))
            == (f.get('auc_consensus') or {}).get('auc_sum_of_z', float('nan'))
            for f in result['folds'])
        any_auc_u = any(
            (f.get('auc_consensus') or {}).get('auc_ucell', float('nan'))
            == (f.get('auc_consensus') or {}).get('auc_ucell', float('nan'))
            for f in result['folds'])
        if any_auc_z:
            cols.append(Column('AUC (z)', '_auc_z', '.3f',
                "Held-out case/control AUC using sum-of-signed-z module "
                "score on the consensus gene list. Higher = consensus "
                "genes distinguish disease from control in the held-out "
                "study. AUC < 0.5 = inverse predictor (sign-flip pathology)."))
        if any_auc_u:
            cols.append(Column('AUC (UCell)', '_auc_u', '.3f',
                "Held-out case/control AUC using rank-based UCell-style "
                "module score. Provided for direct parity with SumRank "
                "paper's validation. Non-parametric, ignores effect "
                "magnitude."))

        # Flatten fold dicts to a DataFrame so the widget can iterate.
        for f in result['folds']:
            au = f.get('auc_consensus') or {}
            f['_auc_z'] = au.get('auc_sum_of_z', float('nan'))
            f['_auc_u'] = au.get('auc_ucell', float('nan'))
        folds_df = pd.DataFrame(result['folds'])

        def _decorate(item, row, col_i):
            # Stash fold index on UserRole+1 of the first column so the
            # drill-down handler can look it up after sorting.
            if col_i == 0 and 'fold_index' in row.index:
                item.setData(Qt.ItemDataRole.UserRole + 1,
                             int(row['fold_index']))
            # AUC cols carry '-' for missing instead of the default 'N/A'
            # (folds without pseudobulk can't compute AUC).
            if cols[col_i].key in ('_auc_z', '_auc_u') \
                    and item.text() == 'N/A':
                item.setText('-')

        self.fold_table.set_schema(cols)
        self.fold_table.set_data(folds_df, decorate=_decorate)
        # Auto-select the first row to populate the drill-down.
        if self.fold_table.rowCount() > 0:
            self.fold_table.selectRow(0)

    def populate_benchmark(self, results: list) -> None:
        """
        Render the per-config comparison table: one row per method
        configuration, columns for each replication mode's avg / range
        and the consensus gene count.
        """

        def _format_states(states):
            if not states:
                return ''
            return ' / '.join(
                states.get(b, '-') for b in _BENCHMARK_BASE_METHODS)

        rows_data = []
        for r in results:
            lr = r['loo_result']
            n_folds = max(1, len(lr['folds']))
            avg_n = sum(f['consensus_size'] for f in lr['folds']) / n_folds
            rows_data.append({
                'config_label': r['config_label'],
                'n_methods':    r.get('n_methods', len(r.get('method_keys', []))),
                'states_str':   _format_states(r.get('states')),
                'avg_genes':    avg_n,
                'strict_avg':   lr['avg_replication_pct_strict'],
                'strict_range': (
                    f"{lr['min_replication_pct_strict']:.0f}-"
                    f"{lr['max_replication_pct_strict']:.0f}%"),
                'loose_avg':    lr['avg_replication_pct_loose'],
                'direction_avg': lr['avg_replication_pct_direction'],
            })

        self.benchmark_table.set_schema([
            Column('Config',          'config_label',  's',
                   "Method configuration tested in this row."),
            Column('N',               'n_methods',     'd',
                   "Number of methods in this configuration."),
            Column('States',          'states_str',    's',
                   "Per-method state: '-' (off), 'all' (all-studies), "
                   "'top50'. Order: DL / SumRank / gwOP."),
            Column('Avg #genes',      'avg_genes',     'd',
                   "Mean consensus size (genes passing FDR < 0.05) "
                   "across the K folds."),
            Column('Strict avg %',    'strict_avg',    '.1f',
                   "Mean held-out replication under STRICT mode "
                   "(held-out FDR<0.05 + same direction). Testable "
                   "denominator excludes missing genes."),
            Column('Strict range',    'strict_range',  's',
                   "min-max across the K folds under strict mode."),
            Column('Loose avg %',     'loose_avg',     '.1f',
                   "Mean held-out replication under LOOSE mode "
                   "(held-out p<0.05 + same direction)."),
            Column('Direction avg %', 'direction_avg', '.1f',
                   "Mean held-out replication under DIRECTION-ONLY mode "
                   "(same sign of log2FC, no p-value test)."),
        ])
        self.benchmark_table.set_data(pd.DataFrame(rows_data))

    # -- internal: drill-down ------------------------------------------

    def _on_fold_selected(self) -> None:
        """When a fold row is clicked, show its per-gene replication."""
        if self._result is None:
            return
        rows = self.fold_table.selectionModel().selectedRows()
        if not rows:
            return
        row = rows[0].row()
        # Pull the fold index back from UserRole+1 on column 0.
        item = self.fold_table.item(row, 0)
        if item is None:
            return
        fold_idx = item.data(Qt.ItemDataRole.UserRole + 1)
        if fold_idx is None:
            return
        try:
            fold = self._result['folds'][int(fold_idx)]
        except (IndexError, TypeError, ValueError):
            return
        self._populate_gene_table(fold)

    def _populate_gene_table(self, fold: dict) -> None:
        """Render one fold's per-gene replication table."""
        pg = fold['per_gene_df']
        self.gene_label.setText(
            f"Held out: {fold['study_held_out']}  -- "
            f"{fold['n_replicated']}/{fold['n_testable']} "
            f"replicated ({fold['replication_pct']:.1f}%, "
            f"{fold['n_missing']} missing excluded)")

        self.gene_table.set_schema([
            Column('Gene',         'gene',             's',
                   "Gene symbol."),
            Column('Cons. log2FC', 'consensus_log2FC', '.3f',
                   "Pooled log2 fold-change from this fold's K-1 consensus "
                   "(does NOT include the held-out study)."),
            Column('Held log2FC',  'study_log2FC',     '.3f',
                   "Single-study log2 fold-change from the held-out study's "
                   "standalone DE result."),
            Column('Held p',       'study_pval',       '.2e',
                   "Held-out study's nominal per-gene p-value (used by Loose "
                   "mode)."),
            Column('Held FDR',     'study_padj',       '.2e',
                   "Held-out study's BH-adjusted p-value (used by Strict "
                   "mode)."),
            Column('Replicated',   'replicated',       'bool',
                   "Did the gene pass the selected primary mode's "
                   "replication test?"),
            Column('Reason',       'reason',           's',
                   "replicated / wrong_direction / not_significant / "
                   "missing. 'wrong_direction' suggests a sign-convention or "
                   "normalisation issue in the held-out study's standalone "
                   "DE; 'not_significant' typically reflects per-study "
                   "underpowering (the meta-analysis pools many such "
                   "underpowered signals)."),
        ])
        self.gene_table.set_data(pg)


__all__ = ["LOOTab"]
