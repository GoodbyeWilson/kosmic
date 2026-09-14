# Home workspace -- KOSMIC's landing page.
#
# JetBrains-style hub:
#   1. Open / New / Recent project actions at the top (always visible).
#   2. Expandable "How does KOSMIC work?" walkthroughs (collapsed by
#      default) for users who want orientation.
#   3. Tutorial / help quick links in the footer.
# Reached via the DNA logo at the top of the nav rail.
from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import QSize, Qt, pyqtSignal
from PyQt6.QtGui import QCursor, QFont
from PyQt6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QPushButton, QScrollArea,
    QToolButton, QVBoxLayout, QWidget,
)

from kosmic.gui.shared.icon_provider import make_icon
from kosmic.gui.shared.theme import get_color
from kosmic.gui.shared.widgets import SecondaryLabel


_WORKFLOWS = (
    {
        'id': 'A',
        'title': 'Run each study end-to-end, then compare',
        'summary': (
            'Process each study separately through scRNA + DE, then pool '
            'the DE results in Meta-Analysis. Best when studies are '
            'small / different tissues, or you trust each study\'s '
            'annotations.'
        ),
        'flow': (
            ('folder-open', 'Project'),
            ('umap',        'scRNA\n(each study)'),
            ('volcano',     'DE\n(each study)'),
            ('forest-plot', 'Meta-Analysis'),
            ('image',       'Figures'),
        ),
        'steps': (
            'Open or create a project (one folder per study inside).',
            'For each study: run the scRNA pipeline (QC, clustering, '
            'annotation) end-to-end.',
            'For each study: run pseudobulk DE in the DE workspace.',
            'Run Meta-Analysis across the per-study DE result CSVs to '
            'get consensus DEGs.',
            'Render publication-quality figures from Figure Export.',
        ),
    },
    {
        'id': 'B',
        'title': 'Unified cell-type labels across studies  (recommended)',
        'summary': (
            'Concatenate all studies into a master dataset, cluster + '
            'annotate once so every cell gets a label from the same model, '
            'split the labels back to each study, then per-study DE with '
            'harmonised cell types. Best for cross-study DE where you '
            'need like-vs-like comparisons.'
        ),
        'flow': (
            ('folder-open', 'Project'),
            ('merge',       'Combine'),
            ('cells',       'Cluster\nmaster'),
            ('split',       'Propagate\nlabels'),
            ('volcano',     'DE\n(each study)'),
            ('forest-plot', 'Meta-Analysis'),
            ('image',       'Figures'),
        ),
        'steps': (
            'Open or create a project with one folder per study.',
            'Click "Combine studies..." on the project page. Select '
            'studies, map condition values to disease / control roles, '
            'build the master h5ad.',
            'Activate the _master row and run the scRNA pipeline on it '
            '-- one clustering + annotation pass covers all studies.',
            'Back on the project page, propagate the master\'s labels '
            'to every per-study h5ad ("Propagate labels..." in the '
            '_master row\'s menu).',
            'For each study: run pseudobulk DE -- now with harmonised '
            'cell-type labels across studies.',
            'Meta-Analysis pools the DE results; Figure Export renders.',
        ),
    },
    {
        'id': 'C',
        'title': 'Meta-analysis from existing DE results',
        'summary': (
            'Already have DE result CSVs from R / Seurat / DESeq2? Drop '
            'them into per-study folders and use only KOSMIC\'s '
            'Meta-Analysis workspace.'
        ),
        'flow': (
            ('folder-open', 'Project'),
            ('dataset',     'Import DE\nCSVs'),
            ('forest-plot', 'Meta-Analysis'),
        ),
        'steps': (
            'Open or create a project; one folder per study.',
            'Drop each study\'s DE result CSV into '
            '{study}/results/de_analysis/statistics/ with the naming '
            'convention {accession}_DE_{method}.csv.',
            'Open Meta-Analysis. Select studies, choose pooling methods '
            '(DL / REML / SumRank / gwOP / Fisher / Stouffer + HKSJ), '
            'run consensus + LOO validation.',
        ),
    },
    {
        'id': 'D',
        'title': 'Make publication figures from existing data',
        'summary': (
            'Already finished analysis elsewhere? Use Figure Export to '
            'render publication-quality volcanos, UMAPs, heatmaps and '
            'forest plots from existing processed h5ads + DE results.'
        ),
        'flow': (
            ('folder-open', 'Project'),
            ('image',       'Figure Export'),
        ),
        'steps': (
            'Open or create a project with your processed h5ads in '
            '{study}/processed_data/ and DE results in '
            '{study}/results/de_analysis/statistics/.',
            'Activate a study and open Figure Export. Pick a figure '
            'type from the sidebar; controls vary per figure.',
        ),
    },
)


class HomeWorkspace(QWidget):
    """KOSMIC landing page: workflow overview + tutorial / help links."""

    tutorial_clicked = pyqtSignal()
    help_clicked = pyqtSignal()

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self) -> None:
        # Scrollable content so smaller windows still reach the footer.
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)

        body = QWidget()
        scroll.setWidget(body)

        outer = QVBoxLayout(body)
        outer.setContentsMargins(48, 36, 48, 36)
        outer.setSpacing(16)

        # ---- Header --------------------------------------------------------
        title = QLabel("KOSMIC")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title_font = QFont("Segoe UI", 22)
        title_font.setBold(True)
        title.setFont(title_font)
        outer.addWidget(title)

        subtitle = SecondaryLabel(
            "Single-cell RNA-seq pipeline + cross-study meta-analysis."
        )
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        outer.addWidget(subtitle)
        outer.addSpacing(20)

        # ---- "How does KOSMIC work?" expandable walkthroughs --------------
        wf_header = QLabel("How does KOSMIC work?")
        wf_header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        wh_font = wf_header.font()
        wh_font.setBold(True)
        wh_font.setPointSize(12)
        wf_header.setFont(wh_font)
        outer.addWidget(wf_header)

        wf_hint = SecondaryLabel(
            "Pick the workflow that fits your data, click to expand "
            "for the step-by-step."
        )
        wf_hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        outer.addWidget(wf_hint)
        outer.addSpacing(4)

        for entry in _WORKFLOWS:
            outer.addWidget(self._build_walkthrough(entry))

        outer.addStretch(1)

        # ---- Footer quick links --------------------------------------------
        footer = QHBoxLayout()
        footer.addStretch()
        self._tutorial_link = QPushButton("Take the tutorial")
        self._tutorial_link.setFlat(True)
        self._tutorial_link.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._tutorial_link.clicked.connect(self.tutorial_clicked.emit)
        footer.addWidget(self._tutorial_link)
        self._help_link = QPushButton("Open help")
        self._help_link.setFlat(True)
        self._help_link.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._help_link.clicked.connect(self.help_clicked.emit)
        footer.addWidget(self._help_link)
        footer.addStretch()
        outer.addLayout(footer)

        # Mount the scroller as this widget's only child.
        wrap = QVBoxLayout(self)
        wrap.setContentsMargins(0, 0, 0, 0)
        wrap.addWidget(scroll)

    # ---- Walkthrough row ---------------------------------------------------

    def _build_walkthrough(self, entry: dict) -> QWidget:
        """Header row (clickable) + collapsible body with diagram + steps."""
        container = QFrame()
        container.setProperty("role", "walkthrough_card")
        col = QVBoxLayout(container)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(0)

        # Header: chevron + title + one-line summary (truncated by space).
        header = QToolButton()
        header.setProperty("role", "walkthrough_header")
        header.setText(f"  {entry['title']}")
        header.setCheckable(True)
        header.setChecked(False)
        header.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        header.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        header.setIcon(make_icon('chevron-right', get_color('fg_secondary'), 14))
        header.setIconSize(QSize(14, 14))
        header.setFixedHeight(34)
        header.setMinimumWidth(0)
        header.setSizePolicy(header.sizePolicy().horizontalPolicy(),
                             header.sizePolicy().verticalPolicy())
        col.addWidget(header)

        # Body (hidden by default).
        body = QWidget()
        body.setVisible(False)
        b_lay = QVBoxLayout(body)
        b_lay.setContentsMargins(28, 8, 16, 14)
        b_lay.setSpacing(10)

        # Summary
        summary = SecondaryLabel(entry['summary'])
        summary.setWordWrap(True)
        b_lay.addWidget(summary)

        # Flow diagram
        flow = self._build_flow_diagram(entry['flow'])
        b_lay.addWidget(flow)

        # Numbered steps
        steps = QVBoxLayout()
        steps.setSpacing(4)
        for i, step_text in enumerate(entry['steps'], 1):
            row = QHBoxLayout()
            row.setSpacing(8)
            n_lbl = QLabel(f"{i}.")
            n_lbl.setFixedWidth(20)
            n_lbl.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)
            n_font = n_lbl.font()
            n_font.setBold(True)
            n_lbl.setFont(n_font)
            row.addWidget(n_lbl, 0, Qt.AlignmentFlag.AlignTop)
            txt = SecondaryLabel(step_text)
            txt.setWordWrap(True)
            row.addWidget(txt, 1)
            steps.addLayout(row)
        b_lay.addLayout(steps)

        col.addWidget(body)

        # Toggle: rotate chevron + show/hide body.
        def _toggled(checked):
            body.setVisible(checked)
            header.setIcon(make_icon(
                'chevron-down' if checked else 'chevron-right',
                get_color('fg_secondary'), 14))
        header.toggled.connect(_toggled)

        return container

    def _build_flow_diagram(self, steps) -> QWidget:
        wrap = QWidget()
        flow = QHBoxLayout(wrap)
        flow.setSpacing(6)
        flow.setContentsMargins(0, 4, 0, 4)
        flow.setAlignment(Qt.AlignmentFlag.AlignLeft)
        accent = get_color('accent_primary')
        fg2 = get_color('fg_secondary')
        for i, (icon_name, label) in enumerate(steps):
            if i > 0:
                arrow = QLabel("→")
                arrow.setStyleSheet(f"color: {fg2}; font-size: 18px;")
                arrow.setAlignment(Qt.AlignmentFlag.AlignCenter)
                flow.addWidget(arrow)
            node = QWidget()
            n = QVBoxLayout(node)
            n.setContentsMargins(0, 0, 0, 0)
            n.setSpacing(2)
            ic = QLabel()
            ic.setFixedSize(32, 32)
            ic.setAlignment(Qt.AlignmentFlag.AlignCenter)
            ic.setPixmap(make_icon(icon_name, accent, 28).pixmap(QSize(28, 28)))
            n.addWidget(ic, 0, Qt.AlignmentFlag.AlignHCenter)
            name = QLabel(label)
            name.setAlignment(Qt.AlignmentFlag.AlignCenter)
            nf = name.font()
            nf.setPointSize(8)
            name.setFont(nf)
            n.addWidget(name, 0, Qt.AlignmentFlag.AlignHCenter)
            flow.addWidget(node)
        flow.addStretch()
        return wrap

    def refresh_theme(self) -> None:
        pass

    def iter_completed_steps(self):
        return iter(())
