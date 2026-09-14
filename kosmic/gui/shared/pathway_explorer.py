# Gene-set selector shared by the DE and Meta workspaces. Takes a
# context flag ('de' or 'meta') for copy + follow-up button variations.

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTreeWidget, QTreeWidgetItem, QHeaderView, QFileDialog,
    QFrame, QDialog, QLineEdit,
    QListWidget, QListWidgetItem, QStackedWidget,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor, QBrush, QFont, QCursor

from kosmic.gui.shared.theme import NoScrollComboBox
from kosmic.gui.shared.widgets import BaseWorker, PrimaryButton, SecondaryButton, SimplePage, TabButton, SecondaryLabel, HintLabel, SectionHeader, CaptionLabel

from kosmic.reference.pathways import builtin_registry
from kosmic.gui.shared import dialogs, run_worker

# Built-in gene set catalogue: (display name, registry key, description).
# Pathways are looked up lazily so missing files surface as a clear error
# instead of an import-time crash.
_BUILTIN_SET_SPECS = [
    ("Energy metabolism -- Core", "metabolic_comprehensive",
     "12 energy pathways: the central flux (glycolysis, PPP, lactate, PDH, FAO, "
     "TCA, OXPHOS with full ETC) plus mitochondrial transport/biogenesis, NAD "
     "metabolism, central amino acid metabolism, and fatty acid synthesis."),
    ("Energy metabolism -- Core (no MT-)", "metabolic_comprehensive_no_mt",
     "As Core, but the mitochondrially-encoded MT- genes are dropped from OXPHOS "
     "so the signal reflects nuclear-encoded genes, not ambient/QC mito contamination."),
    ("Metabolic -- Kirk", "metabolic_kirk",
     "10 pathways including glutamine, ketone body, and nucleotide metabolism."),
    ("Metabolic -- KEGG full (scMetabolism)", "metabolic_kegg_scmetabolism",
     "85 KEGG metabolic pathways (Xiao et al. 2019 / scMetabolism). Select a "
     "focused subset in the next step, or score them all for an unbiased scan."),
    ("GPCRs (HGNC)", "gpcr",
     "95 GPCR subfamilies from HGNC group 139."),
    ("Ion Channels (HGNC)", "ion_channels",
     "41 ion channel families from HGNC group 177."),
    ("Nicotinic Receptors (HGNC)", "nicotinic",
     "3 nicotinic acetylcholine receptor subunit groups."),
    ("Mechanosensing (GO)", "mechanosensing",
     "12 mechanosensing gene groups from Gene Ontology."),
    ("Calcium Signaling (KEGG 2026)", "calcium_signaling",
     "KEGG hsa04020, 250 genes. Flagged as near-significant on preranked "
     "GSEA; committing it here restricts hypothesis-mode FDR to this "
     "pathway alone instead of the whole transcriptome."),
]


def _builtin_sets():
    reg = builtin_registry()
    sets = [(name, reg.get_collection(key), desc)
            for name, key, desc in _BUILTIN_SET_SPECS]

    # Contamination marker panels (from kosmic/reference/contamination) as
    # scoreable gene sets -- run alongside the metabolic pathways as a
    # positive control: if the dominant cell type's markers track disease the
    # same way the metabolic genes do, that points to ambient contamination
    # rather than cell-type-intrinsic biology. Single source of truth (the
    # panel JSON); no gene-list duplication.
    from kosmic.reference.contamination import builtin_panels
    for key, panel in builtin_panels().items():
        genes = list(panel.get('genes', []))
        if not genes:
            continue
        cell = panel.get('dominant_cell_type', key)
        sets.append((
            f"{cell} markers (contamination control)",
            {f"{cell} markers": genes},
            f"{len(genes)}-gene {cell}-specific panel. Score it alongside the "
            "metabolic pathways -- if it moves disease-vs-control the same way, "
            "the metabolic signal may be contamination, not biology.",
        ))
    return sets


BUILTIN_SETS = _builtin_sets()

ENRICHR_CATEGORIES = {
    "KEGG": [
        ("KEGG 2026", "KEGG_2026"),
        ("KEGG 2019 Human", "KEGG_2019_Human"),
    ],
    "Reactome": [
        ("Reactome 2024", "Reactome_Pathways_2024"),
        ("Reactome 2022", "Reactome_2022"),
    ],
    "Gene Ontology": [
        ("GO Biological Process 2025", "GO_Biological_Process_2025"),
        ("GO Cellular Component 2025", "GO_Cellular_Component_2025"),
        ("GO Molecular Function 2025", "GO_Molecular_Function_2025"),
    ],
    "WikiPathways": [
        ("WikiPathways 2024 Human", "WikiPathways_2024_Human"),
    ],
    "MSigDB": [
        ("MSigDB Hallmark 2020", "MSigDB_Hallmark_2020"),
    ],
    "Other": [
        ("BioPlanet 2019", "BioPlanet_2019"),
        ("Panther 2016", "Panther_2016"),
        ("ChEA 2022", "ChEA_2022"),
    ],
}


class _EnrichrFetchWorker(BaseWorker):
    """Background worker to fetch a gene set library from Enrichr.

    Emits 'finished_ok' with a tuple '(display_name, pathways_dict)'.
    """

    def __init__(self, library_id, display_name):
        super().__init__()
        self.library_id = library_id
        self.display_name = display_name

    def _run(self):
        from kosmic.reference.pathways.enrichr import fetch_library
        pathways = fetch_library(self.library_id)
        return (self.display_name, pathways)


class PathwayExplorer(SimplePage):
    """Gene set selector with always-visible coverage table.

    Signals
    -------
    pathways_changed
        Emitted when the active gene set changes. Payload: {pathway: [genes]}.
    """

    help_id = "de/pathway_explorer"

    CONTENT_MARGINS = (20, 20, 20, 20)
    CONTENT_SPACING = 12

    pathways_changed = pyqtSignal(dict)
    confirmed = pyqtSignal()  # emitted when user clicks Confirm (meta context)

    def __init__(self, parent=None, workspace=None, context='de'):
        super().__init__(parent)
        self._workspace = workspace
        self._context = context   # 'de' or 'meta'
        self._current_name = ""
        self._current_pathways = {}
        self._var_names = None
        self._species = 'human'   # detected when set_var_names is called
        self._fetch_worker = None

        self._setup_ui()

        # Load default
        self._select_gene_set(
            "Metabolic -- Comprehensive",
            builtin_registry().get_collection("metabolic_comprehensive"),
        )

    def _setup_ui(self):
        layout = self.body_layout

        # --- Header: current set name + change button ---
        header = QHBoxLayout()
        header.setSpacing(12)

        name_col = QVBoxLayout()
        name_col.setSpacing(2)

        self._set_label = SectionHeader("Active Gene Set")
        name_col.addWidget(self._set_label)

        self._name_label = QLabel("None")
        font = QFont("Segoe UI", 14)
        font.setBold(True)
        self._name_label.setFont(font)
        # Color inherited from global stylesheet
        name_col.addWidget(self._name_label)

        self._summary_label = CaptionLabel("")
        name_col.addWidget(self._summary_label)

        header.addLayout(name_col, 1)

        self._change_btn = PrimaryButton("  Change Gene Set...")
        self._change_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._change_btn.setMinimumHeight(32)
        self._change_btn.setFixedWidth(180)
        self._change_btn.clicked.connect(self._show_browser)
        header.addWidget(self._change_btn, 0, Qt.AlignmentFlag.AlignTop)

        layout.addLayout(header)

        # --- Separator ---
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        layout.addWidget(sep)

        # --- Coverage table (always visible) ---
        # Tree decoration is toggled per-set: flat for collections without
        # an inferred hierarchy, expandable for GPCR / ion-channel sets
        # where parent groups (Adhesion_GPCRs, Calcium_channels, ...) are
        # auto-derived in 'kosmic.reference.pathways.hierarchy'.
        self._coverage_tree = QTreeWidget()
        self._coverage_tree.setHeaderLabels(["Pathway", "Available", "Total", "Coverage", "Missing"])
        self._coverage_tree.setRootIsDecorated(False)
        self._coverage_tree.setAlternatingRowColors(True)
        cov_header = self._coverage_tree.header()
        cov_header.setStretchLastSection(False)
        cov_header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for i in range(1, 5):
            cov_header.setSectionResizeMode(i, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self._coverage_tree, 1)

        # No-data hint (shown when no dataset loaded)
        if self._context == 'meta':
            no_data_text = (
                "Import studies to see gene coverage against "
                "the selected pathways.")
            info_text = (
                "Meta-analysis will run on genes within the "
                "selected pathways only.")
        else:
            no_data_text = (
                "Load a dataset to see gene coverage against "
                "the selected pathways.")
            info_text = (
                "DE runs on all genes in the selected set. "
                "Pathway DE groups results by individual pathway.")

        self._no_data_label = QLabel(no_data_text)
        self._no_data_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._no_data_label.setWordWrap(True)
        self._no_data_label.setProperty("role", "empty_placeholder")
        layout.addWidget(self._no_data_label)

        # Confirm button (meta context: marks workflow step complete)
        if self._context == 'meta':
            self._confirm_btn = PrimaryButton("Confirm Pathways")
            self._confirm_btn.setMinimumHeight(32)
            self._confirm_btn.clicked.connect(self._on_confirm)
            confirm_row = QHBoxLayout()
            confirm_row.addStretch()
            confirm_row.addWidget(self._confirm_btn)
            confirm_row.addStretch()
            layout.addLayout(confirm_row)

        info = HintLabel(info_text)
        info.setWordWrap(True)
        layout.addWidget(info)

    # Core: select a gene set (replaces current)

    def _select_gene_set(self, name, pathways):
        """Set the active gene set. Harmonise gene names, then replace."""
        # Harmonise gene names to current HGNC symbols
        pathways, n_renamed = self._harmonise_gene_set(pathways)

        self._current_name = name
        self._current_pathways = dict(pathways)

        # Update header
        self._name_label.setText(name)
        n_pw = len(pathways)
        n_genes = len({g for genes in pathways.values() for g in genes})
        summary = f"{n_pw} pathways, {n_genes} unique genes"
        if n_renamed > 0:
            summary += f" ({n_renamed} gene names updated to HGNC)"
        self._summary_label.setText(summary)

        # Update coverage
        self._refresh_coverage()

        # Emit
        self.pathways_changed.emit(dict(self._current_pathways))

    # Harmonisation

    _hgnc_lookup = None  # class-level cache

    def _harmonise_gene_set(self, pathways):
        """Map gene names in a pathway dict to current HGNC symbols.

        Returns (harmonised_pathways, n_renamed).
        """
        if PathwayExplorer._hgnc_lookup is None:
            try:
                from kosmic.scrna.inspect.gene_names import load_hgnc_lookup
                PathwayExplorer._hgnc_lookup, _ = load_hgnc_lookup()
            except Exception:
                PathwayExplorer._hgnc_lookup = {}

        lookup = PathwayExplorer._hgnc_lookup
        if not lookup:
            return pathways, 0

        n_renamed = 0
        harmonised = {}
        for pw_name, genes in pathways.items():
            new_genes = []
            for g in genes:
                if g in lookup:
                    new_genes.append(lookup[g])
                    n_renamed += 1
                else:
                    new_genes.append(g)
            harmonised[pw_name] = new_genes
        return harmonised, n_renamed

    # Coverage

    def set_var_names(self, var_names):
        """Set the dataset gene names for coverage calculation.

        Detects the dataset's species from the var names so that human
        pathway lists can be matched against mouse/rat datasets via the
        same naming convention the gene-DE engine uses.
        """
        if var_names is not None:
            self._var_names = set(var_names)
            from kosmic.scrna.inspect.detection import detect_species
            self._species = detect_species(list(var_names))
        else:
            self._var_names = None
            self._species = 'human'
        self._refresh_coverage()

    def _refresh_coverage(self):
        """Update the coverage tree.

        If 'kosmic.reference.pathways.hierarchy.derive_hierarchy' recognises
        the current pathway names, the tree is rendered two-tier (parent
        groups expandable to subfamilies). Otherwise a flat list.
        """
        from kosmic.reference.pathways import derive_hierarchy

        self._coverage_tree.clear()

        if not self._current_pathways:
            self._no_data_label.setText("No gene set selected.")
            self._no_data_label.setVisible(True)
            self._coverage_tree.setVisible(False)
            return

        hierarchy = derive_hierarchy(self._current_pathways)
        self._coverage_tree.setRootIsDecorated(hierarchy is not None)

        if self._var_names is None:
            self._no_data_label.setVisible(True)
            self._coverage_tree.setVisible(True)
            self._populate_without_coverage(hierarchy)
            return

        self._no_data_label.setVisible(False)
        self._coverage_tree.setVisible(True)
        self._populate_with_coverage(hierarchy)

    def _populate_without_coverage(self, hierarchy):
        """Pre-dataset rendering: pathway names + total counts only."""
        if hierarchy is None:
            for pw_name, genes in sorted(self._current_pathways.items()):
                self._coverage_tree.addTopLevelItem(
                    QTreeWidgetItem([pw_name, "—", str(len(genes)), "—", ""])
                )
            return

        grouped_keys = {child for children in hierarchy.values() for child in children}
        # Parent rows + nested children, alphabetised at each level.
        for parent in sorted(hierarchy):
            parent_genes = {
                g for child in hierarchy[parent]
                for g in self._current_pathways.get(child, [])
            }
            parent_item = QTreeWidgetItem([parent, "—", str(len(parent_genes)), "—", ""])
            self._make_bold(parent_item)
            for child in sorted(hierarchy[parent]):
                parent_item.addChild(QTreeWidgetItem(
                    [child, "—", str(len(self._current_pathways[child])), "—", ""]
                ))
            self._coverage_tree.addTopLevelItem(parent_item)

        # Anything not picked up by the hierarchy stays as a flat top-level row.
        for pw_name in sorted(set(self._current_pathways) - grouped_keys):
            genes = self._current_pathways[pw_name]
            self._coverage_tree.addTopLevelItem(
                QTreeWidgetItem([pw_name, "—", str(len(genes)), "—", ""])
            )

    def _populate_with_coverage(self, hierarchy):
        """Rendering with available / total / % / missing per row."""
        if hierarchy is None:
            for pw_name, genes in sorted(self._current_pathways.items()):
                self._coverage_tree.addTopLevelItem(
                    self._make_coverage_item(pw_name, genes)
                )
            return

        grouped_keys = {child for children in hierarchy.values() for child in children}
        for parent in sorted(hierarchy):
            children = hierarchy[parent]
            parent_genes = sorted({
                g for child in children for g in self._current_pathways.get(child, [])
            })
            parent_item = self._make_coverage_item(parent, parent_genes)
            self._make_bold(parent_item)
            for child in sorted(children):
                parent_item.addChild(
                    self._make_coverage_item(child, self._current_pathways[child])
                )
            self._coverage_tree.addTopLevelItem(parent_item)

        for pw_name in sorted(set(self._current_pathways) - grouped_keys):
            self._coverage_tree.addTopLevelItem(
                self._make_coverage_item(pw_name, self._current_pathways[pw_name])
            )

    def _make_coverage_item(self, name, genes):
        """Build a QTreeWidgetItem with available/total/% columns + tooltip.

        Red flags rows the scorer will drop (below the hard floor on
        either gene count or coverage). Orange flags rows that pass but
        sit below the borderline -- still scoreable, but fragile.
        """
        from kosmic import (
            PATHWAY_MIN_GENES, PATHWAY_MIN_COVERAGE,
            PATHWAY_BORDERLINE_COVERAGE,
        )

        var_set = self._var_names
        total = len(genes)
        gene_set = set(genes)
        # Match in the dataset's species namespace (mirrors prepare_gene_coverage
        # in kosmic.de.de_analysis), so human pathway symbols can match a
        # mouse/rat dataset's Title-cased gene names.
        if self._species == 'mouse':
            from kosmic.scrna.inspect.detection import format_gene_for_species
            available_genes = {g for g in gene_set
                               if format_gene_for_species(g, 'mouse') in var_set}
        else:
            available_genes = gene_set & var_set
        available = len(available_genes)
        missing = sorted(gene_set - available_genes)
        coverage = (available / total) if total > 0 else 0.0
        pct = coverage * 100

        pct_str = f"{pct:.0f}%"
        missing_str = ", ".join(missing[:5])
        if len(missing) > 5:
            missing_str += f" (+{len(missing) - 5} more)"

        item = QTreeWidgetItem([name, str(available), str(total), pct_str, missing_str])
        if available < PATHWAY_MIN_GENES or coverage < PATHWAY_MIN_COVERAGE:
            colour = QBrush(QColor("#cc0000"))
        elif coverage < PATHWAY_BORDERLINE_COVERAGE:
            colour = QBrush(QColor("#cc6600"))
        else:
            colour = None
        if colour is not None:
            for col in range(5):
                item.setForeground(col, colour)
        item.setToolTip(4, ", ".join(missing) if missing else "All genes found")
        return item

    @staticmethod
    def _make_bold(item):
        font = item.font(0)
        font.setBold(True)
        for col in range(5):
            item.setFont(col, font)

    # Browser dialog

    def _show_browser(self):
        """Open the gene set browser dialog."""
        dlg = _GeneSetBrowserDialog(self, workspace=self._workspace)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            name, pathways = dlg.result()
            if pathways:
                self._select_gene_set(name, pathways)

    # File import (can also be called from browser dialog)

    def _import_file(self):
        """Import gene sets from a file."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Import Gene Sets", "",
            "Gene set files (*.gmt *.json *.csv);;All files (*)",
        )
        if not path:
            return None, None
        try:
            from pathlib import Path as P
            from kosmic.reference.pathways import formats
            pathways = formats.load_file(path)
            name = P(path).stem
            return name, pathways
        except Exception as e:
            dialogs.warning(self, "Import Failed", f"Could not load file:\n{e}")
            return None, None

    # External API

    def get_flat_pathways(self):
        return dict(self._current_pathways)

    def _on_confirm(self):
        """Confirm the current pathway selection (meta context)."""
        self.confirmed.emit()

# Pathway subset picker (for large fetched / imported libraries)

class _PathwaySubsetDialog(QDialog):
    """Tick a subset of pathways from a large collection (e.g. a fetched KEGG /
    Reactome library) so only the chosen pathways are scored. Returns the
    selected {pathway: genes} subset via 'result_subset'."""

    # Name fragments the "Metabolism preset" button ticks. A heuristic starting
    # point the user reviews -- the authoritative source is the fetched library.
    _METABOLIC_KEYWORDS = (
        'glycolysis', 'gluconeogen', 'citrate cycle', 'tca', 'tricarboxylic',
        'oxidative phosphoryl', 'pentose phosphate', 'fatty acid',
        'pyruvate metabolism', 'ketone', 'valine, leucine', 'branched-chain',
        'glutathione', 'arginine and proline', 'nicotinate', 'nicotinamide',
    )

    def __init__(self, parent, name, pathways):
        super().__init__(parent)
        self._name = name
        self._pathways = dict(pathways)
        self._result = None
        self.setWindowTitle(f"Select pathways -- {name}")
        self.resize(560, 560)
        self._setup_ui()

    def _setup_ui(self):
        lay = QVBoxLayout(self)
        info = SecondaryLabel(
            f"{len(self._pathways):,} pathways in '{self._name}'. Tick the ones "
            "to score. Filter by name and use 'Tick shown', or the Metabolism "
            "preset as a starting point.")
        info.setWordWrap(True)
        lay.addWidget(info)

        self._search = QLineEdit()
        self._search.setPlaceholderText("Filter pathways...")
        self._search.textChanged.connect(self._apply_filter)
        lay.addWidget(self._search)

        self._list = QListWidget()
        for pw in sorted(self._pathways):
            it = QListWidgetItem(f"{pw}  ({len(self._pathways[pw])})")
            it.setData(Qt.ItemDataRole.UserRole, pw)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Unchecked)
            self._list.addItem(it)
        self._list.itemChanged.connect(lambda _i: self._update_count())
        lay.addWidget(self._list, 1)

        row = QHBoxLayout()
        for label, fn in (("Tick shown", self._tick_shown),
                          ("Clear", self._clear_all),
                          ("Metabolism preset", self._metabolic_preset)):
            b = SecondaryButton(label)
            b.clicked.connect(fn)
            row.addWidget(b)
        row.addStretch()
        self._count_label = CaptionLabel("0 selected")
        row.addWidget(self._count_label)
        lay.addLayout(row)

        ok_row = QHBoxLayout()
        ok_row.addStretch()
        cancel = SecondaryButton("Cancel")
        cancel.clicked.connect(self.reject)
        use = PrimaryButton("Use selected")
        use.clicked.connect(self._accept)
        ok_row.addWidget(cancel)
        ok_row.addWidget(use)
        lay.addLayout(ok_row)
        self._update_count()

    def _items(self):
        return (self._list.item(i) for i in range(self._list.count()))

    def _apply_filter(self, text):
        t = text.strip().lower()
        for it in self._items():
            it.setHidden(bool(t) and t not in it.data(Qt.ItemDataRole.UserRole).lower())

    def _tick_shown(self):
        for it in self._items():
            if not it.isHidden():
                it.setCheckState(Qt.CheckState.Checked)

    def _clear_all(self):
        for it in self._items():
            it.setCheckState(Qt.CheckState.Unchecked)

    def _metabolic_preset(self):
        for it in self._items():
            pw = it.data(Qt.ItemDataRole.UserRole).lower()
            if any(k in pw for k in self._METABOLIC_KEYWORDS):
                it.setCheckState(Qt.CheckState.Checked)

    def _update_count(self):
        n = sum(1 for it in self._items()
                if it.checkState() == Qt.CheckState.Checked)
        self._count_label.setText(f"{n} selected")

    def _accept(self):
        chosen = {it.data(Qt.ItemDataRole.UserRole):
                  self._pathways[it.data(Qt.ItemDataRole.UserRole)]
                  for it in self._items()
                  if it.checkState() == Qt.CheckState.Checked}
        if not chosen:
            dialogs.warning(self, "No pathways", "Tick at least one pathway.")
            return
        self._result = chosen
        self.accept()

    def result_subset(self):
        return self._result


# Gene Set Browser Dialog

class _GeneSetBrowserDialog(QDialog):
    """Popup dialog to browse and select a gene set."""

    def __init__(self, explorer: PathwayExplorer, workspace=None):
        super().__init__(explorer)
        self._explorer = explorer
        self._workspace = workspace
        self._selected_name = None
        self._selected_pathways = None
        self._fetch_worker = None

        self.setWindowTitle("Choose Gene Set")
        self.setMinimumSize(550, 450)
        self.resize(600, 500)

        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(0)

        # --- Segmented tab bar ---
        tab_bar = QHBoxLayout()
        tab_bar.setSpacing(0)
        tab_bar.setContentsMargins(0, 0, 0, 0)

        self._tab_buttons = []
        self._stack = QStackedWidget()

        tab_names = ["Built-in", "Import File", "Online (Enrichr)"]
        for i, name in enumerate(tab_names):
            btn = TabButton(name)
            btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
            btn.setCheckable(True)
            btn.clicked.connect(lambda checked, idx=i: self._switch_tab(idx))
            self._tab_buttons.append(btn)
            tab_bar.addWidget(btn)

        self._update_tab_styles()
        tab_bar.addStretch()
        layout.addLayout(tab_bar)
        layout.addSpacing(12)

        # --- Page 0: Built-in ---
        builtin_page = QWidget()
        bl = QVBoxLayout(builtin_page)
        bl.setContentsMargins(0, 0, 0, 0)
        bl.setSpacing(8)

        self._builtin_list = QListWidget()
        self._builtin_list.setSpacing(2)
        self._builtin_list.setAlternatingRowColors(True)
        self._builtin_list.itemDoubleClicked.connect(self._on_builtin_double_click)

        for display_name, gene_dict, description in BUILTIN_SETS:
            n_pw = len(gene_dict)
            n_genes = len({g for genes in gene_dict.values() for g in genes})
            item = QListWidgetItem(
                f"{display_name}\n"
                f"{n_pw} pathways, {n_genes} genes -- {description}"
            )
            item.setData(Qt.ItemDataRole.UserRole, (display_name, gene_dict))
            self._builtin_list.addItem(item)

        self._builtin_list.setCurrentRow(0)
        bl.addWidget(self._builtin_list)

        select_btn = PrimaryButton("Select")
        select_btn.setMinimumHeight(32)
        select_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        select_btn.clicked.connect(self._select_builtin)
        bl.addWidget(select_btn)

        self._stack.addWidget(builtin_page)

        # --- Page 1: Import File ---
        import_page = QWidget()
        il = QVBoxLayout(import_page)
        il.setContentsMargins(0, 0, 0, 0)

        il.addStretch()

        import_info = SecondaryLabel(
            "Import gene sets from a file.\n\n"
            "Supported formats:\n"
            "  GMT — standard pathway format (MSigDB, GSEA, Enrichr exports)\n"
            "  JSON — {pathway: [genes]} dictionary\n"
            "  CSV — two columns: pathway, gene"
        )
        import_info.setWordWrap(True)
        import_info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        il.addWidget(import_info)

        il.addSpacing(16)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        import_btn = PrimaryButton("  Browse for File...")
        import_btn.setMinimumHeight(32)
        import_btn.setFixedWidth(200)
        import_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        import_btn.clicked.connect(self._import_and_accept)
        btn_row.addWidget(import_btn)
        btn_row.addStretch()
        il.addLayout(btn_row)

        il.addStretch()
        self._stack.addWidget(import_page)

        # --- Page 2: Online (Enrichr) ---
        online_page = QWidget()
        ol = QVBoxLayout(online_page)
        ol.setContentsMargins(0, 0, 0, 0)
        ol.setSpacing(8)

        ol.addWidget(QLabel("Download from Enrichr (no login required):"))

        self._enrichr_combo = NoScrollComboBox()
        for category, libs in ENRICHR_CATEGORIES.items():
            for display, lib_id in libs:
                self._enrichr_combo.addItem(f"{category}: {display}", lib_id)
        ol.addWidget(self._enrichr_combo)

        fetch_row = QHBoxLayout()
        self._fetch_btn = PrimaryButton("Download && Use")
        self._fetch_btn.setMinimumHeight(32)
        self._fetch_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self._fetch_btn.clicked.connect(self._fetch_enrichr)
        fetch_row.addWidget(self._fetch_btn)

        self._fetch_status = CaptionLabel("")
        fetch_row.addWidget(self._fetch_status)
        fetch_row.addStretch()
        ol.addLayout(fetch_row)

        enrichr_note = CaptionLabel(
            "Downloads gene set libraries with gene symbols from the Enrichr database.\n"
            "Includes KEGG, Reactome, Gene Ontology, WikiPathways, MSigDB Hallmark, and more.\n\n"
            "No login or API key required."
        )
        enrichr_note.setWordWrap(True)
        ol.addWidget(enrichr_note)

        ol.addStretch()
        self._stack.addWidget(online_page)

        layout.addWidget(self._stack, 1)

        layout.addSpacing(12)

        # Cancel button
        cancel_row = QHBoxLayout()
        cancel_row.addStretch()
        cancel_btn = SecondaryButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        cancel_row.addWidget(cancel_btn)
        layout.addLayout(cancel_row)

        # Select first tab
        self._switch_tab(0)

    def _switch_tab(self, index):
        """Switch to the given tab page."""
        self._stack.setCurrentIndex(index)
        self._update_tab_styles(index)

    def _update_tab_styles(self, active=0):
        """Mark the active tab button.  Styling comes from theme.py."""
        for i, btn in enumerate(self._tab_buttons):
            btn.setChecked(i == active)
            btn.set_active(i == active)

    def result(self):
        """Return (name, pathways) or (None, None) if cancelled."""
        return self._selected_name, self._selected_pathways

    # --- Built-in ---

    def _select_builtin(self):
        item = self._builtin_list.currentItem()
        if item is None:
            return
        name, pathways = item.data(Qt.ItemDataRole.UserRole)
        self._offer_subset(name, pathways)

    def _on_builtin_double_click(self, item):
        name, pathways = item.data(Qt.ItemDataRole.UserRole)
        self._offer_subset(name, pathways)

    # --- Import ---

    def _import_and_accept(self):
        name, pathways = self._explorer._import_file()
        if pathways:
            self._offer_subset(name, pathways)

    def _offer_subset(self, name, pathways):
        """Let the user pick a subset from a large fetched / imported library,
        then accept with just that subset. Small collections pass straight
        through."""
        if len(pathways) <= 12:
            self._selected_name = name
            self._selected_pathways = pathways
            self.accept()
            return
        dlg = _PathwaySubsetDialog(self, name, pathways)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return  # cancelled -> stay in the browser
        subset = dlg.result_subset()
        self._selected_name = (name if len(subset) == len(pathways)
                               else f"{name} ({len(subset)} selected)")
        self._selected_pathways = subset
        self.accept()

    # --- Enrichr ---

    def _fetch_enrichr(self):
        if self._fetch_worker is not None and self._fetch_worker.isRunning():
            return

        lib_id = self._enrichr_combo.currentData()
        display = self._enrichr_combo.currentText()
        if not lib_id:
            return

        self._fetch_btn.setEnabled(False)
        self._fetch_status.setText("Downloading...")

        self._fetch_worker = _EnrichrFetchWorker(lib_id, display)
        run_worker(
            self._fetch_worker,
            on_finished=self._on_enrichr_done,
            on_failed=self._on_enrichr_error,
        )

    def _on_enrichr_done(self, payload):
        display_name, pathways = payload
        self._fetch_btn.setEnabled(True)
        self._fetch_status.setText("")
        self._offer_subset(display_name, pathways)

    def _on_enrichr_error(self, error_msg):
        self._fetch_btn.setEnabled(True)
        self._fetch_status.setText("Download failed")
        dialogs.warning(self, "Download Failed", f"Could not fetch from Enrichr:\n{error_msg}")

