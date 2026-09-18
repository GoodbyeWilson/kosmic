##################################################
## KOSMIC -- Kirk's Open-source Single-cell Meta-analysis Integration and Comparison
##
## PyQt6 application; one window, five workspaces (Project, scRNA
## Analysis, Differential Expression, Meta-Analysis, Figures) plus
## Home. See CODEBASE.md for the layout and CONTRIBUTING.md for the
## rules. Finished work is recorded in git and CODEBASE.md, not here.
##################################################

# TODO (open items only; add and remove as things land)
#  1. Release prerequisites: OSI licence (pyproject.toml), CITATION.cff,
#     Zenodo-archived release, the tutorial studies table in README.md.
#  2. Help pages still unwritten: meta/enrichment, meta/methods,
#     meta/validation (tests/test_help_content.py KNOWN_MISSING); the
#     DE, Meta and Figures sections need the same concept-first pass
#     the Project and scRNA sections had.
#  3. Unit tests for the large GUI pages (gene_de_page, cluster_tab,
#     gene_ma_page) beyond construct-and-show.
#  4. GWAS overlap: replace nearest-gene with MAGMA; the Validation
#     step's GWAS tab is a placeholder until then (dev/bench/gwas_overlap.py).
#  5. Meta vs mega comparison view. The arms' gene universes differ
#     (each applies its own detection filter), so intersect before
#     counting. Decides the paper's venue.
#  6. ADR-001 Phase B: per-module study selection replacing the global
#     active study (docs/ADR-001-project-source-of-truth.md).
#  7. Pathway-level CC permutation (kosmic/meta_analysis/cc_permutation.py
#     pathway_cc_permutation) is tested but not wired into the pathway
#     meta-analysis page.
#  8. DE detection default (config de.detection_min_pct = 0.05) is high
#     for snRNA-seq; consider 0.01. The Filters card exposes six
#     controls where only min-cells-per-donor is worth varying.
#  9. Gene Groups tab: import / Enrichr picker and harmonisation of the
#     chosen set; module-score UMAP overlay.
# 10. Nothing surfaces gene-merge asymmetry outside the Combine dialog;
#     a Project-page line would help.
# 11. Run the remaining DCM studies (chaffin, koenig, reichart) through
#     DE to reach k = 5 for the meta-analysis.

# System Imports
import os
# Set BLAS thread caps before anything that might pull numpy.
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import contextlib
import sys
from pathlib import Path
from typing import ClassVar

try:
    import psutil
except ImportError:
    psutil = None

# Force Agg backend for matplotlib for macOS compatibility
import matplotlib
matplotlib.use('Agg')

from kosmic.paths import scaffold_project, scaffold_study

# PyQt Imports
from PyQt6.QtCore import Qt, QTimer, QSettings, pyqtSignal
from PyQt6.QtGui import QGuiApplication, QAction
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel,
    QHBoxLayout, QStackedWidget, QFileDialog, QSplitter,
)

# KOSMIC!!
from kosmic.gui.shared import borderless, dialogs
from kosmic.gui.shared.theme import apply_theme, get_current_mode, ThemeMode
from kosmic.gui.shared.icon_provider import make_icon
from kosmic.gui.shared.splash import SplashScreen
from kosmic.gui.shared.explorer import ExplorerPane
from kosmic.gui.shared.nav_rail import NavigationRail
from kosmic.gui.shared.output_panel import OutputPanel

# Used for accurate scaling across monitors with different DPI.
QGuiApplication.setHighDpiScaleFactorRoundingPolicy(
    Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
)

# Stack page indices. Visual order in the nav rail is independent
# (defined by NavigationRail._MODE_ITEMS); the rail emits a mode key
# and _MODE_DISPATCH maps it to one of these.
PAGE_SCRNA = 0
PAGE_DE = 1
PAGE_FIGURES = 2
PAGE_META = 3
PAGE_PROJECT = 4
PAGE_HOME = 5


# Our main app. All other files moved out of here to clean up
class AppWindow(QMainWindow):
    """Main application window with navigation rail and stacked tabs/workspaces."""

    # Project state signals: workspaces subscribe to react to user picks.
    project_changed = pyqtSignal(object)        # Path or None
    study_changed = pyqtSignal(object)          # Path or None
    study_list_changed = pyqtSignal()

    def __init__(self, splash_status_cb=None):
        super().__init__()
        self._splash_status_cb = splash_status_cb
        self.settings = QSettings("KOSMIC", "KOSMIC")
        self._is_dark = get_current_mode() == ThemeMode.DARK
        self._scrna_loaded = False
        self._de_loaded = False
        self._figures_loaded = False
        self._meta_loaded = False
        self._project_loaded = False
        self._home_loaded = False
        self._scrna_workspace = None
        self._de_workspace = None
        self._figures_workspace = None
        self._meta_workspace = None
        self._project_workspace = None
        self._home_workspace = None

        self.current_project_dir: Path | None = None
        self.current_accession: str | None = None

        self._setup_ui()
        self._build_menu_bar()
        self._load_geometry()

    @property
    def current_study_path(self) -> Path | None:
        """Active study folder, or None if no project / study is set."""
        if self.current_project_dir is None or self.current_accession is None:
            return None
        return self.current_project_dir / self.current_accession

    def set_project_directory(self, path) -> None:
        """Activate ``path`` as the project. Scaffolds meta_analysis/; clears study selection."""
        new_path = Path(path) if path else None
        if new_path == self.current_project_dir:
            return
        self.current_project_dir = new_path
        self.current_accession = None
        if new_path is not None:
            scaffold_project(new_path)
            self.settings.setValue("last_project_dir", str(new_path))
            self._push_recent_project(new_path)
        self.project_changed.emit(new_path)
        self.study_changed.emit(None)

    _RECENT_PROJECTS_LIMIT = 8

    def _push_recent_project(self, path: Path) -> None:
        """Move ``path`` to the front of the recent-projects list."""
        recents = self._load_recent_projects()
        s = str(path)
        recents = [p for p in recents if p != s]
        recents.insert(0, s)
        recents = recents[: self._RECENT_PROJECTS_LIMIT]
        self.settings.setValue("kosmic/recent_projects", recents)
        if self._project_workspace is not None:
            self._project_workspace.set_recent_projects(recents)

    def _load_recent_projects(self) -> list[str]:
        raw = self.settings.value("kosmic/recent_projects", []) or []
        if isinstance(raw, str):
            return [raw] if raw else []
        return [str(p) for p in raw if p]

    def add_study(self, accession: str) -> None:
        """Scaffold ``accession`` as a study folder under the current project."""
        if self.current_project_dir is None:
            raise RuntimeError("No project directory set")
        if not accession:
            raise ValueError("Empty accession")
        scaffold_study(self.current_project_dir / accession)
        self.study_list_changed.emit()

    def set_active_study(self, accession) -> None:
        """Activate ``accession`` as the study read by per-study workspaces."""
        new = accession if accession else None
        if new == self.current_accession:
            return
        self.current_accession = new
        self.study_changed.emit(self.current_study_path)

    def _setup_ui(self):
        self.setWindowTitle("KOSMIC — Kirk's Open-source Single-cell Meta-analysis Integration Consensus")
        self.setMinimumSize(1000, 700)

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        self._nav_rail = NavigationRail()
        self._nav_rail.mode_selected.connect(self._on_mode_selected)
        self._nav_rail.theme_toggled.connect(self._toggle_theme)
        self._nav_rail.tutorial_clicked.connect(self._start_tutorial)
        self._nav_rail.home_clicked.connect(
            lambda: self._on_mode_selected("home"))
        main_layout.addWidget(self._nav_rail)

        self._explorer = ExplorerPane()
        self._explorer.folder_opened.connect(self._on_folder_opened)

        self._stack = QStackedWidget()
        self._stack.setMinimumHeight(200)

        self._output_panel = OutputPanel()
        self._output_panel.setMinimumHeight(60)
        self._vsplitter = QSplitter(Qt.Orientation.Vertical)
        self._vsplitter.addWidget(self._stack)
        self._vsplitter.addWidget(self._output_panel)
        self._vsplitter.setStretchFactor(0, 3)
        self._vsplitter.setStretchFactor(1, 1)
        self._vsplitter.setChildrenCollapsible(False)
        self._vsplitter.setHandleWidth(6)

        saved_vsizes = self.settings.value("output_splitter")
        if saved_vsizes:
            self._vsplitter.setSizes([int(s) for s in saved_vsizes])
        else:
            self._vsplitter.setSizes([700, 120])

        self._splitter = QSplitter(Qt.Orientation.Horizontal)
        self._splitter.addWidget(self._explorer)
        self._splitter.addWidget(self._vsplitter)
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setCollapsible(0, False)
        self._splitter.setCollapsible(1, False)
        self._splitter.setHandleWidth(2)

        saved_sizes = self.settings.value("explorer_splitter")
        if saved_sizes:
            self._splitter.setSizes([int(s) for s in saved_sizes])
        else:
            self._splitter.setSizes([280, 800])

        main_layout.addWidget(self._splitter)

        # Reserve a stack slot per workspace at its fixed page index.
        # Workspaces are built eagerly below and swap in immediately --
        # users never see these blank slots.
        self._placeholders: dict = {}
        for mode_key, (page_index, *_) in self._MODE_DISPATCH.items():
            slot = QWidget()
            self._placeholders[mode_key] = slot
            self._stack.insertWidget(page_index, slot)

        # Build every workspace eagerly during startup. Each nav-rail
        # click is then just a stack-page swap (no per-click build hang).
        # The splash's set_status callback (if any) shows progress.
        for mode_key, (_, label, _, loader, _) in self._MODE_DISPATCH.items():
            self._splash_status(f"Building {label}...")
            getattr(self, loader)()

        # Land on Project if a project was auto-reopened, else Home.
        if self.current_project_dir is not None:
            self._nav_rail.set_active("project")
            self._activate_project()
            self._output_panel.log("Project")
        else:
            self._stack.setCurrentIndex(PAGE_HOME)
            self._output_panel.log("Home")
            self._explorer.set_workflow_steps([])

        # Auto-restore the last opened project, if it still exists.
        last_project = self.settings.value("last_project_dir", "") or ""
        if last_project and os.path.isdir(str(last_project)):
            self.set_project_directory(str(last_project))

        self._status_bar = self.statusBar()
        self._resource_label = QLabel("CPU: —  RAM: —")
        self._resource_label.setProperty("role", "status_metric")
        self._status_bar.addPermanentWidget(self._resource_label)

        # Prime psutil's CPU counter so the first real tick has a baseline
        # to diff against (otherwise cpu_percent always reads 0).
        if psutil is not None:
            psutil.cpu_percent(interval=None)
        self._resource_timer = QTimer(self)
        self._resource_timer.timeout.connect(self._update_resource_usage)
        self._resource_timer.start(2000)
        self._update_resource_usage()

        saved_root = self.settings.value("explorer_root", "")
        if not saved_root or not os.path.isdir(saved_root):
            saved_root = self.settings.value("last_directory", "")
        if saved_root and os.path.isdir(saved_root):
            self._explorer.set_root(saved_root)
        else:
            self._explorer.set_root(os.path.dirname(os.path.abspath(__file__)))

        # Help panel docks to the right; hidden until F1 / Help menu is hit.
        from kosmic.gui.help.help_browser import HelpManager
        HelpManager.instance().attach_main_window(self)

    def _splash_status(self, message: str) -> None:
        """Forward a startup-progress message to the splash if one is attached."""
        cb = getattr(self, '_splash_status_cb', None)
        if cb is not None:
            cb(message)

    def _build_menu_bar(self):
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("File")

        open_action = QAction("Open Project...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._open_project_from_menu)
        file_menu.addAction(open_action)

        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        view_menu = menu_bar.addMenu("View")

        plot_settings_action = QAction("Plot Settings...", self)
        plot_settings_action.setShortcut("Ctrl+,")
        plot_settings_action.triggered.connect(self._open_plot_settings)
        view_menu.addAction(plot_settings_action)

        ma_settings_action = QAction("Meta-Analysis Settings...", self)
        ma_settings_action.triggered.connect(self._open_ma_settings)
        view_menu.addAction(ma_settings_action)

        from kosmic.gui.help.about_dialogs import show_about
        help_menu = menu_bar.addMenu("Help")

        help_action = QAction("Help on This Page", self)
        help_action.setShortcut("F1")
        help_action.setToolTip("Open the documentation for this screen.")
        help_action.triggered.connect(self._show_help)
        help_menu.addAction(help_action)
        offline_action = QAction("Offline Help Browser", self)
        offline_action.setToolTip(
            "The same pages, bundled with the app, for when there is no "
            "internet connection.")
        offline_action.triggered.connect(self._show_offline_help)
        help_menu.addAction(offline_action)

        help_menu.addSeparator()

        about_action = QAction("About KOSMIC", self)
        about_action.triggered.connect(lambda: show_about(self))
        help_menu.addAction(about_action)

        help_menu.addSeparator()

        tutorial_action = QAction("Start Tutorial", self)
        tutorial_action.triggered.connect(self._start_tutorial)
        help_menu.addAction(tutorial_action)

    def _show_help(self):
        """F1: the documentation page for the screen in front."""
        from kosmic.gui.help.help_browser import HelpManager
        HelpManager.instance().open(self._resolve_help_id())

    def _show_offline_help(self):
        from kosmic.gui.help.help_browser import HelpManager
        HelpManager.instance().show_full(self._resolve_help_id())

    def _resolve_help_id(self) -> str:
        """Find the help_id of the page currently in front.

        Prefers the focused widget, since that pins down which inner tab
        the user is actually in -- but only when focus is still inside the
        current stack page. After switching workspace or tab the focus
        widget is often left behind on the page you came from, which would
        resolve the wrong topic.

        Falls back to the visible page itself, then to the active
        workspace's overview, then to the index.
        """
        from PyQt6.QtWidgets import QApplication

        page = self._stack.currentWidget()

        def _walk(widget):
            while widget is not None:
                help_id = getattr(widget, "help_id", None)
                if help_id:
                    return help_id
                widget = widget.parent()
            return None

        focused = QApplication.focusWidget()
        if focused is not None and page is not None and page.isAncestorOf(focused):
            found = _walk(focused)
            if found:
                return found

        found = _walk(self._current_page_widget())
        if found:
            return found

        # Fall back to a workspace overview based on the current stack page.
        idx = self._stack.currentIndex()
        return {
            PAGE_SCRNA:   "scrna/index",
            PAGE_DE:      "de/index",
            PAGE_FIGURES: "figures/index",
            PAGE_META:    "meta/index",
        }.get(idx, "index")

    def _current_page_widget(self):
        """The innermost visible widget of the page in front.

        Descends through the workspaces' own stacks / tab widgets so a
        per-tab help_id is found rather than the workspace's overview.
        """
        from PyQt6.QtWidgets import QStackedWidget, QTabWidget

        widget = self._stack.currentWidget()
        for _ in range(6):                    # depth guard, not a loop count
            nested = None
            for attr in ('stack', '_stack', 'tabs', '_tabs'):
                candidate = getattr(widget, attr, None)
                if isinstance(candidate, (QStackedWidget, QTabWidget)):
                    nested = candidate.currentWidget()
                    break
            if nested is None:
                return widget
            widget = nested
        return widget

    def _open_plot_settings(self):
        current = self._stack.currentIndex()

        if current == PAGE_META and self._meta_workspace is not None:
            from kosmic.gui.meta_analysis.dialogs.ma_plot_settings import (
                MAPlotSettingsDialog)
            if (not hasattr(self, '_ma_plot_dlg')
                    or self._ma_plot_dlg is None):
                self._ma_plot_dlg = MAPlotSettingsDialog(
                    self._meta_workspace.gene_ma_page, parent=self)
                self._ma_plot_dlg.finished.connect(
                    lambda: setattr(self, '_ma_plot_dlg', None))
            self._ma_plot_dlg.show()
            self._ma_plot_dlg.raise_()
            self._ma_plot_dlg.activateWindow()
            return

        if self._de_workspace is None:
            dialogs.info(self, "Plot Settings",
                         "Load the DE workspace first.")
            return
        from kosmic.gui.de_analysis.plot_settings import PlotSettingsDialog
        if not hasattr(self, '_plot_settings_dlg') or self._plot_settings_dlg is None:
            self._plot_settings_dlg = PlotSettingsDialog(self._de_workspace, parent=self)
            self._plot_settings_dlg.finished.connect(
                lambda: setattr(self, '_plot_settings_dlg', None))
        else:
            self._plot_settings_dlg._auto_select_tab()
        self._plot_settings_dlg.show()
        self._plot_settings_dlg.raise_()
        self._plot_settings_dlg.activateWindow()

    def _open_ma_settings(self):
        from kosmic.gui.meta_analysis.dialogs.ma_settings_dialog import (
            MASettingsDialog, load_de_method)
        dlg = MASettingsDialog(self.settings, parent=self)
        if dlg.exec():
            ws = self._meta_workspace
            if ws is not None:
                ws.gene_ma_page.refresh_method_controls()
                ws.de_method = load_de_method(self.settings)

    # Mode registry: stack page index, log label, "is this workspace
    # built yet" attr name, loader method name, activator method name.
    _MODE_DISPATCH: ClassVar[dict[str, tuple]] = {
        "home":    (PAGE_HOME,    "Home", "_home_loaded",
                    "_load_home_mode", "_activate_home"),
        "project": (PAGE_PROJECT, "Project", "_project_loaded",
                    "_load_project_mode", "_activate_project"),
        "scrna":   (PAGE_SCRNA,   "sc/snRNA-seq Analysis", "_scrna_loaded",
                    "_load_scrna_mode", "_activate_scrna"),
        "de":      (PAGE_DE,      "Differential Expression", "_de_loaded",
                    "_load_de_mode", "_activate_de"),
        "figures": (PAGE_FIGURES, "Figure Export", "_figures_loaded",
                    "_load_figures_mode", "_activate_figures"),
        "meta":    (PAGE_META,    "Meta-Analysis", "_meta_loaded",
                    "_load_meta_mode", "_activate_meta"),
    }

    # Workbenches that require state set on the Project page first.
    _REQUIRES_PROJECT = {"meta"}                         # needs project dir
    _REQUIRES_ACTIVE_STUDY = {"scrna", "de", "figures"}  # needs active study

    def _on_mode_selected(self, mode: str):
        spec = self._MODE_DISPATCH.get(mode)
        if spec is None:
            return

        # Gate workbenches that need project / study state.
        if mode in self._REQUIRES_PROJECT and self.current_project_dir is None:
            dialogs.warning(
                self, "Open a project first",
                "This workspace needs a project folder. Open or create one "
                "from the Project page, then try again.")
            self._nav_rail.set_active("project")
            self._on_mode_selected("project")
            return
        if mode in self._REQUIRES_ACTIVE_STUDY:
            # The study highlighted in the Project table is the one
            # per-study workspaces open on; nobody has to 'set' it.
            # Leaving the Project page adopts the selection (same as
            # its 'Open in...' buttons); switching between other
            # workspaces keeps whatever study they already hold.
            ws = self._project_workspace
            selected = ws.selected_accession if ws is not None else None
            on_project = self._stack.currentIndex() == PAGE_PROJECT
            if selected and (on_project or self.current_study_path is None):
                self.set_active_study(selected)
            elif self.current_study_path is None:
                dialogs.warning(
                    self, "Add a study first",
                    "This workspace works on one study at a time. Add a "
                    "study on the Project page and select it in the table, "
                    "then open this workspace.")
                self._nav_rail.set_active("project")
                self._on_mode_selected("project")
                return

        page_index, label, loaded_attr, loader, activator = spec

        self._nav_rail.set_active(mode)
        with contextlib.suppress(TypeError):
            self._explorer.step_clicked.disconnect()

        # All workspaces are built eagerly during startup, so this is
        # now just a stack-page swap + the workspace's own activation.
        self._stack.setCurrentIndex(page_index)
        getattr(self, activator)()
        self._output_panel.log(label)

    def _activate_scrna(self):
        ws = self._scrna_workspace
        if not self._activate_workspace_common(
                PAGE_SCRNA, ws, ws.WORKFLOW_STEPS,
                on_step_clicked=ws.on_sidebar_step):
            return
        active_tab = ws.stack.currentIndex()
        self._explorer.set_active_step(active_tab)
        self._update_sidebar_status(active_tab)

    def _activate_de(self):
        ws = self._de_workspace

        # Cross-push latest adata / project dir from scRNA.
        if ws and self._scrna_workspace is not None:
            scrna = self._scrna_workspace
            self._sync_scrna_to_de(ws, scrna)
            if scrna.current_project_dir:
                ws.set_project_directory(str(scrna.current_project_dir))

        steps = ws.steps_for_mode(ws.analysis_mode) if ws else []

        if not self._activate_workspace_common(
                PAGE_DE, ws, steps,
                on_step_clicked=ws.on_sidebar_step if ws else None):
            return

        self._sync_de_sidebar_to_stack()
        # Now that the workspace is on screen, let the Dataset page load
        # the study if scRNA did not hand one over. set_project_directory
        # above deliberately does not do this while DE is hidden.
        if ws.current_adata is None:
            ws.switch_tab(0)
        if ws.project_dir:
            self._explorer.update_data_status(
                f"<b>Project:</b> {ws.project_dir.name}")
        else:
            self._explorer.update_data_status("Set a project directory to begin")
        if ws.progress_bar:
            ws.progress_bar.setValue(0)

    def _activate_figures(self):
        # No need for workflow steps in Figure Export.
        ws = self._figures_workspace
        if not self._activate_workspace_common(PAGE_FIGURES, ws, []):
            return
        ws.set_scrna_workspace(self._scrna_workspace)
        ws.set_de_workspace(self._de_workspace)
        ws.on_activated()

    def _activate_meta(self):
        ws = self._meta_workspace
        steps = ws.steps_for_mode(ws.consensus_mode) if ws else []
        if not self._activate_workspace_common(
                PAGE_META, ws, steps,
                on_step_clicked=self._on_meta_step_clicked):
            return
        self._explorer.set_active_step(ws.current_stack_page())
        if ws.project_folder:
            self._explorer.set_root(ws.project_folder)

    def _sync_scrna_to_de(self, ws, scrna) -> None:
        """Push scRNA's adata into DE, gated on readiness and staleness.

        Only pushes once scRNA's data actually passes 'de_readiness'
        (normalised + Control/Disease roles assigned) -- switching to DE
        mid-pipeline shows a clear reason instead of silently accepting
        raw/unlabelled data. Re-pushes only when scRNA's '_adata_version'
        has advanced since the last sync, and never overwrites a dataset
        the user loaded into DE manually (see 'DEWorkspace._manual_load').
        """
        if ws is None or scrna is None or scrna.current_adata is None:
            return
        if getattr(ws, '_manual_load', False):
            return
        if getattr(ws, '_last_scrna_adata_version', None) == scrna._adata_version:
            return

        from kosmic.scrna.inspect.roles import de_readiness
        ready, reason = de_readiness(scrna.current_adata)
        if not ready:
            ws.show_scrna_not_ready(reason)
            return

        if hasattr(scrna, 'current_h5ad_path'):
            ws.h5ad_path = scrna.current_h5ad_path
        ws.set_adata(scrna.current_adata)
        ws._last_scrna_adata_version = scrna._adata_version

    def _sync_de_sidebar_to_stack(self):
        ws = self._de_workspace
        if ws is None:
            return
        stack_page = ws.stack.currentIndex()
        page_map = ws.get_page_map()
        for step_idx, mapped_page in enumerate(page_map):
            if mapped_page == stack_page:
                self._explorer.set_active_step(step_idx)
                return
        self._explorer.set_active_step(0)

    def _update_de_sidebar_status(self, tab_index: int):
        ws = self._de_workspace
        if ws is None:
            return
        text, state = ws.compute_sidebar_status(tab_index)
        self._explorer.update_data_status(text, state)

    def _on_de_mode_changed(self, mode: str):
        ws = self._de_workspace
        if ws is None:
            return

        with contextlib.suppress(TypeError):
            self._explorer.step_clicked.disconnect()

        self._explorer.set_workflow_steps(ws.steps_for_mode(mode))
        self._explorer.step_clicked.connect(ws.on_sidebar_step)

        for i, done in ws.iter_completed_steps():
            if done:
                self._explorer.set_step_complete(i)

        self._sync_de_sidebar_to_stack()

    # One dataset in memory at a time. A study is tens of gigabytes; two
    # workspaces each holding a different one is what emptied the machine
    # when the atlas was opened for its DE while scRNA still held a study.
    def _on_de_dataset_loaded(self, path: str):
        # Same study or not: DE reads its own counts-only copy from disk,
        # so whatever scRNA holds is redundant from here on. (The same
        # file was the case that reached 48 GB: the atlas in scRNA plus
        # its counts in DE.)
        scrna = self._scrna_workspace
        if scrna is None or scrna.current_adata is None:
            return
        held = scrna.current_h5ad_path or "its dataset"
        scrna.release_dataset(
            f"Released {Path(held).name} from memory: DE is loading {Path(path).name}.")
        self._output_panel.log(
            f"Released {Path(held).name} from the scRNA workspace: DE is loading "
            f"{Path(path).name}. One dataset is kept in memory at a time.")

    def _on_scrna_dataset_loaded(self, path: str):
        de = self._de_workspace
        if de is None or de.current_adata is None or not getattr(de, '_manual_load', False):
            return
        held = getattr(de, 'h5ad_path', None) or "its dataset"
        de.release_dataset()
        self._output_panel.log(
            f"Released {Path(held).name} from the DE workspace: scRNA loaded "
            f"{Path(path).name}. One dataset is kept in memory at a time.")

    def _on_de_open_scrna_requested(self):
        """DE's Setup page asked to jump to scRNA -> Inspect (to assign
        condition roles, or to finish an earlier pipeline step)."""
        self._on_mode_selected("scrna")
        if self._scrna_workspace is not None:
            self._scrna_workspace.switch_tab(2)  # Inspect tab

    def _on_meta_step_clicked(self, index: int):
        if self._meta_workspace is not None:
            self._meta_workspace.switch_page(index)
            self._explorer.set_active_step(index)

    def _on_meta_mode_changed(self, mode: str):
        ws = self._meta_workspace
        if ws is None:
            return
        self._explorer.set_workflow_steps(ws.steps_for_mode(mode))
        # Re-mark completed steps so ticks survive the rebuild.
        for i, done in ws.iter_completed_steps():
            if done:
                self._explorer.set_step_complete(i)
        # Mode selection always lands on sidebar index 3 (the first
        # post-mode step). The workspace's switch_page(3) fires after
        # this signal, so we can't reverse-lookup yet.
        self._explorer.set_active_step(3)

    def _connect_workspace_signals(
        self,
        workspace,
        *,
        wire_log_message: bool = True,
        wire_project_dir: bool = True,
        wire_step_completed: bool = True,
        progress_bar=None,
    ) -> None:
        """Wire workspace signals (status_message, log_message,
        project_directory_changed, step_completed) to output-panel +
        explorer slots. Optional kwargs let a workspace skip a signal
        it doesn't emit."""
        workspace.status_message.connect(self._output_panel.log)
        if wire_log_message and hasattr(workspace, 'log_message'):
            workspace.log_message.connect(self._output_panel.log)
        if (wire_project_dir
                and hasattr(workspace, 'project_directory_changed')):
            workspace.project_directory_changed.connect(
                self._update_explorer_root)
        if wire_step_completed and hasattr(workspace, 'step_completed'):
            workspace.step_completed.connect(
                self._explorer.set_step_complete)
        if progress_bar is not None:
            workspace.progress_bar = progress_bar

    def _swap_placeholder_into_stack(self, mode_key: str, widget) -> None:
        page_index = self._MODE_DISPATCH[mode_key][0]
        placeholder = self._placeholders.pop(mode_key)
        self._stack.removeWidget(placeholder)
        placeholder.deleteLater()
        self._stack.insertWidget(page_index, widget)

    @staticmethod
    def _build_error_page(message: str) -> QWidget:
        """Centred, word-wrapped fallback for failed workspace imports."""
        page = QWidget()
        layout = QVBoxLayout(page)
        label = QLabel(message)
        label.setWordWrap(True)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(label)
        return page

    def _activate_workspace_common(
        self,
        page_index: int,
        workspace,
        steps: list,
        on_step_clicked=None,
    ) -> bool:
        """Switch to ``page_index``, install ``steps`` in the explorer,
        wire ``on_step_clicked``, and replay completed-step state.
        Returns False (and the caller bails) if the workspace isn't
        built yet."""
        self._stack.setCurrentIndex(page_index)
        self._explorer.set_workflow_steps(steps)
        if on_step_clicked is not None:
            self._explorer.step_clicked.connect(on_step_clicked)
        if workspace is None:
            return False
        # Step-less workspaces (e.g. Figure Export) have no completed-step
        # state to replay; skip rather than crash the activation.
        if hasattr(workspace, 'iter_completed_steps'):
            for i, done in workspace.iter_completed_steps():
                if done:
                    self._explorer.set_step_complete(i)
        return True

    def _load_scrna_mode(self):
        from kosmic.gui.scrna.workspace import ScRNAWorkspace

        workspace = ScRNAWorkspace(embedded=True)
        self._scrna_workspace = workspace
        self._connect_workspace_signals(workspace, wire_log_message=False)
        workspace.step_completed.connect(self._refresh_sidebar_status)
        workspace.steps_reset.connect(self._on_steps_reset)
        workspace.tab_changed.connect(self._on_scrna_tab_changed)
        workspace.dataset_loaded.connect(self._on_scrna_dataset_loaded)

        workspace.download_tab.data_status_changed.connect(
            self._output_panel.log)
        panel = self._explorer.progress_panel
        for tab in (workspace.download_tab, workspace.gene_names_tab,
                    workspace.qc_tab, workspace.cluster_tab,
                    workspace.annotate_tab, workspace.filter_tab):
            tab.progress_bar = panel.progress_bar
            tab.log_message.connect(self._output_panel.log)
        workspace.download_tab.set_status_label(panel.progress_label)
        workspace.inspect_tab.log_message.connect(self._output_panel.log)

        self._swap_placeholder_into_stack("scrna", workspace)
        self._scrna_loaded = True

        # Project workspace already loaded and may have auto-reopened a
        # project; catch this workspace up to the active study.
        if self.current_study_path is not None:
            self._sync_scrna_to_active_study(self.current_study_path)
        elif self.current_project_dir is not None:
            workspace.set_project_directory(str(self.current_project_dir))

    def _load_de_mode(self):
        try:
            from kosmic.gui.de_analysis.workspace import DEWorkspace
            workspace = DEWorkspace()
            self._connect_workspace_signals(
                workspace, progress_bar=self._explorer.progress_panel.progress_bar)
            workspace.tab_changed.connect(self._explorer.set_active_step)
            workspace.tab_changed.connect(self._update_de_sidebar_status)
            workspace.step_completed.connect(
                lambda _: self._update_de_sidebar_status(0))
            workspace.mode_changed.connect(self._on_de_mode_changed)
            workspace.study_change_requested.connect(self.set_active_study)
            workspace.open_scrna_requested.connect(self._on_de_open_scrna_requested)
            workspace.dataset_loading.connect(self._on_de_dataset_loaded)
            workspace.dataset_loaded.connect(self._on_de_dataset_loaded)
            self._de_workspace = workspace

            # Catch up if a study was activated before DE lazy-loaded.
            # The AppWindow contract is the source of truth; the cross-
            # push from scRNA only adds the in-memory adata if scRNA had
            # already loaded one and it's ready to hand off (see
            # '_sync_scrna_to_de').
            if self.current_study_path is not None:
                workspace.set_project_directory(str(self.current_study_path))
            self._sync_scrna_to_de(workspace, self._scrna_workspace)

        except Exception as e:
            workspace = self._build_error_page(
                f"Differential Expression module could not be loaded:\n{e}")

        self._swap_placeholder_into_stack("de", workspace)
        self._de_loaded = True

    def _load_figures_mode(self):
        try:
            from kosmic.gui.figure_export import FigureExportWorkspace
            workspace = FigureExportWorkspace(
                scrna_workspace=self._scrna_workspace,
                de_workspace=self._de_workspace,
            )
            # No workflow steps; reads project dir from upstream workspaces.
            self._connect_workspace_signals(
                workspace, wire_project_dir=False,
                wire_step_completed=False)
            self._figures_workspace = workspace
        except Exception as e:
            workspace = self._build_error_page(
                f"Figure Export module could not be loaded:\n{e}")

        self._swap_placeholder_into_stack("figures", workspace)
        self._figures_loaded = True

    def _open_combine_dialog(self) -> None:
        """Open the Combine modal from the Project workspace."""
        from kosmic.gui.combine import CombineDialog
        dlg = CombineDialog(self, parent=self)
        dlg.log_message.connect(self._output_panel.log)
        dlg.master_built.connect(self._on_master_built)
        dlg.open_master_requested.connect(self._on_open_master_in_scrna)
        dlg.exec()

    @staticmethod
    def _master_study_names(master_path) -> set:
        """Accessions whose cells are actually in the master.

        Read from the master's own 'study' column, which concat_studies
        writes. Returns an empty set when it cannot be read, so the
        caller falls back to offering every study rather than refusing.
        """
        try:
            from kosmic.scrna.load.h5ad_meta import read_obs
            obs = read_obs(master_path)
            if 'study' not in obs.columns:
                return set()
            return {str(v) for v in obs['study'].unique()}
        except Exception:
            return set()

    def _open_propagate_dialog(self) -> None:
        """Run label propagation from the master to each per-study h5ad."""
        from kosmic.gui.combine import _PropagateWorker
        from kosmic.paths import master_h5ad_path, processed_data_dir, list_studies, MASTER_ACCESSION

        if self.current_project_dir is None:
            return
        master = master_h5ad_path(self.current_project_dir)
        if not master.exists():
            dialogs.warning(self, "No master found",
                            "Build the master first via Combine studies...")
            return

        # Only the studies actually in the master. It used to offer every
        # study in the project, so a five-study folder with a two-study
        # master proposed writing those labels onto three studies whose
        # cells were never in it -- which, with no barcodes in common,
        # falls through to sc.tl.ingest and projects labels by kNN onto a
        # cohort the master has never seen.
        members = self._master_study_names(master)
        study_paths = []
        for accession in list_studies(self.current_project_dir):
            if accession == MASTER_ACCESSION:
                continue
            if members and accession not in members:
                continue
            pdir = processed_data_dir(self.current_project_dir / accession)
            if not pdir.is_dir():
                continue
            h5ads = sorted(pdir.glob("*.h5ad"),
                           key=lambda p: p.stat().st_mtime, reverse=True)
            if h5ads:
                study_paths.append(h5ads[0])

        if not study_paths:
            dialogs.warning(
                self, "No studies",
                "No per-study h5ads found for the studies in this master."
                if members else "No per-study h5ads found.")
            return

        if not dialogs.confirm(
                self, "Propagate labels",
                "Write the master's cell types onto "
                + ", ".join(sorted(p.parent.parent.name for p in study_paths))
                + "?" + chr(10) + chr(10)
                + "They land in new columns -- 'cell_type_atlas' and "
                "'leiden_atlas'. Each study keeps its own 'cell_type' and "
                "'leiden' from when it was processed on its own, so you can "
                "compare the two." + chr(10) + chr(10)
                + "The study files are rewritten to add the columns."):
            return

        self._output_panel.log("Propagating master labels to studies...")
        worker = _PropagateWorker(
            master, study_paths,
            label_cols=('leiden', 'cell_type'),
            parent=self,
        )
        from kosmic.gui.shared import run_worker as run_worker_helper
        # Hold a ref so Qt doesn't GC it mid-run.
        self._propagate_worker = worker
        run_worker_helper(
            worker,
            on_finished=lambda res: self._output_panel.log(
                f"Propagated to {len(res)} studies "
                f"({sum(n for _, n in res):,} cells)."),
            on_failed=lambda msg: self._output_panel.log(
                f"Propagation failed: {msg}"),
            on_progress=self._output_panel.log,
        )

    def _on_master_built(self, master_path: str) -> None:
        """A new master.h5ad just landed -- refresh the project study table."""
        self.study_list_changed.emit()
        self._output_panel.log(f"Master ready: {master_path}")

    def _on_open_master_in_scrna(self) -> None:
        """Set _master as the active study and jump to the scRNA workspace."""
        from kosmic.paths import MASTER_ACCESSION
        self.set_active_study(MASTER_ACCESSION)
        self._on_mode_selected("scrna")

    def _load_meta_mode(self):
        try:
            from kosmic.gui.meta_analysis.workspace import MetaAnalysisWorkspace
            workspace = MetaAnalysisWorkspace()
            self._connect_workspace_signals(
                workspace, progress_bar=self._explorer.progress_panel.progress_bar)
            workspace.data_status.connect(self._explorer.update_data_status)
            workspace.mode_changed.connect(self._on_meta_mode_changed)
            workspace.tab_changed.connect(self._explorer.set_active_step)
            self._meta_workspace = workspace
            self.project_changed.connect(self._sync_meta_to_active_project)
            # Pull the current project directory if one is already loaded.
            if self.current_project_dir is not None:
                workspace.set_project_directory_external(str(self.current_project_dir))
            page = QWidget()
            borderless(QVBoxLayout, page).addWidget(workspace)
        except Exception as e:
            page = self._build_error_page(
                f"Meta-Analysis module could not be loaded:\n{e}")

        self._swap_placeholder_into_stack("meta", page)
        self._meta_loaded = True

    def _sync_meta_to_active_project(self, project_path):
        """Retarget the Meta workspace whenever the active project changes."""
        if project_path is None or self._meta_workspace is None:
            return
        self._meta_workspace.set_project_directory_external(str(project_path))

    def _load_home_mode(self):
        try:
            from kosmic.gui.home import HomeWorkspace
            workspace = HomeWorkspace()
            workspace.tutorial_clicked.connect(self._start_tutorial)
            workspace.help_clicked.connect(self._open_help)
            self._home_workspace = workspace
        except Exception as e:
            workspace = self._build_error_page(
                f"Home workspace could not be loaded:\n{e}")
        self._swap_placeholder_into_stack("home", workspace)
        self._home_loaded = True

    def _activate_home(self):
        self._stack.setCurrentIndex(PAGE_HOME)
        self._explorer.set_workflow_steps([])

    def _open_help(self):
        """Slot for the Home page's 'Open help' link."""
        from kosmic.gui.help.help_browser import HelpManager
        HelpManager.instance().open()

    def _load_project_mode(self):
        try:
            from kosmic.gui.project import ProjectWorkspace
            last = self.settings.value("last_directory", "") or ""
            workspace = ProjectWorkspace(last_directory=str(last))
            workspace.open_requested.connect(self._on_project_open_requested)
            workspace.create_requested.connect(self._on_project_create_requested)
            workspace.recent_project_clicked.connect(self._on_project_open_requested)
            workspace.study_selected.connect(self.set_active_study)
            workspace.studies_add_requested.connect(self._on_studies_add_requested)
            workspace.download_requested.connect(self._on_project_download_requested)
            workspace.refresh_requested.connect(self._refresh_project_workspace)
            workspace.combine_requested.connect(self._open_combine_dialog)
            workspace.study_action_requested.connect(self._on_study_action)
            workspace.open_in_requested.connect(self._on_open_study_in)
            workspace.import_requested.connect(self._on_project_import_requested)
            workspace.step_changed.connect(self._on_project_step_changed)
            workspace.progress_changed.connect(self._sync_project_explorer)
            self.project_changed.connect(self._refresh_project_workspace)
            self.study_list_changed.connect(self._refresh_project_workspace)
            # Study selection only needs to update the active-row marker,
            # not re-scan the filesystem to recompute every status badge.
            self.study_changed.connect(self._on_active_study_changed_project)
            self.study_changed.connect(self._sync_scrna_to_active_study)
            self.study_changed.connect(self._sync_de_to_active_study)
            self._project_workspace = workspace
            workspace.set_recent_projects(self._load_recent_projects())
            self._refresh_project_workspace()
            # Auto-reopen the last project on first load (subsequent
            # mode switches see current_project_dir set, no-op).
            if self.current_project_dir is None:
                last_project = self.settings.value("last_project_dir", "") or ""
                if last_project and Path(last_project).is_dir():
                    self._on_project_open_requested(last_project)
        except Exception as e:
            workspace = self._build_error_page(
                f"Project workspace could not be loaded:\n{e}")

        self._swap_placeholder_into_stack("project", workspace)
        self._project_loaded = True

    def _activate_project(self):
        ws = self._project_workspace
        steps = ws.WORKFLOW_STEPS if ws is not None else []
        if not self._activate_workspace_common(
                PAGE_PROJECT, ws, steps,
                on_step_clicked=ws.on_sidebar_step if ws else None):
            return
        self._explorer.set_active_step(ws.current_step)
        self._explorer.update_data_status(*ws.compute_sidebar_status())
        if self.current_project_dir is not None:
            ws.set_last_directory(str(self.current_project_dir.parent))

    def _on_project_step_changed(self, index: int) -> None:
        if self._stack.currentIndex() == PAGE_PROJECT:
            self._explorer.set_active_step(index)

    def _sync_project_explorer(self) -> None:
        """Replay step completion + STATUS after the project page
        re-scans the folder. Only while the Project page is showing;
        other workspaces own the explorer the rest of the time."""
        ws = self._project_workspace
        if ws is None or self._stack.currentIndex() != PAGE_PROJECT:
            return
        for i, done in ws.iter_completed_steps():
            self._explorer.set_step_complete(i, done)
        self._explorer.set_active_step(ws.current_step)
        self._explorer.update_data_status(*ws.compute_sidebar_status())

    def _refresh_project_workspace(self, *args):
        ws = self._project_workspace
        if ws is None:
            return
        if self.current_project_dir is None:
            ws.show_empty()
        else:
            ws.show_project(self.current_project_dir, self.current_accession)

    def _on_active_study_changed_project(self, *args):
        """Update the Project tab when the active study toggles.

        Cheap update path: only the row markers (bold name + Active
        badge) change. Full re-renders still happen on project /
        study-list changes.
        """
        ws = self._project_workspace
        if ws is None:
            return
        ws.update_active_study(self.current_accession)
        # Confirmation feedback so the user sees something happened.
        if self.current_accession:
            msg = f"Active study: {self.current_accession}"
        else:
            msg = "No active study."
        self._output_panel.log(msg)
        self.statusBar().showMessage(msg, 4000)

    def _on_project_download_requested(self):
        """Jump to the scRNA Load Data tab so the user can download into the active study."""
        self._on_mode_selected("scrna")
        if self._scrna_workspace is not None:
            self._scrna_workspace.switch_tab(0)  # Load Data is tab 0

    def _on_project_import_requested(self, accession: str, kind: str) -> None:
        """Project 'Add Data...' menu: run the intake route for the study
        selected in the table (ADR-001 Phase E). The study becomes active
        first, since the conversion that follows and the scRNA sync both
        act on the active study."""
        if accession:
            self.set_active_study(accession)
        if self.current_study_path is None:
            dialogs.warning(
                self, "Pick a study first",
                "Data is registered to a study. Select one in the study "
                "table, then import.")
            return

        if kind in ('geo', 'assemble'):
            from kosmic.gui.intake.dialog import AddDataDialog
            dlg = AddDataDialog(
                self.current_study_path,
                log_cb=self._output_panel.log,
                on_data_changed=lambda: None,  # refreshed after close
                parent=self,
                start=kind,
            )
            dlg.exec()
            if dlg.data_changed:
                self._refresh_project_workspace()
                # If scRNA is loaded on this study, let its Dataset
                # screen re-detect (shows the Convert banner for raw).
                scrna = self._scrna_workspace
                if scrna is not None:
                    scrna.download_tab.on_tab_activated()
            return

        # Local file imports complete on the Project page: pick, copy into
        # raw_data/, then convert straight into processed_data/ so the study
        # is analysis-ready without a detour through another workspace.
        from kosmic.gui.intake.local_import import import_local_file
        copied = import_local_file(
            self, self.current_study_path, kind, self._output_panel.log)
        if copied is None:
            return
        self._refresh_project_workspace()
        self._start_study_import(copied, kind)

    def _choose_counts_matrix(self, raw_path):
        """Return '(slot, described)' for an h5ad, asking when it is a real choice.

        Returns 'slot=None' to let the worker pick (one obvious candidate),
        a slot token when the user chose one, or 'False' if they cancelled.
        Inspection loads no matrices, so this is instant even on a 13 GB file.
        """
        from PyQt6.QtWidgets import QInputDialog

        from kosmic.scrna.load.matrices import (
            RAW_COUNTS, describe_count_matrices_path,
        )

        try:
            described = describe_count_matrices_path(raw_path)
        except Exception as e:
            self._output_panel.log(f"(could not inspect matrices: {e})")
            return None, []

        counts = [d for d in described if d['verdict'] == RAW_COUNTS]
        if len(counts) < 2:
            return None, described        # nothing to decide; worker picks

        # Say what actually differs rather than guessing at intent: the
        # matrix with fewer stored values has had something removed, which
        # is usually ambient correction or extra filtering. Which is which
        # is the depositor's business, so point at their documentation.
        fullest = max(counts, key=lambda d: d['n_stored'])
        lines = [f"  {d['label']:26s} {d['n_stored']:>15,} stored values, "
                 f"max {d['max']:,.0f}" for d in counts]
        deltas = [
            f"{d['label']} holds {100 * (fullest['n_stored'] - d['n_stored']) / fullest['n_stored']:.1f}%"
            f" fewer values than {fullest['label']}"
            for d in counts
            if d is not fullest and d['n_stored'] < fullest['n_stored']
        ]
        detail = ("\n\n" + "; ".join(deltas)
                  + " — so it has had something removed, typically "
                    "ambient-RNA correction or extra filtering. Which is "
                    "which is the depositor's choice: check the dataset's "
                    "README or the paper's methods.") if deltas else ""

        items = [f"{d['label']}  —  max {d['max']:,.0f}" for d in counts]
        default = next((i for i, d in enumerate(counts) if d['slot'] == 'X'), 0)
        choice, ok = QInputDialog.getItem(
            self, "Which counts matrix?",
            f"{Path(raw_path).name} holds more than one matrix of raw "
            f"integer counts.\n\n" + "\n".join(lines) + detail
            + "\n\nPick the matrix to analyse from. The choice, and the "
              "options rejected, are recorded in the study's provenance.",
            items, default, False)
        if not ok:
            self._output_panel.log("Import cancelled at matrix selection.")
            return False, described
        return counts[items.index(choice)]['slot'], described

    def _start_study_import(self, raw_path, kind: str) -> None:
        """Convert a just-imported raw file into the study's processed_data/.

        Runs on the Project page so setup finishes where the user is, rather
        than deferring to the scRNA Dataset screen. For h5ad the worker also
        picks the raw-counts matrix -- depositors disagree about which slot
        holds it, and a normalised matrix reaches DE without erroring.
        """
        from kosmic.gui.intake.workers import (
            CSVConvertWorker, H5adImportWorker, RDSConvertWorker,
        )
        from kosmic.gui.shared import run_worker
        from kosmic.paths import processed_data_dir

        worker_cls = {'h5ad': H5adImportWorker,
                      'rds': RDSConvertWorker,
                      'csv': CSVConvertWorker}.get(kind)
        if worker_cls is None:
            return

        processed = processed_data_dir(self.current_study_path)
        processed.mkdir(parents=True, exist_ok=True)
        study_name = self.current_study_path.name
        study_dir = self.current_study_path

        # An h5ad can hold more than one valid counts matrix (a CellBender
        # -corrected X alongside the uncorrected CellRanger counts, say).
        # That is a scientific choice, not a technical one, so ask -- and
        # record the answer. Inspection reads no matrices, so it is instant.
        slot = None
        described: list = []
        if kind == 'h5ad':
            slot, described = self._choose_counts_matrix(raw_path)
            if slot is False:            # user cancelled the choice
                return

        def _done(payload):
            _out, message = payload
            self._output_panel.log(message)
            try:
                from kosmic.provenance import record_stage
                record_stage(
                    study_dir, study_name, 'import',
                    params={
                        'source_file': Path(raw_path).name,
                        'matrix_slot': slot or 'X',
                        'candidates': [
                            {'slot': d['slot'], 'verdict': d['verdict'],
                             'max': round(d['max'], 3)} for d in described
                        ],
                    },
                    source=str(raw_path),
                )
            except Exception as e:  # provenance must never block an import
                self._output_panel.log(f"(provenance not recorded: {e})")
            self._refresh_project_workspace()
            scrna = self._scrna_workspace
            if scrna is not None:
                scrna.download_tab.on_tab_activated()
            dialogs.info(
                self, "Import Complete",
                f"{message}\n\nRegistered to '{study_name}'. Open it in "
                "scRNA to set sample/condition columns and run QC.")

        def _failed(message: str):
            self._output_panel.log(f"Import failed: {message}")
            dialogs.warning(self, "Import Failed", message)

        if kind == 'h5ad':
            self._study_import_worker = worker_cls(
                str(raw_path), str(processed), slot=slot)
        else:
            self._study_import_worker = worker_cls(str(raw_path), str(processed))
        run_worker(
            self._study_import_worker,
            on_finished=_done,
            on_failed=_failed,
            on_progress=self._output_panel.log,
        )
        self._output_panel.log(
            f"Importing {raw_path.name} into {study_name}...")

    def _on_open_study_in(self, accession: str, mode: str) -> None:
        """Project details panel: activate a study and open it in a workspace."""
        if not accession:
            return
        if mode in ("scrna", "de"):
            self.set_active_study(accession)
        self._nav_rail.set_active(mode)
        self._on_mode_selected(mode)

    def _on_study_action(self, accession: str, action: str) -> None:
        """Handle a per-row action button or kebab-menu pick from the project page."""
        if self.current_project_dir is None:
            return
        study_path = self.current_project_dir / accession

        if action == "set_active":
            self.set_active_study(accession)
            return

        if action == "download_raw":
            self.set_active_study(accession)
            self._on_project_download_requested()
            return

        if action == "process":
            self.set_active_study(accession)
            self._on_mode_selected("scrna")
            if self._scrna_workspace is not None:
                self._scrna_workspace.switch_tab(2)  # Inspect is tab 2
            return

        if action == "reprocess":
            self.set_active_study(accession)
            self._on_mode_selected("scrna")
            if self._scrna_workspace is not None:
                self._scrna_workspace.switch_tab(3)  # QC is tab 3
            return

        if action == "open_folder":
            self._reveal_in_finder(study_path)
            return

        if action == "propagate_labels":
            self._open_propagate_dialog()
            return

        if action == "remove":
            if not dialogs.confirm(
                self, f"Remove '{accession}'",
                f"Delete the study folder and all its contents?\n\n"
                f"{study_path}\n\nThis cannot be undone."
            ):
                return
            try:
                import shutil
                shutil.rmtree(study_path)
            except Exception as exc:
                dialogs.warning(self, "Remove Failed",
                                f"Could not remove '{accession}':\n{exc}")
                return
            if self.current_accession == accession:
                self.set_active_study(None)
            self.study_list_changed.emit()

    def _reveal_in_finder(self, path) -> None:
        """Open the OS file manager at ``path``."""
        import subprocess
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(path)])
        elif sys.platform == "win32":
            subprocess.Popen(["explorer", str(path)])
        else:
            subprocess.Popen(["xdg-open", str(path)])

    def _on_studies_add_requested(self, accessions):
        """Scaffold a batch of accession folders under the current project."""
        if self.current_project_dir is None:
            return
        added: list[str] = []
        skipped: list[str] = []
        for accession in accessions:
            target = self.current_project_dir / accession
            already_existed = target.is_dir()
            try:
                self.add_study(accession)
            except Exception as exc:
                self._output_panel.log(f"Skipped '{accession}': {exc}")
                skipped.append(accession)
                continue
            if already_existed:
                skipped.append(accession)
            else:
                added.append(accession)
        if added:
            self._output_panel.log(
                f"Scaffolded {len(added)} new "
                f"{'study' if len(added) == 1 else 'studies'}: {', '.join(added)}"
            )
        if skipped:
            self._output_panel.log(
                f"Already existed (left as-is): {', '.join(skipped)}"
            )

    def _sync_scrna_to_active_study(self, study_path):
        """Retarget the scRNA workspace whenever the active study changes.

        The cascade through ``ScRNAWorkspace.set_project_directory`` ends
        in ``download_tab._run_auto_detect``, which auto-loads the most
        recent h5ad via its own background ``_H5adLoadWorker``. No
        AppWindow-level loader -- a second thread would just double the
        IO contention on the same file.
        """
        if study_path is None or self._scrna_workspace is None:
            return
        self._scrna_workspace.set_project_directory(str(study_path))

    def _sync_de_to_active_study(self, study_path):
        """Retarget the DE workspace whenever the active study changes.

        Figure Export reads ``de.project_dir`` and ``scrna.current_project_dir``
        in its page ``on_activated`` hooks, so syncing the upstream
        workspaces keeps it in sync automatically.
        """
        if study_path is None or self._de_workspace is None:
            return
        self._de_workspace.set_project_directory(str(study_path))

    def _on_project_open_requested(self, folder: str):
        self.set_project_directory(folder)
        self.settings.setValue("last_directory", folder)
        self._explorer.set_root(folder)
        self._output_panel.log(f"Project: {folder}")

    def _on_project_create_requested(self, folder: str):
        self.set_project_directory(folder)
        self.settings.setValue("last_directory", str(Path(folder).parent))
        self._explorer.set_root(folder)
        self._output_panel.log(f"Created project: {folder}")

    def _update_explorer_root(self, directory: str):
        self._explorer.set_root(directory)
        self.settings.setValue("explorer_root", directory)

    def _on_folder_opened(self, folder: str):
        self.settings.setValue("explorer_root", folder)
        self._output_panel.log(f"Project: {folder}")

        current_page = self._stack.currentIndex()
        if current_page == PAGE_SCRNA and self._scrna_workspace is not None:
            self._scrna_workspace.set_project_directory(folder)
            self._output_panel.log(f"Analysis folder: {folder}")
        elif current_page == PAGE_DE and self._de_workspace is not None:
            self._de_workspace.set_project_directory(folder)
            self._output_panel.log(f"DE project folder: {folder}")
        elif current_page == PAGE_META and self._meta_workspace is not None:
            self._meta_workspace.set_project_directory_external(folder)
            self._output_panel.log(f"Meta-analysis project: {folder}")

    def _refresh_sidebar_status(self, _step_index=None):
        if self._scrna_workspace:
            active_tab = self._scrna_workspace.stack.currentIndex()
            self._update_sidebar_status(active_tab)

    def _on_steps_reset(self):
        # Clear every step -- raw reload restarts the pipeline from the top.
        n_steps = len(self._scrna_workspace.WORKFLOW_STEPS)
        for i in range(n_steps):
            self._explorer.set_step_complete(i, False)

    def _on_scrna_tab_changed(self, tab_index: int):
        self._explorer.set_active_step(tab_index)
        self._update_sidebar_status(tab_index)

    def _update_sidebar_status(self, tab_index: int):
        ws = self._scrna_workspace
        if ws is None:
            return
        text, state = ws.compute_sidebar_status(tab_index)
        self._explorer.update_data_status(text, state)

    def _open_project_from_menu(self):
        if self._stack.currentIndex() == PAGE_SCRNA and self._scrna_workspace is not None:
            self._scrna_workspace.open_project_dialog()
        elif self._stack.currentIndex() == PAGE_DE and self._de_workspace is not None:
            self._de_workspace.open_project_dialog()
        elif self._stack.currentIndex() == PAGE_META and self._meta_workspace is not None:
            folder = QFileDialog.getExistingDirectory(
                self, "Open Project Folder",
                self.settings.value("last_directory", ""),
            )
            if folder:
                self._explorer.set_root(folder)
                self._meta_workspace.set_project_directory_external(folder)
                self.settings.setValue("explorer_root", folder)

    def _toggle_theme(self):
        self._is_dark = not self._is_dark
        if self._is_dark:
            apply_theme(QApplication.instance(), ThemeMode.DARK)
        else:
            apply_theme(QApplication.instance(), ThemeMode.LIGHT)
        self._nav_rail.refresh_theme()
        self._explorer.refresh_theme()
        self._output_panel.refresh_theme()
        from kosmic.gui.help.help_browser import HelpManager
        HelpManager.instance().refresh_theme()
        if self._scrna_workspace is not None:
            ws = self._scrna_workspace
            ws.refresh_theme()
            for tab in (ws.qc_tab, ws.cluster_tab, ws.gene_names_tab,
                        ws.inspect_tab, ws.annotate_tab):
                tab.refresh_theme()
        if self._de_workspace is not None:
            self._de_workspace.refresh_theme()
        if self._figures_workspace is not None:
            self._figures_workspace.refresh_theme()
        if self._meta_workspace is not None:
            self._meta_workspace.refresh_theme()

    def keyPressEvent(self, event):
        # Ctrl+1..6 switch scRNA tabs when that workspace is active.
        # TODO: Do we want this active on every page? Probably better to remove
        if event.modifiers() == Qt.KeyboardModifier.ControlModifier:
            tab_keys = {
                Qt.Key.Key_1: 0, Qt.Key.Key_2: 1, Qt.Key.Key_3: 2,
                Qt.Key.Key_4: 3, Qt.Key.Key_5: 4, Qt.Key.Key_6: 5,
            }
            if event.key() in tab_keys and self._scrna_workspace is not None:
                if self._stack.currentIndex() == PAGE_SCRNA:
                    tab_index = tab_keys[event.key()]
                    self._scrna_workspace.switch_tab(tab_index)
                    return

        super().keyPressEvent(event)

    def _start_tutorial(self):
        if not self._scrna_loaded:
            self._on_mode_selected("scrna")
        if self._scrna_workspace is not None:
            self._stack.setCurrentIndex(PAGE_SCRNA)
            self._scrna_workspace.start_tutorial()

    def _load_geometry(self):
        geometry = self.settings.value("geometry")
        if not geometry:
            return
        self.restoreGeometry(geometry)
        # Snap back to the primary screen.
        if QApplication.screenAt(self.geometry().center()) is None:
            avail = QApplication.primaryScreen().availableGeometry()
            self.resize(min(self.width(), avail.width()),
                        min(self.height(), avail.height()))
            self.move(avail.center() - self.rect().center())

    def _update_resource_usage(self):
        if psutil is None:
            self._resource_label.setText("CPU/RAM: psutil not installed")
            return
        # System-wide CPU usage
        cpu = psutil.cpu_percent(interval=None)
        vm = psutil.virtual_memory()
        used_gb = (vm.total - vm.available) / (1024 ** 3)
        total_gb = vm.total / (1024 ** 3)
        self._resource_label.setText(
            f"CPU: {cpu:4.1f}%  RAM: {used_gb:.1f}/{total_gb:.0f} GB"
        )

    def closeEvent(self, event):
        self._resource_timer.stop()
        self.settings.setValue("geometry", self.saveGeometry())
        self.settings.setValue("explorer_splitter", self._splitter.sizes())
        self.settings.setValue("output_splitter", self._vsplitter.sizes())
        event.accept()
        # Closing the main window is closing KOSMIC. Qt only quits when
        # the last window closes, and a message box ("Annotations saved
        # to ...") or a dialog left open behind the window is a window,
        # so the process stayed alive, invisible, holding the dataset.
        app = QApplication.instance()
        if app is not None:
            app.closeAllWindows()
            QTimer.singleShot(0, app.quit)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    apply_theme(app, ThemeMode.DARK)
    app.setApplicationName("KOSMIC")
    app.setApplicationVersion("0.1.0")
    app.setWindowIcon(make_icon("dna", "#007acc", 64))

    splash = SplashScreen()
    splash.show()
    app.processEvents()

    def status_cb(msg: str) -> None:
        splash.set_status(msg)
        app.processEvents()

    window = AppWindow(splash_status_cb=status_cb)

    splash.set_status("Ready")
    app.processEvents()
    splash.dismiss_when_ready(window.show)

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
