# Help. F1 opens the documentation website at the page for the screen in
# front; when no site is configured (or it cannot be reached) the same
# page opens in the offline HelpBrowser dialog below. Both read the
# Markdown under 'content/' -- the website is built from that folder by
# mkdocs.yml -- keyed by each page's 'help_id' ('de/gene_de' ->
# 'content/de/gene_de.md' and '<site>/de/gene_de/').


from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QDesktopServices, QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QDialog, QLineEdit, QSplitter, QTextBrowser,
    QToolBar, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from kosmic import HELP_SITE_URL

CONTENT_DIR = Path(__file__).parent / "content"
# The documentation site ('[help] site_url' in config.toml, trailing
# slash). None means 'no site yet': F1 opens the offline browser.
ONLINE_BASE: Optional[str] = HELP_SITE_URL
# Where a source checkout serves the docs to itself when the site is not
# reachable (see HelpManager._start_local_docs).
LOCAL_DOCS_URL = "http://127.0.0.1:8001/"


# Content registry

@dataclass
class HelpEntry:
    """One page in the help tree."""
    help_id: str          # e.g. "de/gene_de"
    title: str            # first H1 of the markdown
    path: Path            # absolute path to the .md file
    body: str             # raw markdown


@dataclass
class HelpSection:
    """A directory in the content tree (e.g. 'de/', 'meta/')."""
    name: str
    title: str            # human-readable
    entries: list[HelpEntry] = field(default_factory=list)


# Section title overrides (otherwise we title-case the directory name).
_SECTION_TITLES = {
    "scrna":   "sc/snRNA-seq",
    "de":      "Differential Expression",
    "meta":    "Meta-Analysis",
    "figures": "Figure Export",
}


def _read_title(md_path: Path) -> str:
    """Extract the first '# Heading' from a markdown file."""
    try:
        with md_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("# "):
                    return line[2:].strip()
    except OSError:
        pass
    return md_path.stem.replace("_", " ").title()


def discover_content(content_dir: Path = CONTENT_DIR
                     ) -> tuple[dict[str, HelpEntry], list[HelpSection]]:
    """Walk 'content_dir' and build the registry + section list.

    Returns
    -------
    entries : dict
        'help_id -> HelpEntry' for every '.md' file found.
    sections : list of HelpSection
        Top-level sections in display order. Files at the top level
        ('index.md', 'glossary.md' etc.) appear in a synthetic
        first section called "Getting started".
    """
    entries: dict[str, HelpEntry] = {}
    sections: list[HelpSection] = []
    if not content_dir.is_dir():
        return entries, sections

    top_level = HelpSection("_top", "Getting started")
    for md_path in sorted(content_dir.glob("*.md")):
        help_id = md_path.stem
        body = md_path.read_text(encoding="utf-8")
        entry = HelpEntry(help_id, _read_title(md_path), md_path, body)
        entries[help_id] = entry
        top_level.entries.append(entry)
    if top_level.entries:
        sections.append(top_level)

    for sub in sorted(p for p in content_dir.iterdir() if p.is_dir()):
        title = _SECTION_TITLES.get(sub.name, sub.name.replace("_", " ").title())
        section = HelpSection(sub.name, title)
        for md_path in sorted(sub.rglob("*.md")):
            rel = md_path.relative_to(content_dir).with_suffix("")
            help_id = rel.as_posix()
            body = md_path.read_text(encoding="utf-8")
            entry = HelpEntry(help_id, _read_title(md_path), md_path, body)
            entries[help_id] = entry
            section.entries.append(entry)
        if section.entries:
            sections.append(section)
    return entries, sections


# Browser dialog

class HelpBrowser(QDialog):
    """Help dialog: TOC tree on the left, markdown viewer on the right.

    Features
    --------
    - Back / Forward / Home navigation with a history stack
    - Search box that filters the tree by title + body substring
    - "View online" link out to the published docs site
    - Singleton-style usage via 'HelpManager'
    """

    DEFAULT_HELP_ID = "index"

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("KOSMIC Help")
        self.resize(900, 640)

        self._entries, self._sections = discover_content()
        self._current_id: Optional[str] = None
        self._history: list[str] = []
        self._history_pos: int = -1

        self._build_ui()
        self._populate_tree()

    # -- UI scaffold --------------------------------------------------------

    def _build_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)

        toolbar = QToolBar()
        self._back_btn = toolbar.addAction("← Back", self.go_back)
        self._fwd_btn = toolbar.addAction("Forward →", self.go_forward)
        toolbar.addSeparator()
        toolbar.addAction("Home", lambda: self.show_id(self.DEFAULT_HELP_ID))
        toolbar.addSeparator()
        self._online_btn = toolbar.addAction(
            "View online", self._open_online)
        # No site yet: the action stays but cannot fire (see ONLINE_BASE).
        self._online_btn.setVisible(ONLINE_BASE is not None)
        outer.addWidget(toolbar)

        self._search = QLineEdit()
        self._search.setPlaceholderText("Search help…")
        self._search.textChanged.connect(self._filter_tree)
        outer.addWidget(self._search)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter, 1)

        self._tree = QTreeWidget()
        self._tree.setHeaderHidden(True)
        self._tree.itemSelectionChanged.connect(self._on_tree_selected)
        splitter.addWidget(self._tree)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        self._browser = QTextBrowser()
        self._browser.setOpenExternalLinks(True)
        self._browser.anchorClicked.connect(self._on_anchor_clicked)
        # We handle internal links manually so they navigate within the dialog.
        self._browser.setOpenLinks(False)
        right_layout.addWidget(self._browser)
        splitter.addWidget(right)

        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([240, 660])

        # Esc closes the dialog.
        QShortcut(QKeySequence(Qt.Key.Key_Escape), self, activated=self.accept)

        self._update_nav_buttons()

    def _populate_tree(self) -> None:
        self._tree.clear()
        for section in self._sections:
            top = QTreeWidgetItem([section.title])
            top.setData(0, Qt.ItemDataRole.UserRole, None)
            top.setExpanded(True)
            self._tree.addTopLevelItem(top)
            for entry in section.entries:
                item = QTreeWidgetItem([entry.title])
                item.setData(0, Qt.ItemDataRole.UserRole, entry.help_id)
                top.addChild(item)
            top.setExpanded(True)
        self._tree.expandAll()

    # -- Navigation ---------------------------------------------------------

    def show_id(self, help_id: str) -> None:
        """Display the entry for 'help_id', recording the visit in history."""
        from kosmic.gui.help.help_render import render_cached
        entry = self._entries.get(help_id)
        if entry is None:
            entry = self._entries.get(self.DEFAULT_HELP_ID)
        if entry is None:
            self._browser.setHtml(render_cached(
                f"# Help unavailable\n\nNo help content found for '{help_id}'."))
            return

        # Truncate any forward history if we navigated away from a back-state.
        if self._history_pos < len(self._history) - 1:
            self._history = self._history[: self._history_pos + 1]
        if not self._history or self._history[-1] != entry.help_id:
            self._history.append(entry.help_id)
            self._history_pos = len(self._history) - 1

        self._render_entry(entry)
        self._select_in_tree(entry.help_id)
        self._update_nav_buttons()

    def _render_entry(self, entry: HelpEntry) -> None:
        from kosmic.gui.help.help_render import render_cached
        self._current_id = entry.help_id
        # search-paths so embedded relative images resolve against the .md
        self._browser.setSearchPaths([str(entry.path.parent)])
        self._browser.setHtml(render_cached(entry.body, self._image_width()))

    def _image_width(self) -> int:
        """Usable width for embedded screenshots, in 40px steps so a
        few pixels of resize do not defeat the render cache."""
        w = self._browser.viewport().width() - 24
        return max(240, (w // 40) * 40)

    def go_back(self) -> None:
        if self._history_pos <= 0:
            return
        self._history_pos -= 1
        help_id = self._history[self._history_pos]
        entry = self._entries.get(help_id)
        if entry is not None:
            self._render_entry(entry)
            self._select_in_tree(help_id)
        self._update_nav_buttons()

    def go_forward(self) -> None:
        if self._history_pos >= len(self._history) - 1:
            return
        self._history_pos += 1
        help_id = self._history[self._history_pos]
        entry = self._entries.get(help_id)
        if entry is not None:
            self._render_entry(entry)
            self._select_in_tree(help_id)
        self._update_nav_buttons()

    def _update_nav_buttons(self) -> None:
        self._back_btn.setEnabled(self._history_pos > 0)
        self._fwd_btn.setEnabled(
            self._history_pos < len(self._history) - 1)
        self._online_btn.setEnabled(
            self._current_id is not None and ONLINE_BASE is not None)

    # -- Slots --------------------------------------------------------------

    def _on_tree_selected(self) -> None:
        items = self._tree.selectedItems()
        if not items:
            return
        help_id = items[0].data(0, Qt.ItemDataRole.UserRole)
        if help_id and help_id != self._current_id:
            self.show_id(help_id)

    def _on_anchor_clicked(self, url: QUrl) -> None:
        """Resolve internal markdown links to other help ids."""
        s = url.toString()
        if url.scheme() in {"http", "https", "mailto"}:
            QDesktopServices.openUrl(url)
            return
        # Treat the link as a relative .md path next to the current entry.
        current = self._entries.get(self._current_id) if self._current_id else None
        if current is None:
            return
        try:
            target = (current.path.parent / s).resolve()
            rel = target.relative_to(CONTENT_DIR.resolve()).with_suffix("")
            self.show_id(rel.as_posix())
        except (OSError, ValueError):
            pass

    def _filter_tree(self, query: str) -> None:
        query = query.strip().lower()
        for i in range(self._tree.topLevelItemCount()):
            top = self._tree.topLevelItem(i)
            any_visible = False
            for j in range(top.childCount()):
                child = top.child(j)
                help_id = child.data(0, Qt.ItemDataRole.UserRole)
                entry = self._entries.get(help_id) if help_id else None
                visible = (
                    not query
                    or (entry is not None
                        and (query in entry.title.lower()
                             or query in entry.body.lower()))
                )
                child.setHidden(not visible)
                any_visible = any_visible or visible
            top.setHidden(not any_visible)

    def _select_in_tree(self, help_id: str) -> None:
        self._tree.blockSignals(True)
        try:
            for i in range(self._tree.topLevelItemCount()):
                top = self._tree.topLevelItem(i)
                for j in range(top.childCount()):
                    child = top.child(j)
                    if child.data(0, Qt.ItemDataRole.UserRole) == help_id:
                        self._tree.setCurrentItem(child)
                        return
        finally:
            self._tree.blockSignals(False)

    def _open_online(self) -> None:
        if self._current_id is not None:
            HelpManager.instance().open_online(self._current_id)


# Manager

class HelpManager:
    """Process-wide singleton: routes help requests to the site or the
    offline browser, and owns the one HelpBrowser dialog."""

    _instance: Optional["HelpManager"] = None

    def __init__(self) -> None:
        self._browser: Optional[HelpBrowser] = None
        self._main_window = None  # type: ignore[assignment]

    @classmethod
    def instance(cls) -> "HelpManager":
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def attach_main_window(self, main_window) -> None:
        """Remember the window the offline browser should be parented to."""
        self._main_window = main_window

    # -- entry points ----------------------------------------------------

    def open(self, help_id: Optional[str] = None) -> None:
        """F1: the page for 'help_id' on the website, else offline."""
        help_id = help_id or HelpBrowser.DEFAULT_HELP_ID
        if not self.open_online(help_id):
            self.show_full(help_id)

    @staticmethod
    def online_url(help_id: str, base: Optional[str] = None) -> Optional[str]:
        base = base or ONLINE_BASE
        if not base:
            return None
        # 'scrna/index' is the section page, served at '<site>/scrna/'.
        if help_id == HelpBrowser.DEFAULT_HELP_ID:
            page = ""
        elif help_id.endswith("/index"):
            page = help_id[:-len("index")]
        else:
            page = f"{help_id}/"
        return f"{base}{page}"

    def open_online(self, help_id: str) -> bool:
        """Open the page in the system browser: on the documentation
        site if it answers, else on a local 'mkdocs serve' when this is
        a source checkout that can run one. False when neither applies,
        so the caller shows the bundled copy instead of a browser error
        page."""
        url = self.online_url(help_id)
        if url is not None and _host_reachable(url):
            return bool(QDesktopServices.openUrl(QUrl(url)))
        local = self.online_url(help_id, base=LOCAL_DOCS_URL)
        if local is not None and self._start_local_docs(local):
            return bool(QDesktopServices.openUrl(QUrl(local)))
        return False

    # -- developer convenience: serve the docs locally when the site is
    #    not reachable and this is a source checkout ---------------------

    _docs_proc = None

    def _start_local_docs(self, url: str) -> bool:
        """Start 'mkdocs serve' at the loopback 'url' if it is not already
        answering, and wait for it. Only possible in a source checkout
        with mkdocs installed -- a packaged install has neither, so this
        returns False there and the bundled copy is shown. The server
        dies with the app."""
        if _host_reachable(url):
            return True
        import atexit
        import importlib.util
        import subprocess
        import sys
        import time
        from urllib.parse import urlparse

        parsed = urlparse(url)
        if parsed.hostname not in ("127.0.0.1", "localhost", "::1"):
            return False
        repo = Path(__file__).resolve().parents[3]
        if not (repo / "mkdocs.yml").is_file():
            return False
        if importlib.util.find_spec("mkdocs") is None:
            return False
        if self._docs_proc is not None and self._docs_proc.poll() is None:
            return _host_reachable(url)   # started earlier, still coming up?

        addr = f"{parsed.hostname}:{parsed.port or 8000}"
        creation = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        # '--livereload' is the documented default, but with click >= 8.2
        # mkdocs 1.6 ends up with it False unless passed explicitly, and
        # then never rebuilds on edit.
        self._docs_proc = subprocess.Popen(
            [sys.executable, "-m", "mkdocs", "serve", "--livereload",
             "-a", addr],
            cwd=str(repo), stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, creationflags=creation)
        atexit.register(self._stop_local_docs)
        # A cold build of the site takes a second or two.
        for _ in range(40):
            if _host_reachable(url, timeout=0.25):
                return True
            if self._docs_proc.poll() is not None:
                return False
            time.sleep(0.25)
        return False

    def _stop_local_docs(self) -> None:
        proc = self._docs_proc
        if proc is not None and proc.poll() is None:
            proc.terminate()

    def show_full(self, help_id: Optional[str] = None,
                  parent: Optional[QWidget] = None) -> None:
        """The offline browser dialog (TOC + search + history)."""
        if self._browser is None:
            self._browser = HelpBrowser(parent or self._main_window)
        if not self._browser.isVisible():
            self._browser.show()
        self._browser.raise_()
        self._browser.activateWindow()
        self._browser.show_id(help_id or HelpBrowser.DEFAULT_HELP_ID)

    def refresh_theme(self) -> None:
        if self._browser is not None:
            self._browser.refresh_theme()


def _host_reachable(url: str, timeout: float = 1.0) -> bool:
    """A TCP connect to the site's host. Cheap enough for a key press and
    avoids handing the user a browser error page when offline."""
    import socket
    from urllib.parse import urlparse
    parsed = urlparse(url)
    host = parsed.hostname
    if not host:
        return False
    port = parsed.port or (443 if parsed.scheme == "https" else 80)
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except OSError:
        return False
