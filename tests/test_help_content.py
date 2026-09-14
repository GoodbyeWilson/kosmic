"""The F1 help content is also the documentation site (mkdocs.yml), so
every page the app can ask for must exist, every page must be in the
site's nav, and links between pages must resolve."""
from __future__ import annotations

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONTENT = ROOT / "kosmic" / "gui" / "help" / "content"
GUI = ROOT / "kosmic" / "gui"

_HELP_ID_RE = re.compile(r'help_id\s*=\s*"([a-z0-9_/]+)"')


def _help_ids_in_code() -> set[str]:
    ids = set()
    for py in GUI.rglob("*.py"):
        ids.update(_HELP_ID_RE.findall(py.read_text(encoding="utf-8")))
    # The Project workspace serves one id per step through a property.
    from kosmic.gui.project.workspace import ProjectWorkspace
    ids.update(ProjectWorkspace._STEP_HELP)
    return ids


def _pages() -> set[str]:
    return {p.relative_to(CONTENT).with_suffix("").as_posix()
            for p in CONTENT.rglob("*.md")}


# Pages the app asks for that have not been written yet. Shrink this as
# they land; the test fails if one is written and left listed here, so
# it cannot quietly go stale.
KNOWN_MISSING = {
    "meta/enrichment", "meta/methods", "meta/validation",
}


def test_every_help_id_has_a_page():
    missing = _help_ids_in_code() - _pages()
    new = sorted(missing - KNOWN_MISSING)
    assert not new, f"help_id without a page: {new}"
    written = sorted(KNOWN_MISSING - missing)
    assert not written, f"now written -- remove from KNOWN_MISSING: {written}"


def _nav_pages(node) -> set[str]:
    out = set()
    if isinstance(node, dict):
        for v in node.values():
            out |= _nav_pages(v)
    elif isinstance(node, list):
        for v in node:
            out |= _nav_pages(v)
    elif isinstance(node, str):
        if node.endswith("/"):
            return out      # a generated section (literate-nav), not a file
        out.add(node.removesuffix(".md"))
    return out


def test_every_page_is_in_the_site_nav():
    cfg = yaml.safe_load((ROOT / "mkdocs.yml").read_text(encoding="utf-8"))
    in_nav = _nav_pages(cfg["nav"])
    assert cfg["docs_dir"] == "kosmic/gui/help/content"
    missing = sorted(_pages() - in_nav)
    assert not missing, f"pages not in mkdocs nav: {missing}"
    dangling = sorted(in_nav - _pages())
    assert not dangling, f"nav entries without a page: {dangling}"


_LINK_RE = re.compile(r"\]\(([^)#\s]+\.md)(#[^)]*)?\)")


def test_internal_links_resolve():
    broken = []
    for page in CONTENT.rglob("*.md"):
        for target, _anchor in _LINK_RE.findall(page.read_text(encoding="utf-8")):
            if target.startswith(("http://", "https://")):
                continue
            if not (page.parent / target).resolve().is_file():
                broken.append(f"{page.relative_to(CONTENT)} -> {target}")
    assert not broken, "\n".join(broken)


def test_project_help_follows_the_step():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])  # noqa: F841 -- keep alive
    from kosmic.gui.project.workspace import ProjectWorkspace
    from kosmic.project_workflow import STEP_REVIEW
    ws = ProjectWorkspace()
    assert ws.help_id == "project/project"
    ws._project_dir = Path(".")            # a step needs an open project
    ws.on_sidebar_step(STEP_REVIEW)
    assert ws.help_id == "project/review"
    ws.deleteLater()
