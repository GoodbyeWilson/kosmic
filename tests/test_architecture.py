"""Architecture invariants (see CONTRIBUTING.md, "Architecture").

These tests encode the structural rules from CLAUDE.md and the ADRs so
that every subsequent cleanup phase has a tripwire if it accidentally
drifts:

1. No module in `src/` imports Qt (PyQt6 / PyQt5 / PySide6 / PySide2).
   src/ is the pure-logic layer; any Qt import there is a layering break.

2. No module in `src/` imports from `gui/`.
   Reverse-direction coupling would make `src/` depend on the GUI layer.

3. Every top-level module under `src/`, `kosmic/gui/tabs/`, `kosmic/gui/de_analysis/`,
   `kosmic/gui/meta_analysis/` imports cleanly under the `kirk` env.

4. `kosmic/gui/shared/` may not import from any workspace
   subpackage. The shared toolbox knows nothing about its consumers.

5. Every concrete page/tab class in a workspace subpackage
   must inherit from one of the four page bases (SimplePage,
   TabbedPage, SidebarPage, SidebarTabbedPage) or FigurePage. No page
   may hand-roll its own QSplitter + sidebar geometry.

Parsed via ``ast`` -- docstring mentions and comments do not trigger.
"""
from __future__ import annotations

import ast
import importlib
import os
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_KOSMIC = os.path.join(_ROOT, "kosmic")
_KOSMIC_GUI = os.path.join(_KOSMIC, "gui")

_QT_PREFIXES = ("PyQt6", "PyQt5", "PySide6", "PySide2")
_GUI_PREFIX = "kosmic.gui"


def _iter_python_files(root: str, *, skip_dir: str | None = None):
    """Walk *root* yielding .py files; skip the *skip_dir* subtree if given."""
    for dirpath, dirs, files in os.walk(root):
        if skip_dir is not None and os.path.commonpath(
                [os.path.abspath(dirpath), os.path.abspath(skip_dir)]
                ) == os.path.abspath(skip_dir):
            continue
        for f in files:
            if f.endswith(".py") and not f.startswith("._"):
                yield os.path.join(dirpath, f)


def _iter_pure_kosmic_files():
    """All .py files under kosmic/ EXCEPT those in kosmic/gui/."""
    yield from _iter_python_files(_KOSMIC, skip_dir=_KOSMIC_GUI)


def _imported_modules(path: str):
    """Yield every module name imported by `path`, via AST."""
    with open(path, "r", encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.module is not None:
                yield node.module


# ---------------------------------------------------------------------
# Rule 1: no Qt imports in src/
# ---------------------------------------------------------------------

@pytest.mark.parametrize(
    "filepath",
    list(_iter_pure_kosmic_files()),
    ids=lambda p: os.path.relpath(p, _ROOT),
)
def test_no_qt_imports_in_src(filepath: str):
    """src/ must stay free of Qt imports (CLAUDE.md principle 2)."""
    offenders = [
        m for m in _imported_modules(filepath)
        if any(m == q or m.startswith(q + ".") for q in _QT_PREFIXES)
    ]
    assert not offenders, (
        f"{os.path.relpath(filepath, _ROOT)} imports Qt: {offenders}. "
        "src/ is the pure-logic layer (CLAUDE.md principle 2)."
    )


# ---------------------------------------------------------------------
# Rule 2: no src -> gui imports
# ---------------------------------------------------------------------

@pytest.mark.parametrize(
    "filepath",
    list(_iter_pure_kosmic_files()),
    ids=lambda p: os.path.relpath(p, _ROOT),
)
def test_src_does_not_import_gui(filepath: str):
    """Reverse-direction coupling (src -> gui) is forbidden."""
    offenders = [
        m for m in _imported_modules(filepath)
        if m == _GUI_PREFIX or m.startswith(_GUI_PREFIX + ".")
    ]
    assert not offenders, (
        f"{os.path.relpath(filepath, _ROOT)} imports from gui/: {offenders}. "
        "gui/ depends on src/, never the reverse."
    )


# ---------------------------------------------------------------------
# Rule 3: every module imports cleanly
# ---------------------------------------------------------------------

def _as_module(path: str) -> str:
    """Convert a file path under the repo root into a dotted module name."""
    rel = os.path.relpath(path, _ROOT)
    if rel.endswith("__init__.py"):
        rel = os.path.dirname(rel)
    else:
        rel = rel[:-3]  # strip .py
    return rel.replace(os.sep, ".")


def _collect_import_targets():
    targets = []
    # kosmic/ (excluding kosmic/gui/) -- pure-logic; import should be cheap.
    for p in _iter_pure_kosmic_files():
        if os.path.basename(p) == "__init__.py" and p.endswith(
                os.sep + "kosmic" + os.sep + "__init__.py"):
            continue
        targets.append(_as_module(p))

    # kosmic/gui/* workspace packages -- importing them requires PyQt6 but no QApplication.
    for sub in ("kosmic/gui/scrna", "kosmic/gui/de_analysis",
                "kosmic/gui/meta_analysis", "kosmic/gui/figure_export",
                "kosmic/gui/shared", "kosmic/gui/help"):
        sub_abs = os.path.join(_ROOT, sub)
        if not os.path.isdir(sub_abs):
            continue
        for p in _iter_python_files(sub_abs):
            # Skip empty __init__.py files (they have nothing to test).
            if os.path.basename(p) == "__init__.py":
                with open(p, "r", encoding="utf-8") as fh:
                    if not fh.read().strip():
                        continue
            targets.append(_as_module(p))
    return sorted(set(targets))


# Ensure Qt is available headlessly for gui/* imports.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.mark.parametrize("module_name", _collect_import_targets())
def test_module_imports_cleanly(module_name: str):
    """Every active src/ + gui/{tabs,de_analysis,meta_analysis} module must
    import without raising."""
    try:
        importlib.import_module(module_name)
    except Exception as exc:  # noqa: BLE001 -- any failure is the signal here
        pytest.fail(f"import {module_name!r} raised {type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------
# Rule 4: kosmic/gui/shared/ does not import from workspaces
# ---------------------------------------------------------------------

_SHARED_DIR = os.path.join(_KOSMIC_GUI, "shared")
_WORKSPACE_PREFIXES = (
    "kosmic.gui.scrna",
    "kosmic.gui.de_analysis",
    "kosmic.gui.meta_analysis",
    "kosmic.gui.figure_export",
    "kosmic.gui.help",
)


def _iter_shared_files():
    yield from _iter_python_files(_SHARED_DIR)


@pytest.mark.parametrize(
    "filepath",
    list(_iter_shared_files()),
    ids=lambda p: os.path.relpath(p, _ROOT),
)
def test_shared_does_not_import_workspaces(filepath: str):
    """`kosmic/gui/shared/` is the toolbox -- it must not depend on any
    workspace subpackage, otherwise toolbox additions can pull in
    workspace-specific imports and break the dependency direction
    declared in CONTRIBUTING.md."""
    offenders = [
        m for m in _imported_modules(filepath)
        if any(m == prefix or m.startswith(prefix + ".")
               for prefix in _WORKSPACE_PREFIXES)
    ]
    assert not offenders, (
        f"{os.path.relpath(filepath, _ROOT)} imports from a workspace: "
        f"{offenders}. The shared toolbox may not depend on its consumers "
        "(see CONTRIBUTING.md)."
    )


# ---------------------------------------------------------------------
# Rule 5: pages inherit from a base
# ---------------------------------------------------------------------

# Files that contain a "page" or "tab" class but legitimately don't fit
# the base-class taxonomy. Keep this list short -- every entry is a
# justification you'll need to give a future reviewer.
_PAGE_BASE_EXCEPTIONS = {
    # AnnotateTab is structurally embedded inside ClusterTab. Documented
    # in cluster_tab.py:_setup_ui and annotate_tab.py docstring.
    "kosmic/gui/scrna/tabs/annotate_tab.py",
    # MethodsComparisonPage hosts a SidebarPage via composition because
    # it has a page-title row above the sidebar+content body that the
    # SidebarPage layout doesn't accommodate. The composition is
    # documented in methods_comparison_page._setup_ui.
    "kosmic/gui/meta_analysis/pages/methods_comparison_page.py",
}

# Directories whose contents are page-helpers / sub-widgets, not the
# top-level page class. The split that extracted gene_ma_page's
# tabs into ``pages/gene_ma/`` -- those tab widgets are children of
# ``GeneMAPage`` (which itself still inherits a page base) and don't
# need to inherit one themselves.
_PAGE_BASE_EXCEPTION_DIRS = {
    "kosmic/gui/meta_analysis/pages/gene_ma",
}

_PAGE_DIRS = (
    "kosmic/gui/scrna/tabs",
    "kosmic/gui/de_analysis/pages",
    "kosmic/gui/meta_analysis/pages",
    "kosmic/gui/figure_export/pages",
)


def _iter_page_files():
    """Files that should contain page/tab classes inheriting from a base."""
    for sub in _PAGE_DIRS:
        sub_abs = os.path.join(_ROOT, sub)
        if not os.path.isdir(sub_abs):
            continue
        for p in _iter_python_files(sub_abs):
            if os.path.basename(p) in ("__init__.py", "_base.py"):
                continue
            rel = os.path.relpath(p, _ROOT).replace(os.sep, "/")
            if rel in _PAGE_BASE_EXCEPTIONS:
                continue
            if any(rel.startswith(d + "/") for d in _PAGE_BASE_EXCEPTION_DIRS):
                continue
            yield p


_PAGE_BASES = (
    "SimplePage",
    "TabbedPage",
    "SidebarPage",
    "SidebarTabbedPage",
    "FigurePage",
)


def _class_bases(node: ast.ClassDef) -> list[str]:
    """Return the simple names of a class's base classes."""
    out = []
    for b in node.bases:
        if isinstance(b, ast.Name):
            out.append(b.id)
        elif isinstance(b, ast.Attribute):
            out.append(b.attr)
    return out


def _page_class_bases() -> dict:
    """Map every page-directory class name to its declared bases.

    Lets the check below follow inheritance instead of demanding an
    archetype be named directly: a page that subclasses another KOSMIC
    page inherits that page's archetype, and so is not hand-rolling
    geometry -- which is the thing this test exists to catch.
    """
    out: dict[str, list[str]] = {}
    for path in _iter_page_files():
        with open(path, "r", encoding="utf-8") as fh:
            tree = ast.parse(fh.read(), filename=path)
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                out.setdefault(node.name, _class_bases(node))
    return out


_PAGE_CLASS_BASES = _page_class_bases()


def _reaches_page_base(bases, seen=None) -> bool:
    """True if any base is an archetype, or transitively reaches one."""
    seen = seen or set()
    for b in bases:
        if b in _PAGE_BASES:
            return True
        if b in seen or b not in _PAGE_CLASS_BASES:
            continue
        seen.add(b)
        if _reaches_page_base(_PAGE_CLASS_BASES[b], seen):
            return True
    return False


@pytest.mark.parametrize(
    "filepath",
    list(_iter_page_files()),
    ids=lambda p: os.path.relpath(p, _ROOT),
)
def test_page_inherits_from_base(filepath: str):
    """Every concrete *Tab / *Page class in a workspace's tabs/ or
    pages/ directory inherits from one of the four page bases (or
    ``FigurePage``).

    This catches regressions where someone adds a new page that
    hand-rolls its own QSplitter + sidebar geometry instead of
    subclassing ``SidebarPage`` (or whichever archetype fits).
    """
    with open(filepath, "r", encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=filepath)

    page_classes = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.ClassDef)
        and (node.name.endswith("Page") or node.name.endswith("Tab"))
        and not node.name.startswith("_")
    ]

    if not page_classes:
        return  # Helper / dialog / worker file -- not a page.

    rel = os.path.relpath(filepath, _ROOT)
    for cls in page_classes:
        bases = _class_bases(cls)
        if not _reaches_page_base(bases):
            pytest.fail(
                f"{rel}::{cls.name} does not inherit from a page base. "
                f"Found bases: {bases}. Expected one of {_PAGE_BASES}."
            )
