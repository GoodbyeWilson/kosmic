"""Generate the developer reference for the docs site from the source tree.

Runs inside `mkdocs build` (the mkdocs-gen-files plugin, see mkdocs.yml)
and writes one virtual page per module under `reference/`, plus the
navigation for that section. Nothing is written to disk in the repo; the
pages exist only in the built site, so they cannot go stale.

Each page carries:

- the module's leading comment block, if it has one instead of a
  docstring (most analysis modules explain themselves that way and
  mkdocstrings only renders docstrings);
- the module's API rendered by mkdocstrings from its docstrings.

The analysis package (`kosmic/` minus `kosmic/gui/`) is documented in
full. GUI modules get their header only: their public surface is Qt
slots and widgets, which is not what a reader of the reference wants.
"""
from __future__ import annotations

from pathlib import Path

import mkdocs_gen_files

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "kosmic"

nav = mkdocs_gen_files.Nav()

# Navigation order follows the layering in CODEBASE.md, not the
# alphabet: the core modules, then the pipeline in the order data
# moves through it, then the GUI. Flat modules at the package root are
# grouped under 'core'.
ORDER = [
    "core", "scrna", "de", "meta_analysis", "combine",
    "visualisation", "reference", "gui",
]
CORE_ORDER = ["paths", "manifest", "provenance", "project_workflow", "numerical"]
# Sub-packages in pipeline order where the alphabet would mislead.
SUB_ORDER = {
    "scrna": ["load", "inspect", "qc", "cluster", "annotate"],
    "gui": ["shared", "home", "project", "intake", "combine", "scrna",
            "de_analysis", "meta_analysis", "figure_export", "help"],
}


def _sub_rank(group, sub):
    """Sort key for a path below a group: pipeline order for the first
    component where one is declared, alphabetical otherwise."""
    order = SUB_ORDER.get(group, [])
    first = sub[0] if sub else ""
    return (order.index(first) if first in order else len(order), *sub)
TITLES = {
    "core": "Core (paths, manifest, provenance)",
    "scrna": "scRNA: single-study processing",
    "de": "Differential expression",
    "meta_analysis": "Meta-analysis",
    "combine": "Shared atlas (combine)",
    "visualisation": "Visualisation",
    "reference": "Reference data",
    "gui": "GUI",
}

entries = []   # (sort key, nav key, doc path)


def leading_comment(src: str) -> str:
    """The `#` block at the top of a module, as Markdown paragraphs."""
    lines = []
    for raw in src.splitlines():
        if raw.startswith("#!") or raw.startswith("# -*-"):
            continue
        if raw.startswith("#"):
            lines.append(raw[1:].lstrip(" ") if raw != "#" else "")
            continue
        if raw.strip() == "" and not lines:
            continue
        break
    text = "\n".join(lines).strip()
    return text


def module_has_docstring(src: str) -> bool:
    import ast
    try:
        return bool(ast.get_docstring(ast.parse(src)))
    except SyntaxError:
        return False


for path in sorted(PKG.rglob("*.py")):
    rel = path.relative_to(ROOT)
    if "__pycache__" in rel.parts:
        continue
    parts = list(rel.with_suffix("").parts)          # ['kosmic', 'scrna', 'qc', 'decontx']
    if parts[-1] == "__init__":
        parts = parts[:-1]
        if len(parts) == 1:
            continue                                  # the package root has nothing to say
        doc_path = Path("reference", *parts, "index.md")
    else:
        doc_path = Path("reference", *parts).with_suffix(".md")
    module = ".".join(parts)
    is_gui = "gui" in parts

    src = path.read_text(encoding="utf-8", errors="ignore")
    header = "" if module_has_docstring(src) else leading_comment(src)

    with mkdocs_gen_files.open(doc_path, "w") as f:
        f.write(f"# `{module}`\n\n")
        f.write(f"Source: `{rel.as_posix()}`\n\n")
        if header:
            f.write(header + "\n\n")
        if is_gui:
            # Header only; the API of a Qt page is not reference material.
            if module_has_docstring(src):
                f.write(f"::: {module}\n    options:\n      members: false\n")
        else:
            f.write(f"::: {module}\n")

    sub = parts[1:]                                   # drop the 'kosmic' root
    is_package = doc_path.name == "index.md"
    if is_package:
        # kosmic/de/__init__.py -> the section's own index page
        group, key = sub[0], (TITLES[sub[0]],) + tuple(sub[1:])
        rank = (ORDER.index(group), _sub_rank(group, sub[1:]))
    elif len(sub) == 1:
        # kosmic/paths.py -> the virtual 'core' group, in a fixed order
        group, key = "core", (TITLES["core"], sub[0])
        rank = (ORDER.index(group), (CORE_ORDER.index(sub[0])
                                     if sub[0] in CORE_ORDER else 99, sub[0]))
    else:
        group, key = sub[0], (TITLES[sub[0]],) + tuple(sub[1:])
        rank = (ORDER.index(group), _sub_rank(group, sub[1:]))
    entries.append((rank, key, doc_path))
    mkdocs_gen_files.set_edit_path(doc_path, rel.as_posix())

for _, key, doc_path in sorted(entries, key=lambda e: e[0]):
    nav[key] = doc_path.relative_to("reference").as_posix()

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as f:
    f.writelines(nav.build_literate_nav())

with mkdocs_gen_files.open("reference/index.md", "w") as f:
    f.write(
        "# Developer reference\n\n"
        "Generated from the source tree on every build: one page per module, "
        "carrying the module's own explanation and the API rendered from its "
        "docstrings. It cannot drift from the code because it is not written "
        "by hand.\n\n"
        "`kosmic/` (everything except `kosmic/gui/`) is the analysis package "
        "and is documented in full; it imports no Qt and is what a reader "
        "checking a method wants. `kosmic/gui/` pages show the module's "
        "explanation only.\n\n"
        "For how the pieces fit together — the two-layer split, the page "
        "archetypes, the data hand-offs between workspaces — read "
        "`CONTRIBUTING.md` and `CODEBASE.md` in the repository.\n"
    )
