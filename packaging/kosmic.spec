# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller build for KOSMIC.

    pip install -e ".[build]"
    pyinstaller packaging/kosmic.spec

Produces dist/KOSMIC/ (one-directory build: KOSMIC.exe beside an
_internal/ folder). One-directory rather than one-file because the
scientific stack is large and a one-file build would unpack it to a
temp folder on every launch.

Data files are bundled under their package-relative paths, because the
code locates them from ``__file__`` (``importlib.resources`` for the
reference data, ``Path(__file__).parent`` for icons and help content)
and ``config.toml`` from the package's parent. Keeping the tree shape
means nothing in the code needs to know it is frozen.
"""
import os
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_submodules

# KOSMIC_BUILD_CONSOLE=1 keeps a console window so a startup traceback
# is visible; the release build has none.
CONSOLE = os.environ.get("KOSMIC_BUILD_CONSOLE") == "1"

ROOT = Path(SPECPATH).resolve().parent      # repo root; SPECPATH is packaging/
PKG = ROOT / "kosmic"


def package_data(subdir: str, patterns=("*",)) -> list[tuple[str, str]]:
    """Every non-Python file under kosmic/<subdir>, at the same relative path."""
    out = []
    base = PKG / subdir
    for pat in patterns:
        for p in base.rglob(pat):
            if p.is_file() and p.suffix != ".py" and "__pycache__" not in p.parts:
                out.append((str(p), str(Path("kosmic") / p.relative_to(PKG).parent)))
    return out


datas = [
    (str(ROOT / "config.toml"), "."),
    *package_data("reference"),
    *package_data("gui/shared/icons"),
    *package_data("gui/shared/assets"),
    *package_data("gui/help"),           # content/*.md, images, about.md, tutorial/steps.toml
]

# Packages whose submodules or data PyInstaller's static analysis misses.
# collect_all pulls modules, data and binaries for each.
hiddenimports = []
binaries = []
for pkg in (
    "scanpy", "anndata", "pydeseq2", "gseapy", "goatools", "leidenalg",
    "igraph", "harmonypy", "umap", "pynndescent", "numba", "llvmlite",
    "statsmodels", "formulaic", "formulaic_contrasts", "pyqtgraph",
    "seaborn", "markdown_it", "mdurl", "h5py", "zarr", "numcodecs",
    "skmisc", "sklearn", "scipy", "pydot", "natsort", "session_info",
    "legacy_api_wrap", "array_api_compat", "psutil", "openpyxl",
    "xlsxwriter", "imageio", "tifffile", "skimage", "networkx",
):
    try:
        d, b, h = collect_all(pkg)
    except Exception:
        continue
    datas += d
    binaries += b
    hiddenimports += h

hiddenimports += collect_submodules("kosmic")

# collect_all is generous: drop what no user needs. goatools.test_data is
# ~100 MB of test fixtures shipped as Python modules.
def _keep(entry):
    return "goatools/test_data" not in entry[0].replace("\\", "/")


datas = [d for d in datas if _keep(d)]
hiddenimports = [h for h in hiddenimports if not h.startswith("goatools.test_data")]

excludes = [
    # Never used by the app; some are large.
    "torch", "pyarrow", "goatools.test_data", "tkinter", "PyQt5", "PySide2",
    "PySide6", "IPython",
    "jupyter", "notebook", "pytest", "mkdocs", "mkdocs_material",
    "mkdocstrings", "ruff", "PyInstaller",
]

a = Analysis(
    [str(ROOT / "main.py")],
    pathex=[str(ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    exclude_binaries=True,
    name="KOSMIC",
    debug=False,
    strip=False,
    upx=False,
    console=CONSOLE,        # GUI app: no console window in a release build
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="KOSMIC",
)
