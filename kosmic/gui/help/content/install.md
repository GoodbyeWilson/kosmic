# Installation

KOSMIC is a desktop application for Windows, macOS and Linux. It needs
a graphical session; it does not run on a headless server. Allow 8 GB
of memory as a minimum and 16–32 GB for datasets of several hundred
thousand cells. Windows 11 is the most heavily tested platform, Ubuntu
24.04 moderately, macOS least — please report problems there.

There are two ways to install it.

## Installer

One-click installers for Windows, macOS and Linux will be published on
the project's GitHub Releases page. They include Python and every
dependency; nothing else needs to be installed. Until the first release
is published, install from source as below.

## From source

You need two things installed first. Both are free.

**Python 3.11 or newer.** Check what you have by opening a terminal
(Windows: press the Windows key, type `cmd`, Enter; macOS: Terminal;
Linux: your usual shell) and running:

```sh
python --version
```

If that says 3.11 or higher, you are set. If it says an older version,
or "not recognized" / "not found", install it from
[python.org/downloads](https://www.python.org/downloads/).

> **Note:** On Windows, tick **Add python.exe to PATH** on the first
> screen of the installer. It is unticked by default, and leaving it
> unticked is the most common reason the next step says `python` is not
> recognised. If you already installed without it, re-run the installer,
> choose *Modify*, and tick it there.

On macOS the system `python3` is often too old; the python.org
installer is the reliable option. On Linux, `python3.11` from your
package manager is fine (`sudo apt install python3.11 python3.11-venv`
on Ubuntu/Debian).

**Git**, to download the code, from
[git-scm.com](https://git-scm.com/downloads); the installer's defaults
are fine. If you would rather not, GitHub's *Code → Download ZIP*
button works too — unzip it and skip the `git clone` line below.

### Install

Open a terminal in the folder where you want KOSMIC to live, then:

```sh
git clone https://github.com/GoodbyeWilson/kosmic.git
cd kosmic
python -m venv .venv
```

Then activate the environment. This is the one step that differs by
operating system, and you repeat it every time you open a new terminal
to use KOSMIC:

```sh
.venv\Scripts\activate          # Windows (cmd or PowerShell)
source .venv/bin/activate       # macOS / Linux
```

Your prompt should now start with `(.venv)`. Finally:

```sh
pip install -e ".[dev]"
python main.py
```

The install takes about five minutes and downloads roughly a gigabyte
(the scientific Python stack). No compiler is needed on any platform.

> **Note:** If PowerShell refuses to run the activate script ("running
> scripts is disabled on this system"), run
> `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once, answer Y,
> and try again. This is Windows' default policy, not a KOSMIC problem.

If you already use conda, skip the venv lines: `conda create -n kosmic
python=3.11`, `conda activate kosmic`, then the same
`pip install -e ".[dev]"`. Everything comes from PyPI either way.

## Optional components

**R, only if you will import Seurat `.rds` or `.Robj` files.** KOSMIC
converts them by running R as a separate process, so R must be
installed with the **Seurat** and **Matrix** packages (**SeuratDisk** is
used if present). Get R from [cran.r-project.org](https://cran.r-project.org/),
then in R:

```r
install.packages(c("Seurat", "Matrix"))
```

KOSMIC finds R automatically — on `PATH`, or under `C:\Program Files\R`
on Windows (newest version wins). If you only load `.h5ad`, `.h5`,
`.mtx` or 10x folders, you do not need R.

**Optional Python features**, installed as extras (quote the brackets —
`zsh` treats bare brackets as a pattern):

```sh
pip install -e ".[all]"        # everything below
pip install -e ".[annotate]"   # CellTypist cell-type annotation
pip install -e ".[rds]"        # .rds reading without R (simple objects only)
pip install -e ".[tsne]"       # openTSNE backend
pip install -e ".[monitor]"    # status-bar CPU / RAM monitor
```

Two features depend on packages that are not on PyPI and are installed
separately if you want them: **SoupX**-style ambient-RNA estimation
(the QC step's Ambient stage) and **DecontX** decontamination (the
Decontaminate step). Each is imported only when used, and KOSMIC says
what is missing rather than failing at startup.

## Checking it works

Launch with `python main.py`. The window opens on the Home page; the
Project workspace is where an analysis starts — see
[The Project workspace](project/index.md).

To run the test suite (a few minutes; useful if you are contributing):

```sh
pytest tests/
```
