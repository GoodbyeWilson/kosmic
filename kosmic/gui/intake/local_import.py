# Local-file intake (ADR-001 Phase E): file picker + copy into the
# target study, without leaving the Project page. Every format lands in
# raw_data/ and is imported from the study's Dataset screen -- h5ad
# included, because which slot holds raw counts varies by depositor
# (CellxGene normalises X and keeps counts in .raw) and copying one
# straight into processed_data/ would install log-normalised values as
# the working dataset without anyone looking.
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import QFileDialog

from kosmic.paths import raw_data_dir

_FILTERS = {
    'h5ad': "AnnData Files (*.h5ad);;All Files (*)",
    'csv': "Tabular Files (*.csv *.csv.gz *.tsv *.tsv.gz);;All Files (*)",
    'rds': ("Seurat Files (*.rds *.RDS *.robj *.Robj *.robj.gz "
            "*.Robj.gz);;All Files (*)"),
}


def import_local_file(parent, study_dir, kind: str, log_cb) -> Optional[Path]:
    """Pick a local file and copy it into 'study_dir'.

    Returns the copied path, or None if the picker was cancelled.
    """
    study = Path(study_dir)
    path, _ = QFileDialog.getOpenFileName(
        parent, f"Import {kind.upper()} File", "",
        _FILTERS.get(kind, "All Files (*)"))
    if not path:
        return None

    src = Path(path)
    raw = raw_data_dir(study)
    raw.mkdir(exist_ok=True)
    dest = raw / src.name
    if not dest.exists():
        shutil.copy2(src, dest)

    log_cb(f"Imported {src.name} into {study.name}/raw_data")
    return dest
