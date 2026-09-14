# Study manifest -- a small per-study JSON cache ('<study>/study.json')
# of dataset identity and semantics (sample / condition / cell-type
# columns, role map), so the Project workspace, DE, and Meta can show
# and pre-fill them without loading the h5ad.
#
# Rules (ADR-001): the manifest is a cache, never a second source of
# truth -- the h5ad and the pseudobulk 'role' column stay authoritative.
# The Inspect tab's save flow is the only writer. Absence is always
# tolerated ("not configured yet"). Asset flags are never stored here;
# they are derived from disk via 'paths.study_status'.
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Optional, Union

from kosmic.paths import study_status

PathLike = Union[str, Path]

MANIFEST_NAME = "study.json"
MANIFEST_VERSION = 1


def manifest_path(study_dir: PathLike) -> Path:
    """Return the manifest path for a study folder."""
    return Path(study_dir) / MANIFEST_NAME


def read_manifest(study_dir: PathLike) -> Optional[dict]:
    """Return the study's manifest dict, or None when absent/unreadable.

    Tolerant by design: a missing, corrupt, or wrong-shaped file means
    "not configured yet" and never raises.
    """
    path = manifest_path(study_dir)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def write_manifest(study_dir: PathLike, *, dataset: Optional[dict] = None,
                   semantics: Optional[dict] = None) -> Path:
    """Write the study manifest, replacing any previous one.

    'dataset' identifies the working file ('file', 'n_cells',
    'n_genes'); 'semantics' carries 'sample_column', 'condition_column',
    'cell_type_column', and 'role_map' (condition value -> role).
    None-valued semantic keys are dropped so consumers can treat
    missing keys uniformly.
    """
    record = {
        "version": MANIFEST_VERSION,
        "updated": datetime.now().isoformat(timespec="seconds"),
    }
    if dataset:
        record["dataset"] = dict(dataset)
    if semantics:
        record["semantics"] = {
            k: v for k, v in semantics.items() if v not in (None, "", {})
        }
    path = manifest_path(study_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return path


def _decode(value) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def _peek_h5ad(study_dir: Path) -> dict:
    """Derive dataset identity and semantics from the newest processed
    h5ad via cheap h5py metadata reads (no matrix data loaded). Used
    when the manifest hasn't recorded them yet (pre-manifest studies)."""
    out = {"dataset": None, "semantics": None}
    processed = Path(study_dir) / "processed_data"
    if not processed.is_dir():
        return out
    h5ads = sorted(processed.glob("*.h5ad"),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    if not h5ads:
        return out
    newest = h5ads[0]
    out["dataset"] = {"file": newest.name}
    try:
        import h5py
        with h5py.File(str(newest), "r") as h5:
            shape = None
            if "X" in h5:
                node = h5["X"]
                shape = node.attrs.get("shape")
                if shape is None and hasattr(node, "shape"):
                    shape = node.shape
            if shape is not None and len(shape) == 2:
                out["dataset"] = {"file": newest.name,
                                  "n_cells": int(shape[0]),
                                  "n_genes": int(shape[1])}

            semantics = {}
            obs_keys = set(h5["obs"].keys()) if "obs" in h5 else set()
            if "sample" in obs_keys:
                semantics["sample_column"] = "sample"
            if "cell_type" in obs_keys:
                semantics["cell_type_column"] = "cell_type"

            uns = h5.get("uns")
            if uns is not None:
                cond = uns.get("role_condition_col")
                if cond is not None:
                    try:
                        semantics["condition_column"] = _decode(cond[()])
                    except (TypeError, AttributeError):
                        pass
                rm = uns.get("role_map")
                if rm is not None and hasattr(rm, "keys"):
                    role_map = {}
                    for key in rm.keys():
                        try:
                            role_map[str(key)] = _decode(rm[key][()])
                        except (TypeError, AttributeError, KeyError):
                            continue
                    if role_map:
                        semantics["role_map"] = role_map

            # Roles baked into obs['_role'] without a recorded role_map
            if "role_map" not in semantics and "_role" in obs_keys:
                try:
                    cats = [_decode(c) for c in
                            h5["obs/_role/categories"][:]]
                    roles = [c for c in cats
                             if c in ("control", "disease", "exclude")]
                    if roles:
                        semantics.setdefault("condition_column", "_role")
                        semantics["role_map"] = {c: c for c in roles}
                except (KeyError, TypeError, OSError):
                    pass

            out["semantics"] = semantics or None
    except (OSError, KeyError, ValueError, ImportError):
        pass
    return out


def study_overview(study_dir: PathLike) -> dict:
    """Merge disk-derived asset flags with the cached manifest.

    Returns {'accession', 'status': <study_status dict>,
    'dataset': dict|None, 'semantics': dict|None}. Works for any study
    folder; every field degrades gracefully when nothing is recorded.
    When the manifest predates a field (studies processed before the
    manifest existed), dataset identity and semantics are derived from
    the newest processed h5ad's header instead.
    """
    study = Path(study_dir)
    manifest = read_manifest(study) or {}
    dataset = manifest.get("dataset")
    semantics = manifest.get("semantics")
    if dataset is None or semantics is None:
        peek = _peek_h5ad(study)
        dataset = dataset or peek["dataset"]
        semantics = semantics or peek["semantics"]
    return {
        "accession": study.name,
        "status": study_status(study),
        "dataset": dataset,
        "semantics": semantics,
    }
