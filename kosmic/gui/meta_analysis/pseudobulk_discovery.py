# Pseudobulk-path discovery for the meta-analysis pages.
#
# Both the consensus-MA page and the methods-comparison page need to
# locate '*_pseudobulk.csv' files alongside (or under) the user's DE
# result CSVs. Discovery is the same algorithm in both places, so it
# lives here as a free function rather than as a method on either page.
#
# Caching: results are cached in-memory keyed on the dataset-path tuple,
# and persisted across sessions via 'QSettings' so re-opening the app
# on a slow/SMB share doesn't re-walk the filesystem.
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Callable, Optional, Sequence

from PyQt6.QtCore import QSettings

from kosmic.paths import pseudobulk_path, study_dir

_QSETTINGS_KEY = "consensus/pb_paths_cache_v2"

# Process-local cache, keyed on tuple-of-dataset-paths. Survives the
# duration of a single app run; the QSettings persistence below covers
# cross-restart caching.
_memory_cache: dict[tuple, list[str]] = {}


def clear_memory_cache() -> None:
    """
    Drop the in-process discovery cache.

    Call this when the meaning of "the same dataset path tuple" has
    changed -- for example, when the user picks a different DE method
    so the project-folder fallback layout shifts. The persisted
    QSettings cache is left untouched (it's keyed on absolute paths).
    """
    _memory_cache.clear()


def _emit(log_cb: Optional[Callable[[str], None]], msg: str) -> None:
    if log_cb is not None:
        log_cb(msg)


def _load_persisted_cache() -> dict:
    settings = QSettings("KOSMIC", "KOSMIC")
    try:
        raw = settings.value(_QSETTINGS_KEY, "{}")
        return json.loads(raw) if raw else {}
    except Exception:
        return {}


def _save_persisted_cache(cache: dict) -> None:
    settings = QSettings("KOSMIC", "KOSMIC")
    try:
        settings.setValue(_QSETTINGS_KEY, json.dumps(cache))
    except Exception:
        pass


def _try_project_folder_fallback(
    project_folder: Optional[str], name: str,
    log_cb: Optional[Callable[[str], None]],
) -> Optional[str]:
    """Try to find '<name>_pseudobulk.csv' under 'project_folder'."""
    if not project_folder:
        return None
    candidates = [
        # Standard layout: <project>/<dataset>/results/de_analysis/statistics/<dataset>_pseudobulk.csv
        str(pseudobulk_path(study_dir(project_folder, name), name)),
        # Maybe project_folder IS the data root one level up
        str(pseudobulk_path(study_dir(Path(project_folder).parent, name), name)),
        # Maybe it's directly inside project_folder
        os.path.join(str(project_folder), f'{name}_pseudobulk.csv'),
    ]
    for candidate in candidates:
        try:
            if os.path.exists(candidate):
                return candidate
        except OSError:
            continue
    _emit(log_cb, f"    fallback tried: {candidates[0]}")
    return None


def find_pseudobulk_paths(
    datasets: Sequence,
    project_folder: Optional[str],
    *,
    log_cb: Optional[Callable[[str], None]] = None,
) -> list[str]:
    """
    Return one '*_pseudobulk.csv' path per entry in 'datasets'.

    Each dataset dict must have 'path' (DE-stats CSV) and 'name' (study
    label). Falls back to 'project_folder' when the pseudobulk isn't
    next to its DE-stats CSV. Returns '[]' if any study can't be
    resolved -- the consensus pipeline needs a complete set.
    """
    cache_key = tuple(
        d.get('path', '') for d in datasets if isinstance(d, dict))

    if cache_key in _memory_cache:
        return list(_memory_cache[cache_key])

    persisted = _load_persisted_cache()
    cache_key_str = "|".join(cache_key)
    if cache_key_str in persisted:
        cached_paths = persisted[cache_key_str]
        # Trust the cache without re-validating with os.path.exists,
        # because exists() returns False on a dead SMB connection
        # which would defeat the cache. load_pseudobulk has its own
        # retry/re-auth logic for the actual reads.
        _emit(log_cb,
              f"Pseudobulk: using {len(cached_paths)} cached paths "
              "from previous session")
        _memory_cache[cache_key] = list(cached_paths)
        return list(cached_paths)

    _emit(log_cb, f"  project_folder = {project_folder}")

    pb_paths: list[str] = []
    for d in datasets:
        if not isinstance(d, dict) or 'path' not in d:
            continue
        name = d.get('name', '?')

        # Prefer project-folder layout; fall back to scanning the
        # original DE output directory.
        fallback = _try_project_folder_fallback(project_folder, name, log_cb)
        if fallback:
            pb_paths.append(fallback)
            continue

        stats_dir = str(Path(str(d['path'])).parent)
        try:
            files = os.listdir(stats_dir)
        except OSError:
            _emit(log_cb,
                  f"  pseudobulk scan failed for {name}: "
                  f"cannot list {stats_dir} and no fallback found")
            continue

        matches = [f for f in files if f.endswith('_pseudobulk.csv')]
        if matches:
            pb_paths.append(os.path.join(stats_dir, matches[0]))
        else:
            _emit(log_cb, f"  no pseudobulk file in {stats_dir}")

    if len(pb_paths) != len(datasets):
        return []

    _memory_cache[cache_key] = list(pb_paths)
    persisted[cache_key_str] = list(pb_paths)
    _save_persisted_cache(persisted)
    return pb_paths
