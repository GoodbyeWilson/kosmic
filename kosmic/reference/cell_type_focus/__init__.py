# Loader for cell-type focus presets used by the Annotate tab to
# pre-select tissue-relevant cell types (e.g. cardiac, stroke).
#
# Bundled presets live next to this file as 'X.json'. Users add their
# own by dropping a JSON in '~/.kosmic/cell_type_focus/'; user files
# shadow bundled ones with the same 'name'.
#
# Each preset has a 'name' (str) and any combination of 'panglaodb',
# 'cellmarker2_human', 'cellmarker2_mouse' (lists of cell-type strings).
from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path

USER_PRESET_DIR = Path.home() / '.kosmic' / 'cell_type_focus'


def load_focus_presets() -> dict[str, dict]:
    """
    Return {preset_name: preset_dict} from bundled + user dirs.

    User-dir presets shadow bundled ones with the same name. Malformed
    JSON in the user dir is silently skipped (a single bad file
    shouldn't break the whole presets list).
    """
    presets: dict[str, dict] = {}

    pkg = files('kosmic.reference.cell_type_focus')
    for entry in pkg.iterdir():
        if entry.name.endswith('.json'):
            with entry.open('rb') as f:
                preset = json.loads(f.read().decode('utf-8'))
            presets[preset['name']] = preset

    if USER_PRESET_DIR.exists():
        for path in sorted(USER_PRESET_DIR.glob('*.json')):
            try:
                preset = json.loads(path.read_text(encoding='utf-8'))
                presets[preset['name']] = preset
            except (OSError, json.JSONDecodeError, KeyError):
                continue

    return presets
