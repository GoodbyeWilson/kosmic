# The question asked when an h5ad holds more than one matrix of raw counts.
#
# A Broad Single Cell Portal file, for instance, keeps CellBender-corrected
# counts in X and the uncorrected Cell Ranger counts in a layer. Which one to
# analyse is a scientific choice, not a technical one, so both import routes
# -- the Project page and Load Data -- ask with this one dialog, and the
# import records the answer in the study's 'load' stage.
from __future__ import annotations

from pathlib import Path
from typing import Callable, Optional


def candidates_summary(described: list) -> list:
    """The matrices considered, as recorded with the import."""
    return [{'slot': d['slot'], 'verdict': d['verdict'], 'max': round(d['max'], 3),
             **({'median_counts_per_cell': round(d['median_total'], 1)}
                if d.get('median_total') is not None else {})}
            for d in described]


def choose_counts_matrix(parent, raw_path, log: Optional[Callable[[str], None]] = None):
    """Return '(slot, described)' for an h5ad, asking when it is a real choice.

    'slot' is None to let the import pick (one obvious candidate), the slot
    token the user chose, or False if they cancelled. 'described' is the
    description of every matrix in the file. Inspection loads no matrices,
    so this is instant even on a 13 GB file.
    """
    from PyQt6.QtWidgets import QInputDialog

    from kosmic.scrna.load.matrices import RAW_COUNTS, describe_count_matrices_path

    say = log or (lambda _msg: None)
    try:
        described = describe_count_matrices_path(raw_path)
    except Exception as e:  # noqa: BLE001 -- the import reports a bad file
        say(f"(could not inspect matrices: {e})")
        return None, []

    counts = [d for d in described if d['verdict'] == RAW_COUNTS]
    if len(counts) < 2:
        return None, described        # nothing to decide; the import picks

    # Say what actually differs rather than guessing at intent: the matrix
    # with fewer stored values, or fewer counts per cell, has had something
    # removed, usually ambient correction or extra filtering. Which is which
    # is the depositor's business, so point at their documentation.
    fullest = max(counts, key=lambda d: d['n_stored'])

    def line(d):
        total = d.get('median_total')
        per_cell = f", median {total:,.0f} counts per cell" if total is not None else ""
        return f"  {d['label']:26s} {d['n_stored']:>15,} stored values{per_cell}"

    deltas = [
        f"{d['label']} holds {100 * (fullest['n_stored'] - d['n_stored']) / fullest['n_stored']:.1f}%"
        f" fewer values than {fullest['label']}"
        for d in counts
        if d is not fullest and d['n_stored'] < fullest['n_stored']
    ]
    detail = ("\n\n" + "; ".join(deltas)
              + " — so it has had something removed, typically ambient-RNA "
                "correction (for example CellBender) or extra filtering. Which "
                "is which is the depositor's choice: check the dataset's README "
                "or the paper's methods.") if deltas else ""

    items = [f"{d['label']}  —  max {d['max']:,.0f}" for d in counts]
    default = next((i for i, d in enumerate(counts) if d['slot'] == 'X'), 0)
    choice, ok = QInputDialog.getItem(
        parent, "Which counts matrix?",
        f"{Path(raw_path).name} holds more than one matrix of raw integer "
        "counts.\n\n" + "\n".join(line(d) for d in counts) + detail
        + "\n\nPick the matrix to analyse from. If ambient RNA is to be removed "
          "in KOSMIC (Decontaminate), import the uncorrected counts so it is "
          "removed once. The choice, and the options not chosen, are recorded "
          "in the study's methods.",
        items, default, False)
    if not ok:
        say("Import cancelled at matrix selection.")
        return False, described
    return counts[items.index(choice)]['slot'], described
