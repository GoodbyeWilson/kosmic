"""Test-wide settings.

pandas copy-on-write is the default from pandas 3, where `to_numpy()` on a
column can return a read-only view; code that writes into such an array
passes with pandas 2 and fails in CI. Switching it on under pandas 2 makes
local runs behave the same way.

Every test waits for the figure previews it started before the next one
runs, so a render from one test does not finish (and draw, or report)
during another.
"""
import sys

import pandas as pd
import pytest

if int(pd.__version__.split('.')[0]) < 3:
    pd.options.mode.copy_on_write = True


@pytest.fixture(autouse=True)
def _wait_for_figure_renders():
    yield
    # Only when a test has loaded the figure pages: importing them here
    # would bring Qt into every test.
    base = sys.modules.get('kosmic.gui.figure_export.pages._base')
    if base is not None:
        base.wait_for_renders()
