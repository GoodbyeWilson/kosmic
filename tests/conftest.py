"""Test-wide settings.

pandas copy-on-write is the default from pandas 3, where `to_numpy()` on a
column can return a read-only view; code that writes into such an array
passes with pandas 2 and fails in CI. Switching it on under pandas 2 makes
local runs behave the same way.
"""
import pandas as pd

if int(pd.__version__.split('.')[0]) < 3:
    pd.options.mode.copy_on_write = True
