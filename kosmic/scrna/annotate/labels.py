# Writing a cell-type label into obs.
#
# obs['cell_type'] is usually categorical, and pandas refuses to write a
# value that is not already one of its categories ("Cannot setitem on a
# Categorical with a new category"). A manual relabel to a type the
# automated pass never used -- most often 'Unknown' -- therefore failed
# silently in the GUI. 'set_cell_type' adds the category first.
from __future__ import annotations

import numpy as np
import pandas as pd


def set_cell_type(obs: pd.DataFrame, mask, label: str,
                  col: str = 'cell_type') -> int:
    """Set ``obs[col]`` to ``label`` for the rows in ``mask``; return the count.

    Works whether the column is categorical, object or absent (created as
    categorical). A categorical column gains the label as a category
    when it lacks it; categories are never removed.
    """
    mask = np.asarray(mask, dtype=bool)
    if col not in obs.columns:
        obs[col] = pd.Categorical([None] * len(obs))
    s = obs[col]
    if isinstance(s.dtype, pd.CategoricalDtype):
        if label not in s.cat.categories:
            obs[col] = s.cat.add_categories([label])
        obs.loc[mask, col] = label
    else:
        obs.loc[mask, col] = label
    return int(mask.sum())
