"""Figure Export workspace -- a top-level sibling to scRNA / DE / Meta.

The workspace presents matplotlib figure generators in a single place. Each
page renders one figure type from data owned by the scRNA / DE workspaces.
"""
from kosmic.gui.figure_export.workspace import FigureExportWorkspace

__all__ = ["FigureExportWorkspace"]
