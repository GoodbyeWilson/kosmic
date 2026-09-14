"""
PanglaoDB Marker Database Integration
=====================================
Load and filter cell type markers from PanglaoDB database.
"""

import pandas as pd
from importlib.resources import files
from pathlib import Path
from typing import Dict, List, Optional


_DEFAULT_DATA_PATH = files("kosmic.reference.marker_dbs") / "PanglaoDB_markers.tsv"


class PanglaoDBMarkers:
    """Load and query PanglaoDB marker database."""

    def __init__(self, tsv_path: Optional[str] = None):
        """Initialize with path to PanglaoDB TSV file."""
        if tsv_path is None:
            tsv_path = _DEFAULT_DATA_PATH

        self.tsv_path = Path(tsv_path)
        self._df = None
        self._load()

    def _load(self):
        """Load the TSV file."""
        if not self.tsv_path.exists():
            raise FileNotFoundError(f"PanglaoDB markers file not found: {self.tsv_path}")

        self._df = pd.read_csv(self.tsv_path, sep='\t')
        # Clean up column names
        self._df.columns = [c.strip().replace(' ', '_') for c in self._df.columns]

    @property
    def organs(self) -> List[str]:
        """Get list of available organs/tissues."""
        organs = self._df['organ'].dropna().unique().tolist()
        organs = [o for o in organs if o and o != 'NA' and o != 'organ']
        return sorted(organs)

    @property
    def cell_types(self) -> List[str]:
        """Get list of all cell types."""
        cell_types = self._df['cell_type'].dropna().unique().tolist()
        cell_types = [c for c in cell_types if c and c != 'cell type']
        return sorted(cell_types)

    def get_cell_types_for_organ(self, organ: str) -> List[str]:
        """Get cell types associated with a specific organ."""
        mask = self._df['organ'] == organ
        cell_types = self._df.loc[mask, 'cell_type'].dropna().unique().tolist()
        return sorted([c for c in cell_types if c and c != 'cell type'])

    def get_markers(
        self,
        cell_types: Optional[List[str]] = None,
        organ: Optional[str] = None,
        species: str = 'both',  # 'human', 'mouse', or 'both'
        min_specificity: float = 0.0,
        canonical_only: bool = False
    ) -> Dict[str, List[str]]:
        """
        Get marker genes for specified cell types.

        Parameters
        ----------
        cell_types : list, optional
            List of cell types to get markers for. If None, uses organ filter.
        organ : str, optional
            Filter to cell types from this organ/tissue.
        species : str
            'human', 'mouse', or 'both'
        min_specificity : float
            Minimum specificity score (0-1) to include marker
        canonical_only : bool
            If True, only include canonical markers

        Returns
        -------
        dict
            Dictionary mapping cell type -> list of marker genes
        """
        df = self._df.copy()

        # Filter by organ if specified
        if organ and organ != "All":
            df = df[df['organ'] == organ]

        # Filter by cell types if specified
        if cell_types:
            df = df[df['cell_type'].isin(cell_types)]

        # Filter by species
        if species == 'human':
            df = df[df['species'].str.contains('Hs', na=False)]
            spec_col = 'specificity_human'
        elif species == 'mouse':
            df = df[df['species'].str.contains('Mm', na=False)]
            spec_col = 'specificity_mouse'
        else:
            # Use average of human and mouse specificity
            df['avg_specificity'] = df[['specificity_human', 'specificity_mouse']].mean(axis=1)
            spec_col = 'avg_specificity'

        # Filter by specificity
        if min_specificity > 0:
            df = df[df[spec_col] >= min_specificity]

        # Filter canonical markers only
        if canonical_only:
            df = df[df['canonical_marker'] == 1]

        # Build marker dictionary
        markers = {}
        for cell_type in df['cell_type'].unique():
            if not cell_type or cell_type == 'cell type':
                continue
            ct_df = df[df['cell_type'] == cell_type]
            # Sort by specificity (highest first) and take gene symbols
            ct_df = ct_df.sort_values(spec_col, ascending=False)
            genes = ct_df['official_gene_symbol'].dropna().unique().tolist()
            if genes:
                markers[cell_type] = genes

        return markers

# Singleton instance for easy access
_instance = None

def get_panglaodb() -> PanglaoDBMarkers:
    """Get singleton instance of PanglaoDB markers."""
    global _instance
    if _instance is None:
        _instance = PanglaoDBMarkers()
    return _instance


def get_organs() -> List[str]:
    """Get list of available organs."""
    return get_panglaodb().organs


def get_cell_types(organ: Optional[str] = None) -> List[str]:
    """Get list of cell types, optionally filtered by organ."""
    if organ:
        return get_panglaodb().get_cell_types_for_organ(organ)
    return get_panglaodb().cell_types


def get_markers(
    cell_types: List[str],
    organ: Optional[str] = None,
    species: str = 'both',
    min_specificity: float = 0.0,
    canonical_only: bool = False
) -> Dict[str, List[str]]:
    """Get marker genes for cell types."""
    return get_panglaodb().get_markers(
        cell_types=cell_types,
        organ=organ,
        species=species,
        min_specificity=min_specificity,
        canonical_only=canonical_only
    )
