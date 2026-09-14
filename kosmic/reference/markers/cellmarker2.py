"""
CellMarker 2.0 Database Integration
====================================
Load and filter cell type markers from CellMarker 2.0 database.
Source: http://bio-bigdata.hrbmu.edu.cn/CellMarker/
"""

import pandas as pd
from importlib.resources import files
from pathlib import Path
from typing import Dict, List, Optional


_DEFAULT_DATA_PATH = files("kosmic.reference.marker_dbs") / "Cell_marker_Seq.xlsx"


class CellMarker2DB:
    """Load and query CellMarker 2.0 marker database."""

    def __init__(self, xlsx_path: Optional[str] = None):
        """Initialize with path to CellMarker 2.0 xlsx file."""
        if xlsx_path is None:
            xlsx_path = _DEFAULT_DATA_PATH

        self.xlsx_path = Path(xlsx_path)
        self._df = None
        self._load()

    def _load(self):
        """Load the xlsx file."""
        if not self.xlsx_path.exists():
            raise FileNotFoundError(f"CellMarker 2.0 file not found: {self.xlsx_path}")

        self._df = pd.read_excel(self.xlsx_path)
        # Drop rows with missing Symbol (no gene marker)
        self._df = self._df.dropna(subset=['Symbol'])

    @property
    def tissues(self) -> List[str]:
        """Get list of available tissues."""
        tissues = self._df['tissue_class'].dropna().unique().tolist()
        return sorted([t for t in tissues if t])

    def get_cell_types(self, tissue: Optional[str] = None, species: Optional[str] = None) -> List[str]:
        """Get cell types, optionally filtered by tissue and species."""
        df = self._df
        if tissue and tissue != "All":
            df = df[df['tissue_class'] == tissue]
        if species:
            df = df[df['species'] == species]
        cell_types = df['cell_name'].dropna().unique().tolist()
        return sorted([c for c in cell_types if c])

    def get_markers(
        self,
        cell_types: List[str],
        tissue: Optional[str] = None,
        species: str = 'Human',
        cancer_type: str = 'Normal'
    ) -> Dict[str, List[str]]:
        """
        Get marker genes for specified cell types.

        Parameters
        ----------
        cell_types : list
            List of cell type names to get markers for.
        tissue : str, optional
            Filter to markers from this tissue.
        species : str
            'Human' or 'Mouse'
        cancer_type : str
            'Normal' to exclude cancer markers, 'All' for everything.

        Returns
        -------
        dict
            Dictionary mapping cell_type -> list of gene symbols
        """
        df = self._df.copy()

        if species:
            df = df[df['species'] == species]
        if tissue and tissue != "All":
            df = df[df['tissue_class'] == tissue]
        if cancer_type and cancer_type != "All":
            df = df[df['cancer_type'] == cancer_type]

        df = df[df['cell_name'].isin(cell_types)]

        markers = {}
        for cell_type in df['cell_name'].unique():
            if not cell_type:
                continue
            ct_df = df[df['cell_name'] == cell_type]
            genes = ct_df['Symbol'].dropna().unique().tolist()
            if genes:
                markers[cell_type] = genes

        return markers


# Singleton instance
_instance = None


def get_cellmarker2() -> CellMarker2DB:
    """Get singleton instance of CellMarker 2.0 database."""
    global _instance
    if _instance is None:
        _instance = CellMarker2DB()
    return _instance


def get_tissues() -> List[str]:
    """Get list of available tissues."""
    return get_cellmarker2().tissues


def get_cell_types(tissue: Optional[str] = None, species: Optional[str] = None) -> List[str]:
    """Get list of cell types, optionally filtered by tissue and species."""
    return get_cellmarker2().get_cell_types(tissue, species)


def get_markers(
    cell_types: List[str],
    tissue: Optional[str] = None,
    species: str = 'Human',
    cancer_type: str = 'Normal'
) -> Dict[str, List[str]]:
    """Get marker genes for cell types."""
    return get_cellmarker2().get_markers(
        cell_types=cell_types,
        tissue=tissue,
        species=species,
        cancer_type=cancer_type
    )
