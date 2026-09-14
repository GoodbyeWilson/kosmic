"""
Gene Set File Formats
======================
Import/export gene sets in GMT, JSON, and CSV formats.
Pure Python — no GUI imports.
"""

import csv
import json
from pathlib import Path


def load_gmt(path):
    """Load gene sets from a GMT (Gene Matrix Transposed) file.

    GMT format: PathwayName<TAB>Description<TAB>GENE1<TAB>GENE2<TAB>...

    Parameters
    ----------
    path : str or Path
        Path to .gmt file.

    Returns
    -------
    dict
        {pathway_name: [gene_list]}.
    """
    pathways = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            # parts[1] is description (often a URL or 'na'), skip it
            genes = [g.strip() for g in parts[2:] if g.strip()]
            pathways[name] = genes
    return pathways


def save_gmt(pathways, path, descriptions=None):
    """Save gene sets to a GMT file.

    Parameters
    ----------
    pathways : dict
        {pathway_name: [gene_list]}.
    path : str or Path
        Output .gmt file path.
    descriptions : dict, optional
        {pathway_name: description_string}. Defaults to 'na'.
    """
    descriptions = descriptions or {}
    with open(path, 'w', encoding='utf-8', newline='') as f:
        for name, genes in pathways.items():
            desc = descriptions.get(name, 'na')
            line = '\t'.join([name, desc] + genes)
            f.write(line + '\n')


def load_json(path):
    """Load gene sets from a JSON file.

    Expected format: {"pathway_name": ["GENE1", "GENE2", ...], ...}

    Parameters
    ----------
    path : str or Path
        Path to .json file.

    Returns
    -------
    dict
        {pathway_name: [gene_list]}.
    """
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("JSON file must contain a top-level dictionary")
    return data


def save_json(pathways, path, indent=2):
    """Save gene sets to a JSON file.

    Parameters
    ----------
    pathways : dict
        {pathway_name: [gene_list]}.
    path : str or Path
        Output .json file path.
    indent : int
        JSON indentation level.
    """
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(pathways, f, indent=indent)


def load_csv(path):
    """Load gene sets from a CSV file.

    Expected format: two columns — pathway_name, gene_symbol.
    Each row maps one gene to one pathway. A pathway appears on multiple rows.

    Parameters
    ----------
    path : str or Path
        Path to .csv file.

    Returns
    -------
    dict
        {pathway_name: [gene_list]}.
    """
    pathways = {}
    with open(path, 'r', encoding='utf-8', newline='') as f:
        reader = csv.reader(f)
        next(reader, None)  # discard header row
        for row in reader:
            if len(row) < 2:
                continue
            name = row[0].strip()
            gene = row[1].strip()
            if name and gene:
                pathways.setdefault(name, []).append(gene)
    return pathways


def load_file(path):
    """Auto-detect format and load gene sets from a file.

    Supported extensions: .gmt, .json, .csv

    Parameters
    ----------
    path : str or Path
        Path to gene set file.

    Returns
    -------
    dict
        {pathway_name: [gene_list]}.

    Raises
    ------
    ValueError
        If file extension is not recognised.
    """
    path = Path(path)
    ext = path.suffix.lower()
    if ext == '.gmt':
        return load_gmt(path)
    elif ext == '.json':
        return load_json(path)
    elif ext == '.csv':
        return load_csv(path)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .gmt, .json, or .csv")
