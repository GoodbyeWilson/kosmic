# Per-dataset condition role assignment.
#
# Every value in the chosen condition column is mapped to
# 'control', 'disease', or 'exclude' via a 'role_map' dict
# saved in 'analysis_settings.json'. This module is the boundary
# between that role model and binary callers that still expect
# 'disease_label' / 'control_label' strings.
#
# Two helpers live here:
#
# - 'parse_role_map' -- reads 'parameters' from a project's
#   'analysis_settings.json'. Returns the 'role_map' field if
#   present; otherwise synthesises one from
#   'disease_group' / 'control_group' fields.
# - 'roles_to_labels' -- picks one 'disease_label' and one
#   'control_label' from a 'role_map' for binary callers. When the
#   role map specifies multiple disease values (e.g. DCM + ICM) the
#   binary projection picks the first by sort order.
#
# Lookup is case-insensitive. The saved role map preserves whatever
# casing exists in the data on first save (round-trip stable).
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


# Legal role values stored in role_map.
_ROLE_DISEASE = 'disease'
_ROLE_CONTROL = 'control'
_ROLE_EXCLUDE = 'exclude'
_VALID_ROLES = frozenset((_ROLE_DISEASE, _ROLE_CONTROL, _ROLE_EXCLUDE))


def parse_role_map(parameters: dict) -> Optional[dict[str, str]]:
    """Extract a role map from a settings 'parameters' dict.

    Reads two formats:

    1. **role_map:** 'parameters['role_map']' is a 'dict' mapping
       condition-column values to one of 'control', 'disease',
       'exclude'. Used directly.
    2. **disease_group/control_group:** single strings synthesised into
       '{disease_group: 'disease', control_group: 'control'}'. Other
       values in the condition column are not added -- they're treated
       as excluded (binary contrast).

    Returns 'None' if neither format is present.
    """
    if not isinstance(parameters, dict):
        return None

    role_map = parameters.get('role_map')
    if isinstance(role_map, dict) and role_map:
        # Validate roles; silently drop any unknown role values.
        cleaned = {
            str(k): v for k, v in role_map.items()
            if isinstance(v, str) and v in _VALID_ROLES
        }
        return cleaned or None

    disease_group = parameters.get('disease_group')
    control_group = parameters.get('control_group')
    if disease_group is not None and control_group is not None:
        return {
            str(disease_group): _ROLE_DISEASE,
            str(control_group): _ROLE_CONTROL,
        }

    return None


def roles_to_labels(role_map: dict[str, str]) -> tuple[Optional[str], Optional[str]]:
    """Pick one 'disease_label' and one 'control_label' from a role map.

    Binary projection for callers that consume single-string labels
    rather than the full role map. Picks the first value (sorted)
    assigned to each role. Returns '(None, None)' if either role is
    unrepresented in the map.

    Examples
    --------
    >>> roles_to_labels({'NF': 'control', 'DCM': 'disease'})
    ('DCM', 'NF')
    >>> roles_to_labels({'NF': 'control', 'DCM': 'disease', 'ICM': 'disease'})
    ('DCM', 'NF')
    >>> roles_to_labels({'NF': 'control'})  # no disease assignment
    (None, 'NF')
    """
    diseases = sorted(k for k, r in role_map.items() if r == _ROLE_DISEASE)
    controls = sorted(k for k, r in role_map.items() if r == _ROLE_CONTROL)
    disease_label = diseases[0] if diseases else None
    control_label = controls[0] if controls else None
    return disease_label, control_label


def de_readiness(adata) -> Tuple[bool, str]:
    """Whether 'adata' is ready to hand off from scRNA to DE.

    Ready means: normalised counts in '.X' (same '< 20' heuristic used
    elsewhere in the scRNA workspace for this check) and a Control/Disease
    role assignment present, either as a 'role_map' in 'adata.uns' or as
    an already-resolved '_role' column in 'adata.obs'.

    Returns '(ready, reason)'; 'reason' is '' when ready, otherwise a
    one-line, user-facing explanation of what's missing.
    """
    if adata is None:
        return False, "No dataset loaded."

    try:
        is_normalised = float(adata.X.max()) < 20
    except (TypeError, ValueError, AttributeError):
        is_normalised = True
    if not is_normalised:
        return False, (
            "Data isn't normalised yet -- finish Quality Control in "
            "scRNA before running DE.")

    role_map = adata.uns.get('role_map') if hasattr(adata, 'uns') else None
    has_role_map = (
        isinstance(role_map, dict)
        and any(v == _ROLE_CONTROL for v in role_map.values())
        and any(v == _ROLE_DISEASE for v in role_map.values())
    )

    has_role_col = False
    if hasattr(adata, 'obs') and '_role' in adata.obs.columns:
        present = set(adata.obs['_role'].astype(str).str.lower().unique())
        has_role_col = _ROLE_CONTROL in present and _ROLE_DISEASE in present

    if not (has_role_map or has_role_col):
        return False, (
            "No Control/Disease roles assigned yet -- set them in "
            "scRNA → Inspect → Setup.")

    return True, ""


def resolve_roles(conditions, role_map: dict[str, str]) -> Tuple[np.ndarray, np.ndarray]:
    """Resolve per-sample condition labels to disease / two-group masks.

    Used by inner-loop callers ('fast_deseq2_perm', 'fast_ttest_de',
    'pydeseq2_reference_de') so they no longer need to know about
    literal 'disease_label' / 'control_label' strings -- they
    receive pre-computed boolean masks instead.

    Lookup is case-insensitive. Values in 'conditions' that are not
    present in 'role_map', or are mapped to 'exclude', get
    'False' in both returned masks.

    Parameters
    ----------
    conditions : array-like of str
        Per-sample condition values, e.g. 'adata.obs[condition_col].values'.
    role_map : dict
        Maps condition-column values to one of 'control',
        'disease', 'exclude'. Keys are matched case-insensitively.

    Returns
    -------
    is_disease : ndarray of bool, shape (n_samples,)
        True iff the sample's role is 'disease'.
    two_group_mask : ndarray of bool, shape (n_samples,)
        True iff the sample's role is 'disease' or 'control'
        (i.e. not excluded and not unmapped).
    """
    conditions_lc = np.array([str(c).lower() for c in np.asarray(conditions)])

    disease_keys = [str(k).lower() for k, v in role_map.items()
                    if v == _ROLE_DISEASE]
    control_keys = [str(k).lower() for k, v in role_map.items()
                    if v == _ROLE_CONTROL]

    is_disease = np.isin(conditions_lc, disease_keys)
    is_control = np.isin(conditions_lc, control_keys)
    two_group_mask = is_disease | is_control

    return is_disease, two_group_mask



# ---------------------------------------------------------------------------
# Honouring 'exclude' upstream of DE
# ---------------------------------------------------------------------------
#
# 'exclude' started as a DE-only idea: the contrast skips those samples.
# Everything before DE ignored it, so an excluded arm still drove HVG
# selection, PCA, the neighbour graph, the clusters and the UMAP. On
# GSE292067 that was 67,707 doxorubicin-cardiomyopathy cells -- 44% of the
# dataset -- shaping an embedding used to answer a DCM-vs-donor question.
#
# The cells are not deleted (that is what the Inspect tab's Filter Dataset
# does, and it is irreversible). Instead a step computes on the included
# cells and scatters its results back, leaving excluded cells present but
# unlabelled: no cluster, no coordinates, no influence.

def included_mask(adata) -> Optional[np.ndarray]:
    """Boolean mask of cells not marked ``_role == 'exclude'``.

    Returns ``None`` when there is nothing to honour -- no ``_role``
    column, or no cell excluded by it -- so callers can skip the
    subset-and-scatter path entirely and keep the fast, simple route.
    Also returns ``None`` when *every* cell is excluded, since running an
    analysis on nothing is never what was meant.
    """
    if not hasattr(adata, 'obs') or '_role' not in adata.obs.columns:
        return None
    keep = (adata.obs['_role'].astype(str).str.lower()
            != _ROLE_EXCLUDE).to_numpy()
    if keep.all() or not keep.any():
        return None
    return keep


def scatter_results(full, sub, mask, obsm_keys=(), obs_cols=(),
                    uns_keys=(), obsp_keys=()):
    """Write a subset analysis back onto the full object.

    Excluded rows get NaN in ``obsm`` and an unset category in ``obs``,
    which is the honest representation: the analysis has no opinion about
    a cell it never saw. ``obsp`` matrices (the neighbour graph) are not
    scattered -- they are cell-by-cell and meaningless off the subset --
    so they are dropped from ``full`` rather than left stale.

    Parameters
    ----------
    full : AnnData
        The object to write onto, unchanged in cell count.
    sub : AnnData
        The analysed subset, in the same order as ``full[mask]``.
    mask : ndarray of bool
        Which rows of ``full`` ``sub`` corresponds to.
    """
    import pandas as pd

    for key in obsm_keys:
        if key not in sub.obsm:
            continue
        out = np.full((full.n_obs, sub.obsm[key].shape[1]), np.nan,
                      dtype=np.float32)
        out[mask] = np.asarray(sub.obsm[key], dtype=np.float32)
        full.obsm[key] = out

    for col in obs_cols:
        if col not in sub.obs.columns:
            continue
        values = pd.Series(pd.NA, index=full.obs_names, dtype='object')
        values.iloc[np.flatnonzero(mask)] = sub.obs[col].astype(str).values
        full.obs[col] = pd.Categorical(values)

    for key in uns_keys:
        if key in sub.uns:
            full.uns[key] = sub.uns[key]

    for key in obsp_keys:
        full.obsp.pop(key, None)

    return full


def cluster_labels(adata, cluster_col: str = 'leiden') -> list:
    """Cluster labels that actually have cells, sorted numerically first.

    Cells left out of the embedding (see :func:`included_mask`) carry no
    cluster, so the column holds nulls. ``Series.unique()`` returns NaN as
    a value but ``series == NaN`` matches nothing, so the naive loop hands
    downstream code an empty group -- ``value_counts().index[0]`` then
    raises ``IndexError: index 0 is out of bounds``. A categorical can
    also retain a category after its last cell is filtered away, which
    fails the same way.

    Returns ``[]`` when the column is absent.
    """
    if not hasattr(adata, 'obs') or cluster_col not in adata.obs.columns:
        return []
    counts = adata.obs[cluster_col].value_counts(dropna=True)
    present = [label for label, n in counts.items() if n > 0]
    return sorted(present,
                  key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))


# Pattern dictionaries the Inspect tab uses to pre-fill role guesses.
# Public so the GUI's tiny "guess role from value" wrapper can import
# them; centralised here so they aren't duplicated across pages.
# Biased toward heart-failure / immunology conventions but cover most
# disease-vs-healthy designs. Always override-able in the Inspect-tab UI.
CONTROL_VALUE_PATTERNS = (
    'control', 'ctrl', 'normal', 'healthy', 'sham', 'wt',
    'wildtype', 'wild_type', 'wild-type', 'untreated', 'donor',
    'nf', 'nondiseased', 'non_failing', 'non-failing',
)
DISEASE_VALUE_PATTERNS = (
    'disease', 'case', 'patient', 'treated', 'mutant', 'ko',
    'knockout', 'dcm', 'icm', 'hcm', 'hf', 'hfref', 'hfpef',
    'tumor', 'cancer', 'failing',
)
