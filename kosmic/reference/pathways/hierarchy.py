"""Auto-derive parent groupings for GPCR / ion-channel pathway collections.

Pattern-matches subfamily-level pathway names (e.g.
"Adhesion_G_protein_coupled_receptors_subfamily_A") to a smaller set of
parent categories ("Adhesion_GPCRs"). Returns None when no real
collapsing occurs, so callers can skip parent-level scoring for
collections without an obvious hierarchy.
"""

import re


_GPCR_PREFIX_RULES = (
    (re.compile(r'^Adhesion_G_protein_coupled_receptors_subfamily_'), 'Adhesion_GPCRs'),
    (re.compile(r'^Olfactory_receptors_family_'), 'Olfactory_receptors'),
)

_GPCR_EXACT = {
    'C_C_motif_chemokine_receptors':       'Chemokine_receptors',
    'C_X_C_motif_chemokine_receptors':     'Chemokine_receptors',
    'C_X_3_C_motif_chemokine_receptors':   'Chemokine_receptors',
    'X_C_motif_chemokine_receptors':       'Chemokine_receptors',
    'Atypical_chemokine_receptors':        'Chemokine_receptors',
    'Taste_1_receptors':                   'Taste_receptors',
    'Taste_2_receptors':                   'Taste_receptors',
    'G_protein_coupled_receptors_Class_A_orphans':   'Orphan_GPCRs',
    'G_protein_coupled_receptors_Class_C_orphans':   'Orphan_GPCRs',
    'G_protein_coupled_receptors_Class_F_frizzled':  'Frizzled_GPCRs',
}

_ION_CHANNEL_PREFIX_RULES = (
    ('Calcium_',  'Calcium_channels'),
    ('Potassium_', 'Potassium_channels'),
    ('Sodium_',   'Sodium_channels'),
    ('Chloride_', 'Chloride_channels'),
)

_ION_CHANNEL_REGEX_RULES = (
    (re.compile(r'^Glutamate_ionotropic'), 'Glutamate_ionotropic_receptors'),
)


def _classify(key):
    """Return the parent category for ``key``, or None if no rule matches."""
    for pattern, parent in _GPCR_PREFIX_RULES:
        if pattern.match(key):
            return parent
    if key in _GPCR_EXACT:
        return _GPCR_EXACT[key]
    for prefix, parent in _ION_CHANNEL_PREFIX_RULES:
        if key.startswith(prefix):
            return parent
    for pattern, parent in _ION_CHANNEL_REGEX_RULES:
        if pattern.match(key):
            return parent
    return None


def derive_hierarchy(pathways_dict):
    """Group pathway keys by their parent category.

    Returns a {parent_name: [child_keys]} mapping, or None if no parent
    ends up collapsing more than one child.
    """
    hierarchy: dict[str, list[str]] = {}
    for key in pathways_dict:
        parent = _classify(key)
        if parent is not None:
            hierarchy.setdefault(parent, []).append(key)
    if not any(len(children) > 1 for children in hierarchy.values()):
        return None
    return hierarchy


def build_parent_gene_sets(flat_pathways, hierarchy):
    """Merge subfamily gene lists into parent-level gene sets."""
    return {
        parent: sorted({g for child in children for g in flat_pathways.get(child, [])})
        for parent, children in hierarchy.items()
    }
