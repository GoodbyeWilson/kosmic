"""A reference discovered in kosmic/reference/atlases/ can be loaded by the
name the Annotate dropdown shows, not only by its file path.

The built-in catalogue is empty, so name resolution has to go through
discovery; before this, the dropdown offered "HeartMap LV broad (Datar
2026)" and load_reference raised "Reference not found" for it.
"""
from kosmic.scrna.annotate.reference import (
    list_available_references, load_reference,
)


def test_every_discovered_reference_loads_by_name():
    available = list_available_references()
    assert available, "the two shipped LV references should be discovered"
    for name in available:
        ref = load_reference(name)
        assert ref['name'] == name
        assert ref['n_cell_types'] > 0 and ref['n_genes'] > 0


def test_unknown_name_lists_what_is_available():
    import pytest
    with pytest.raises(FileNotFoundError, match="Available:"):
        load_reference("no such reference")
