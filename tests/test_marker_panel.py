"""Ordered marker panels for the annotation sense check.

The dot plot only reads as a check when genes stay grouped by the cell
type they mark: annotated types on one axis, their markers on the other,
and a correct annotation shows a diagonal. Alphabetical order, or a gene
appearing twice because two types claim it, destroys that.
"""
import pytest

from kosmic.scrna.inspect.gene_group import marker_panel

MARKERS = {
    "Cardiomyocyte": ["TNNT2", "MYH7", "ACTC1"],
    "Fibroblast": ["DCN", "LUM", "COL1A1"],
    "Endothelial": ["PECAM1", "VWF"],
}


def test_genes_stay_grouped_by_their_cell_type():
    genes, _ = marker_panel(MARKERS, ["Cardiomyocyte", "Fibroblast"])
    assert genes == ["TNNT2", "MYH7", "ACTC1", "DCN", "LUM", "COL1A1"]


def test_requested_order_is_honoured():
    genes, _ = marker_panel(MARKERS, ["Fibroblast", "Cardiomyocyte"])
    assert genes[0] == "DCN"


def test_every_gene_knows_which_type_claimed_it():
    _, owner = marker_panel(MARKERS, ["Cardiomyocyte"])
    assert owner["TNNT2"] == "Cardiomyocyte"


def test_a_shared_gene_appears_once_under_the_earlier_type():
    """Two rows for one gene would double-count it in the plot."""
    markers = {"A": ["SHARED", "A1"], "B": ["SHARED", "B1"]}
    genes, owner = marker_panel(markers, ["A", "B"])
    assert genes == ["SHARED", "A1", "B1"]
    assert owner["SHARED"] == "A"


def test_max_per_type_caps_each_type_not_the_total():
    genes, _ = marker_panel(MARKERS, max_per_type=2)
    assert len(genes) == 6            # 2 each from three types
    assert genes[:2] == ["TNNT2", "MYH7"]      # Cardiomyocyte sorts first


def test_unknown_cell_types_are_skipped():
    genes, _ = marker_panel(MARKERS, ["Cardiomyocyte", "Not a type"])
    assert genes == ["TNNT2", "MYH7", "ACTC1"]


def test_a_type_with_no_markers_is_skipped():
    genes, _ = marker_panel({"A": [], "B": ["B1"]}, ["A", "B"])
    assert genes == ["B1"]


def test_defaults_to_every_type_sorted():
    genes, _ = marker_panel(MARKERS)
    assert genes[0] == "TNNT2"        # Cardiomyocyte sorts first


@pytest.mark.parametrize("cell_types", [[], ["Nope"]])
def test_nothing_selected_gives_an_empty_panel(cell_types):
    assert marker_panel(MARKERS, cell_types) == ([], {})


# --- matching your labels to a marker database ------------------------------

from kosmic.scrna.inspect.gene_group import match_marker_types  # noqa: E402

DB = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Macrophage", "T_cell",
      "B_cell", "Pericyte"]


def test_exact_names_match():
    assert match_marker_types(["Fibroblast"], DB)["Fibroblast"] == "Fibroblast"


def test_a_qualified_label_finds_its_base_type():
    """CellTypist says 'Ventricular Cardiomyocyte'; the database says
    'Cardiomyocyte'."""
    out = match_marker_types(["Ventricular Cardiomyocyte"], DB)
    assert out["Ventricular Cardiomyocyte"] == "Cardiomyocyte"


def test_a_prefixed_label_finds_its_base_type():
    out = match_marker_types(["LYVE1+ Macrophage"], DB)
    assert out["LYVE1+ Macrophage"] == "Macrophage"


def test_punctuation_and_case_are_ignored():
    assert match_marker_types(["CD4+ T Cell"], DB)["CD4+ T Cell"] == "T_cell"


def test_mast_cell_does_not_steal_t_cell_markers():
    """The bug this replaced: 'T_cell' normalises to 'tcell', which is a
    substring of 'mastcell', so Mast Cell silently got T-cell markers."""
    assert "Mast Cell" not in match_marker_types(["Mast Cell"], DB)


def test_an_unmatched_type_is_absent_not_guessed():
    """Better to say a type cannot be checked than to check it wrongly."""
    out = match_marker_types(["Neural Crest", "Adipocyte"], DB)
    assert out == {}


def test_the_most_specific_match_wins():
    out = match_marker_types(["Ventricular Cardiomyocyte"],
                             ["Myocyte", "Ventricular Cardiomyocyte"])
    assert out["Ventricular Cardiomyocyte"] == "Ventricular Cardiomyocyte"


def test_a_shared_base_type_may_serve_two_labels():
    out = match_marker_types(
        ["Endothelial Cell", "Lymphatic Endothelial Cell"], DB)
    assert set(out.values()) == {"Endothelial"}


def test_plurals_match_the_singular():
    """PanglaoDB writes 'Fibroblasts'; an annotation says 'Fibroblast'."""
    assert match_marker_types(["Fibroblast"], ["Fibroblasts"]) == {
        "Fibroblast": "Fibroblasts"}


def test_abbreviations_are_expanded():
    out = match_marker_types(["Smooth Muscle Cell", "NK Cell"],
                             ["SMC", "Natural killer cell"])
    assert out == {"Smooth Muscle Cell": "SMC",
                   "NK Cell": "Natural killer cell"}


def test_the_cardiac_panel_covers_a_heart_atlas():
    """Every broad label the CellTypist heart model gives has markers."""
    from kosmic.reference.marker_panels import builtin_panels
    labels = ["Adipocyte", "Endocardial Cell", "Endothelial Cell",
              "Fibroblast", "Lymphatic Endothelial Cell", "Macrophage",
              "Mast Cell", "Mesothelial Cell", "Neural Crest", "Pericyte",
              "Smooth Muscle Cell", "T Cell", "Ventricular Cardiomyocyte"]
    out = match_marker_types(labels, builtin_panels()["Cardiac (human, curated)"])
    assert set(out) == set(labels)
    assert out["Lymphatic Endothelial Cell"] == "Lymphatic endothelial"


# --- the Marker Check dot plot ------------------------------------------------

@pytest.fixture(scope="module")
def qapp():
    import os
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6.QtWidgets import QApplication
    yield QApplication.instance() or QApplication([])


def test_dot_plot_orders_rows_as_a_diagonal_and_paints(qapp):
    """Rows follow the panel order by expression, as in the exported
    figure, and the widget paints without error."""
    import pandas as pd

    from kosmic.gui.shared.plots import GeneGroupDotPlot
    mean = pd.DataFrame({'Fibroblast': [0.1, 2.0], 'Ventricular CM': [3.0, 0.2]},
                        index=['TNNT2', 'DCN'])
    pct = mean.gt(1).astype(float)
    w = GeneGroupDotPlot()
    w.resize(600, 400)
    w.set_data(mean, pct, {'TNNT2': 'Cardiomyocyte', 'DCN': 'Fibroblast'})
    assert w._rows == ['Ventricular CM', 'Fibroblast']
    assert not w.grab().isNull()
