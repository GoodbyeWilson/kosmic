"""A master built from the studies' DecontX-corrected counts (ADR-011)."""
import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from kosmic import provenance
from kosmic.combine import concat_studies
from kosmic.scrna.qc.decontx import DECONTX_LAYER

GENES = ["TTN", "NPPA", "DCN"]


def _study(path, decontx=True, n_obs=3):
    """Original counts of 5 and DecontX counts of 3 in every cell."""
    path.parent.mkdir(parents=True, exist_ok=True)
    a = ad.AnnData(
        X=csr_matrix(np.full((n_obs, len(GENES)), 0.5, dtype=np.float32)),
        var=pd.DataFrame(index=GENES),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]))
    a.layers['counts'] = csr_matrix(np.full(a.shape, 5, dtype=np.float32))
    if decontx:
        a.layers[DECONTX_LAYER] = csr_matrix(np.full(a.shape, 3, dtype=np.float32))
    a.write_h5ad(path)
    return path


def _paths(tmp_path, decontx=(True, True)):
    return [_study(tmp_path / s / "processed_data" / f"{s}.h5ad", d)
            for s, d in zip(("S1", "S2"), decontx)]


def test_master_is_built_from_the_original_counts_by_default(tmp_path):
    out = concat_studies(_paths(tmp_path), tmp_path / "master.h5ad")
    assert set(np.unique(ad.read_h5ad(out).X.toarray())) == {5}


def test_master_is_built_from_the_decontx_counts_when_asked(tmp_path):
    out = concat_studies(_paths(tmp_path), tmp_path / "master.h5ad",
                         counts_layer=DECONTX_LAYER)
    master = ad.read_h5ad(out)
    assert master.n_obs == 6
    assert set(np.unique(master.X.toarray())) == {3}


def test_a_study_without_the_layer_is_named(tmp_path):
    with pytest.raises(ValueError, match="S2"):
        concat_studies(_paths(tmp_path, decontx=(True, False)),
                       tmp_path / "master.h5ad", counts_layer=DECONTX_LAYER)


def test_the_record_says_which_counts_a_master_was_built_from(tmp_path):
    record_dir = tmp_path / "_master" / "processed_data"
    assert provenance.combined_counts_layer(record_dir) is None
    provenance.record_stage(record_dir, "master", "combine",
                            {"studies": ["S1"], "counts_layer": DECONTX_LAYER})
    assert provenance.combined_counts_layer(record_dir) == DECONTX_LAYER


def test_methods_call_a_decontx_masters_counts_decontaminated(tmp_path):
    """A DE on a DecontX master reads its 'raw' counts, which are the
    corrected ones; the Methods text must say so and cite DecontX."""
    from kosmic.gui.de_analysis.pages.methods_page import MethodsPage
    decontx_master = {"stages": [{"stage": "combine",
                                  "params": {"counts_layer": DECONTX_LAYER}}]}
    plain_master = {"stages": [{"stage": "combine", "params": {"counts_layer": None}}]}
    assert "DecontX" in MethodsPage._count_source_text(None, decontx_master)
    assert "DecontX" in MethodsPage._count_source_text("raw", decontx_master)
    assert MethodsPage._count_source_text("raw", plain_master) == "raw counts"
    assert MethodsPage._count_source_text("raw", None) == "raw counts"
    assert "DecontX" in MethodsPage._count_source_text(DECONTX_LAYER, None)


def test_combine_offers_decontx_only_when_every_study_has_it():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from types import SimpleNamespace
    from PyQt6.QtWidgets import QApplication, QCheckBox
    from kosmic.gui.combine.workspace import CombineDialog
    from kosmic.gui.shared.widgets import SecondaryLabel

    app = QApplication.instance() or QApplication([])  # noqa: F841 (held for the test)
    page = SimpleNamespace(
        _decontx_check=QCheckBox(), _decontx_label=SecondaryLabel(""),
        _study_obs={"S1": {"layers": ["counts", DECONTX_LAYER]},
                    "S2": {"layers": ["counts"]}},
        _selected_accessions=lambda: selected)
    refresh = CombineDialog._refresh_decontx_option

    selected = ["S1", "S2"]
    refresh(page)
    assert not page._decontx_check.isEnabled()
    assert "S2" in page._decontx_label.text()

    selected = ["S1"]
    refresh(page)
    assert page._decontx_check.isEnabled()
    page._decontx_check.setChecked(True)

    selected = ["S1", "S2"]          # adding a study without the layer unticks it
    refresh(page)
    assert not page._decontx_check.isChecked()
