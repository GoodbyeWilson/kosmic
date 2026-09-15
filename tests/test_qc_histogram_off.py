"""A QC bound of 0 means "off" everywhere, including the histogram preview.

The filter (kosmic.scrna.qc.filter) skips max_genes / max_counts when
they are 0, but the histogram's live pass/removed count treated 0 as a
real cap and reported every cell removed. The preview must agree with
what the filter would do.
"""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt6")
pytest.importorskip("pyqtgraph")


@pytest.fixture(scope="module")
def app():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _widget(app):
    from kosmic.gui.scrna.tabs.qc_tab import _QCHistogramWidget
    w = _QCHistogramWidget()
    rng = np.random.default_rng(0)
    obs = pd.DataFrame({
        'n_genes_by_counts': rng.integers(300, 5000, 500),
        'total_counts': rng.integers(1000, 20000, 500),
        'pct_counts_mt': rng.uniform(0, 3, 500),
    })
    w.set_data(obs)
    return w


def test_zero_max_is_no_bound(app):
    w = _widget(app)
    w.set_threshold("n_genes_by_counts", "min", 200.0)
    w.set_threshold("n_genes_by_counts", "max", 0.0)
    w.set_threshold("total_counts", "max", 0.0)
    w.set_threshold("pct_counts_mt", "max", 10.0)
    assert 'max' not in w._thresholds.get("n_genes_by_counts", {})
    assert 'max' not in w._thresholds.get("total_counts", {})
    assert "500 / 500 pass" in w._info_label.text()


def test_positive_max_after_zero_counts_cells(app):
    w = _widget(app)
    w.set_threshold("n_genes_by_counts", "max", 0.0)
    w.set_threshold("n_genes_by_counts", "max", 1000.0)
    assert w._thresholds["n_genes_by_counts"]["max"] == 1000.0
    assert "500 / 500 pass" not in w._info_label.text()
    # and back to off
    w.set_threshold("n_genes_by_counts", "max", 0.0)
    assert "500 / 500 pass" in w._info_label.text()


def test_excluded_count_is_stated_in_the_info_line(app):
    w = _widget(app)
    # Rebuild with an excluded count, as the QC tab does for 'exclude' roles.
    w.set_data(w._obs_df, n_excluded=32_138)
    text = w._info_label.text()
    assert "500 / 500 pass" in text
    assert "32,138 excluded cells not shown" in text
