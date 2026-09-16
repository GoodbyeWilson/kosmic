"""The elbow of the PCA variance curve (kosmic.visualisation.scrna.pca_plots).

Single-cell data never reaches 80% of total variance within the PCs
computed, so the old "PC at 80% cumulative variance" estimate always
returned the last PC and won the max -- "Elbow detected at PC 80" on a
run of 80 PCs. The elbow is now the knee alone, floored at min_pcs, and
a knee on the first or last PC is reported as not found.
"""
import numpy as np

from kosmic.visualisation.scrna.pca_plots import elbow_report, find_elbow


def _curve(n=80, knee=18):
    """Steep fall to the knee, then a slow tail -- a typical scRNA scree."""
    pcs = np.arange(n)
    y = np.where(pcs < knee, 0.06 * np.exp(-pcs / 6), 0.0025 * np.exp(-(pcs - knee) / 200))
    return y / y.sum() * 0.6          # 60% of total variance in these PCs


def test_knee_is_found_and_not_the_ceiling():
    r = elbow_report(_curve())
    assert r['found']
    assert 10 <= r['suggested'] <= 30
    assert r['suggested'] != 80
    assert find_elbow(_curve()) == r['suggested']


def test_floor_applies():
    y = np.r_[0.5, 0.3, 0.1, np.full(37, 0.002)]      # knee at PC 3-4
    r = elbow_report(y, min_pcs=10)
    assert r['elbow'] <= 4 and r['found'] and r['suggested'] == 10


def test_no_knee_is_reported_not_invented():
    y = np.exp(np.linspace(np.log(0.02), np.log(0.001), 40))   # straight on log scale: no elbow
    r = elbow_report(y)
    assert not r['found']
    assert r['suggested'] == 10
