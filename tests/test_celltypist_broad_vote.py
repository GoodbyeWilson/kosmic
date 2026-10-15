"""Broad CellTypist labels and the majority vote over them.

On the DCM atlas a fine-label vote named three clusters after the type
with fewest fine labels while most cells were another type: a T/NK
cluster 'CD4+ T Cell' at 32% (CD8 40%), a macrophage subcluster 'CD16+
Monocyte' at 38% (LYVE1+ macrophage 46%). Lymphoid and myeloid subtypes
are grouped one level up, and with broad labels the vote is over them.
"""
import numpy as np
import pandas as pd
import pytest

from kosmic.scrna.annotate.score import collapse_celltypist_label as broad
from kosmic.scrna.annotate.score import majority_vote


@pytest.mark.parametrize('fine, expected', [
    ('CD4+T_naive', 'T Cell'), ('CD8+T_trans', 'T Cell'), ('MAIT-like', 'T Cell'),
    ('gdT', 'T Cell'), ('T/NK_cycling', 'T Cell'),
    ('NK_CD16hi', 'NK Cell'), ('ILC', 'NK Cell'),
    ('LYVE1+IGF1+MP', 'Macrophage'), ('LYVE1+MP_cycling', 'Macrophage'),
    ('MoMP', 'Macrophage'),
    ('CD14+Mo', 'Monocyte'), ('CD16+Mo', 'Monocyte'), ('DC', 'Dendritic Cell'),
    ('Mast', 'Mast Cell'), ('B_plasma', 'B Cell'),
    ('vCM3_stressed', 'Ventricular Cardiomyocyte'), ('EC7_endocardial', 'Endocardial Cell'),
])
def test_heart_model_labels_group_into_broad_types(fine, expected):
    assert broad(fine) == expected


def test_a_type_split_into_subtypes_wins_the_vote_by_its_cells():
    labels = (['CD16+Mo'] * 38 + ['LYVE1+IGF1+MP'] * 25 + ['LYVE1+TIMD4+MP'] * 21
              + ['DC'] * 16)
    clusters = ['0'] * 100
    assert set(majority_vote(clusters, labels)) == {'CD16+Mo'}          # fine vote
    voted = majority_vote(clusters, labels, broad)
    # The winning broad type, named by its most common fine label.
    assert set(voted) == {'LYVE1+IGF1+MP'} and broad('LYVE1+IGF1+MP') == 'Macrophage'


def test_each_cluster_votes_on_its_own_cells_and_unclustered_cells_get_none():
    clusters = pd.Series(['0'] * 4 + ['1'] * 4 + [np.nan], dtype=object)
    labels = ['vCM1', 'vCM1', 'vCM2', 'FB1', 'FB1', 'FB2', 'FB2', 'FB2', 'vCM1']
    voted = majority_vote(clusters, labels, broad)
    assert list(voted[:8]) == ['vCM1'] * 4 + ['FB2'] * 4
    assert pd.isna(voted.iloc[8])


def test_without_grouping_the_vote_is_unchanged():
    labels = ['a', 'a', 'b'] * 3
    assert set(majority_vote(['0'] * 9, labels)) == {'a'}
