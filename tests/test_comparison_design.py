"""The merged pseudobulk design matched to the per-study designs (ADR-017 §4)."""
from pathlib import Path

import pandas as pd

from kosmic.comparison.design import Design, StudyDesign, differences, resolve_design


def _study(name, covariates=('sex_inferred',), sample='sample', **kw):
    return StudyDesign(study=name, de_path=Path(f'{name}_DE_deseq2.csv'),
                       covariates=tuple(covariates), sample_column=sample,
                       counts=kw.get('counts', 'decontX_counts'),
                       cell_type='Ventricular Cardiomyocyte',
                       cell_type_column='cell_type_atlas', min_cells=10, min_counts=0)


def _obs(donors, sexes, extra=None):
    """Two cells per donor; 'sample' is the donor, 'sex' and 'sex_inferred' agree."""
    rows = []
    for d, s in zip(donors, sexes):
        rows += [{'sample': d, 'sex': s, 'sex_inferred': s, **(extra or {})}] * 2
    return pd.DataFrame(rows)


def _atlas(studies):
    """An atlas obs with each study's donors split across both conditions and sexes."""
    rows = []
    for st in studies:
        for i, (role, sex) in enumerate([('control', 'female'), ('control', 'male'),
                                         ('disease', 'female'), ('disease', 'male')]):
            rows.append({'study': st, 'sample': f'{st}_d{i}', '_role': role, 'sex': sex,
                         'cell_type': 'Ventricular Cardiomyocyte'})
    return pd.DataFrame(rows)


SEXES = ['female', 'male', 'female', 'male']


def _studies_obs(names):
    return {n: _obs([f'{n}_d{i}' for i in range(4)], SEXES) for n in names}


def test_inferred_sex_maps_to_the_atlas_sex_column():
    names = ['A', 'B']
    d = resolve_design([_study(n) for n in names], _studies_obs(names),
                       {'A': 'sex_inferred', 'B': 'sex'}, _atlas(names))
    assert d.mapping == {'sex_inferred': 'sex'}
    assert d.formula == '~sex + study + condition'
    assert d.warnings == [] and d.can_fit
    assert d.cell_type_column == 'cell_type'   # the atlas's own labels


def test_the_formula_holds_every_study_covariate_then_study():
    names = ['A', 'B']
    atlas = _atlas(names).assign(age=[30, 40, 50, 60] * 2)
    obs = {n: o.assign(age=1) for n, o in _studies_obs(names).items()}
    d = resolve_design([_study(n, ('sex_inferred', 'age')) for n in names], obs,
                       {n: 'sex' for n in names}, atlas)
    assert d.formula == '~sex + age + study + condition'


def test_a_designated_sex_column_that_disagrees_is_reported():
    names = ['A', 'B']
    obs = _studies_obs(names)
    obs['A'].loc[obs['A']['sample'] == 'A_d0', 'sex'] = 'male'   # recorded differs from inferred
    d = resolve_design([_study(n) for n in names], obs, {'A': 'sex', 'B': 'sex'}, _atlas(names))
    assert any("A: the model used 'sex_inferred'" in w and '1 donor(s) (A_d0)' in w
               for w in d.warnings)


def test_a_different_sample_column_is_accepted_when_it_names_the_same_donors():
    names = ['A', 'B']
    obs = _studies_obs(names)
    obs['A']['donor_id'] = obs['A']['sample'].str.replace('_d', '_donor')
    d = resolve_design([_study('A', sample='donor_id'), _study('B')], obs,
                       {n: 'sex' for n in names}, _atlas(names))
    assert d.warnings == []


def test_a_sample_column_naming_other_donors_is_reported():
    names = ['A', 'B']
    obs = _studies_obs(names)
    obs['A']['donor_id'] = 'one_donor'   # every cell one donor: 'sample' splits it further
    d = resolve_design([_study('A', sample='donor_id'), _study('B')], obs,
                       {n: 'sex' for n in names}, _atlas(names))
    assert any("A: 'donor_id' and 'sample' identify different donors" in w for w in d.warnings)


def test_studies_with_different_covariates_are_reported():
    names = ['A', 'B']
    d = resolve_design([_study('A'), _study('B', covariates=())], _studies_obs(names),
                       {n: 'sex' for n in names}, _atlas(names))
    assert any('different covariates' in w and 'B: none' in w for w in d.warnings)


def test_a_covariate_missing_from_the_atlas_blocks_fitting():
    names = ['A', 'B']
    d = resolve_design([_study(n, ('sex_inferred', 'batch')) for n in names],
                       _studies_obs(names), {n: 'sex' for n in names}, _atlas(names))
    assert d.unresolved == ['batch'] and not d.can_fit
    assert any("no column for the per-study covariate 'batch'" in w for w in d.warnings)


def test_a_covariate_that_never_varies_within_a_study_is_confounded():
    names = ['A', 'B']
    atlas = _atlas(names)
    atlas['site'] = atlas['study'].map({'A': 'x', 'B': 'y'})
    obs = {n: o.assign(site='x') for n, o in _studies_obs(names).items()}
    d = resolve_design([_study(n, ('sex_inferred', 'site')) for n in names], obs,
                       {n: 'sex' for n in names}, atlas)
    assert any("'site' does not vary within any study" in w for w in d.warnings)
    assert 'site' in d.covariates   # warned about, never dropped


def test_a_covariate_confounded_with_condition_is_reported():
    names = ['A', 'B']
    atlas = _atlas(names)
    atlas['sex'] = atlas['_role'].map({'control': 'female', 'disease': 'male'})
    d = resolve_design([_study(n) for n in names], _studies_obs(names),
                       {n: 'sex' for n in names}, atlas)
    assert any("'sex' is confounded with condition" in w for w in d.warnings)


def _design():
    return Design(covariates=['sex'], mapping={'sex_inferred': 'sex'}, sample_column='sample',
                  cell_type_column='cell_type', cell_type='Ventricular Cardiomyocyte',
                  counts='decontX_counts', min_cells=10, min_counts=0)


def test_a_saved_run_with_sex_and_study_matches():
    record = {'covariates': ['sex', 'study'], 'sample_column': 'sample',
              'min_cells': 10, 'min_counts': 0}
    assert differences(record, _design(), counts='decontX_counts') == []
    # Order of the terms does not matter.
    record['covariates'] = ['study', 'sex']
    assert differences(record, _design(), counts='decontX_counts') == []


def test_a_saved_run_without_study_does_not_match():
    out = differences({'covariates': ['sex']}, _design(), counts='decontX_counts')
    assert out == ["The saved merged run does not include study as a covariate."]


def test_a_saved_run_on_other_counts_or_covariates_does_not_match():
    out = differences({'covariates': ['study', 'age']}, _design(), counts='original_counts')
    assert any('covariates the studies did not use: age' in o for o in out)
    assert any('does not adjust for: sex' in o for o in out)
    assert any('used original_counts' in o for o in out)


def test_a_saved_run_without_a_record_does_not_match():
    assert differences(None, _design()) != []
