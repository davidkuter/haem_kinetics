import numpy as np
import pandas as pd

from haem_kinetics.components.experimental_data import ExperimentalData
from haem_kinetics.components.fit_metrics import dv_total_fe_fg, score_fractionation
from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.model2 import Model2


def test_perfect_match_is_zero_rmse():
    exp = ExperimentalData()
    exp.no_drug_dd2()
    df = pd.DataFrame({
        'conc_hb_dv': exp.data['Hb'],
        'conc_fe3pp': exp.data['Hm'],
        'conc_hz': exp.data['Hz'],
        'conc_fe2pp': 0.0,
        'conc_hb_rbc': 0.0,
    }, index=exp.data.index)
    m = score_fractionation(df, exp)
    for name in ('Hb', 'Hm', 'Hz', 'DV_Fe'):
        assert m[name]['rmse'] < 1e-12
        assert m[name]['chi2_red'] < 1e-12


def test_model2_hz_beats_model1():
    kwargs = dict(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    m1 = Model1()
    m1.run(**kwargs)
    m1.score_vs_experiment()
    m2 = Model2()
    m2.run(**kwargs)
    m2.score_vs_experiment()
    assert m2.fit_metrics['Hz']['rmse'] < m1.fit_metrics['Hz']['rmse']
    assert m2.fit_metrics['DV_Fe']['rmse'] < m1.fit_metrics['DV_Fe']['rmse']
    assert np.isfinite(m1.fit_metrics['Hb']['chi2_red'])


def test_dv_fe_skips_hb_assay_and_free_haem_lumped_cols():
    """Assay Hb and aq+lip Hm lumps must not double-count in DV Fe."""
    idx = [20.0, 24.0]
    df = pd.DataFrame({
        'conc_hb_htv': [1.0, 1.0],
        'conc_hb_dv': [0.0, 0.0],
        'conc_hb_assay': [1.0, 1.0],
        'conc_fe2pp': [0.0, 0.0],
        'conc_fe3pp_aq': [0.5, 0.5],
        'conc_fe3pp_lip': [0.5, 0.5],
        'conc_fe3pp_free': [1.0, 1.0],
        'conc_hz': [0.0, 0.0],
        'conc_hb_rbc': [10.0, 10.0],
    }, index=idx)
    tot = dv_total_fe_fg(df)
    assert abs(float(tot.iloc[0]) - 2.0) < 1e-12


def test_dv_fe_skips_hb_assay_column():
    """Assay Hb is native+globin; must not double-count in DV Fe."""
    idx = [20.0, 24.0]
    df = pd.DataFrame({
        'conc_hb_dv': [1.0, 1.0],
        'conc_hb_globin': [1.0, 1.0],
        'conc_hb_assay': [2.0, 2.0],
        'conc_fe2pp': [0.0, 0.0],
        'conc_fe3pp': [0.0, 0.0],
        'conc_hz': [0.0, 0.0],
        'conc_hb_rbc': [10.0, 10.0],
    }, index=idx)
    tot = dv_total_fe_fg(df)
    assert abs(float(tot.iloc[0]) - 2.0) < 1e-12


def test_hb_assay_map_differs_from_native_only():
    exp = ExperimentalData()
    exp.no_drug_dd2()
    df = pd.DataFrame({
        'conc_hb_dv': 0.5,
        'conc_hb_globin': 0.5,
        'conc_fe3pp': exp.data['Hm'],
        'conc_hz': exp.data['Hz'],
        'conc_fe2pp': 0.0,
        'conc_hb_rbc': 0.0,
    }, index=exp.data.index)
    df['conc_hb_assay'] = df['conc_hb_dv'] + df['conc_hb_globin']
    m_sum = score_fractionation(
        df, exp, species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp', 'Hz': 'conc_hz'}
    )
    m_nat = score_fractionation(
        df, exp, species_map={'Hb': 'conc_hb_dv', 'Hm': 'conc_fe3pp', 'Hz': 'conc_hz'}
    )
    assert m_sum['Hb']['rmse'] != m_nat['Hb']['rmse']
    assert abs(m_sum['DV_Fe']['rmse'] - m_nat['DV_Fe']['rmse']) < 1e-12
