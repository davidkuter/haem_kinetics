from pathlib import Path
import math

import pandas as pd

from haem_kinetics.models.helpers import (
    garnie_dd2_vol_dv_fl,
    garnie_dd2_vol_dv_L,
    garnie_pm_amount_scale,
    variable_dv_volume_L,
    _GARNIE_PM_AGES_H,
    _GARNIE_PM_I_PERCENT,
    _GARNIE_PM_IV_PERCENT,
    _garnie_pm_scale_knots,
)
from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.model2 import Model2a
from haem_kinetics.models.model2b import Model2b
from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.model4a import Model4a
from haem_kinetics.models.model4b import Model4b
from haem_kinetics.models.model5 import Model5
from haem_kinetics.models.model6 import Model6
from haem_kinetics.models.model7 import Model7
from haem_kinetics.models.model8 import Model8
from haem_kinetics.models.model9a import Model9a
from haem_kinetics.models.model9b import Model9b
from haem_kinetics.models.model9c import Model9c
from haem_kinetics.models.model10 import Model10
from haem_kinetics.models.model99 import Model99


def test_garnie_dd2_vol_matches_reconstructed_means():
    """Gompertz+collapse reconstruction vs xlsx/paper landmarks."""
    assert abs(garnie_dd2_vol_dv_fl(20.0) - 0.675) < 0.01
    assert abs(garnie_dd2_vol_dv_fl(24.0) - 2.278) < 0.01
    assert abs(garnie_dd2_vol_dv_fl(32.0) - 3.557) < 0.02
    assert abs(garnie_dd2_vol_dv_fl(36.0) - 0.763) < 0.01
    assert abs(garnie_dd2_vol_dv_fl(46.0)) < 1e-12


def test_garnie_dd2_vol_positive_in_trophozoite_window():
    vol, _dvol = garnie_dd2_vol_dv_L(0.0)
    assert vol > 0.0
    vol_end, _dvol = garnie_dd2_vol_dv_L(1700.0)  # ~44.3 h
    assert vol_end > 0.0


def test_variable_dv_volume_L_matches_garnie_schedule():
    for t in (0.0, 240.0, 960.0, 1700.0):
        assert variable_dv_volume_L(t) == garnie_dd2_vol_dv_L(t)


def _assert_seed_fg_and_total_fe(model):
    df = model.concentrations
    hb_cols = [c for c in (
        'conc_hb_htv', 'conc_hb_dv', 'conc_hb_globin',
        'conc_fe2pp', 'conc_fe3pp', 'conc_fe3pp_aq', 'conc_fe3pp_lip',
        'conc_fe3pp_xtal', 'conc_hz',
    ) if c in df.columns]
    start = float(df.iloc[0][hb_cols].sum())
    assert abs(start - (0.018 + 0.36) * 55.85) < 0.05
    tot = df[[c for c in df.columns
              if c.startswith('conc_') and c not in (
                  'conc_hb_assay', 'conc_hb_dv_obs', 'conc_fe3pp_free',
              )]].sum(axis=1)
    budget = model.const.total_fe_fg_cell
    assert abs(float(tot.iloc[0]) - budget) < 0.5
    assert abs(float(tot.iloc[-1]) - budget) < 0.5
    for col in hb_cols + ['conc_hb_rbc']:
        assert float(df[col].min()) > -1e-6


def test_model1_preserves_seed_fg_and_total_fe():
    model = Model1()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)


def test_model2a_preserves_seed_fg_and_total_fe():
    model = Model2a()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)


def test_model2b_preserves_seed_fg_and_total_fe():
    model = Model2b()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)


def test_model2b_uptake_is_two_phase_remaining_host():
    from haem_kinetics.models.helpers import (
        F_EXP_BREAK_T_MIN,
        fraction_exp_growth_two_phase,
    )
    model = Model2b()
    model.initial_values[model.HOST_KEY] = 0.02
    t_early = F_EXP_BREAK_T_MIN - 60.0
    t_late = F_EXP_BREAK_T_MIN + 60.0
    for t in (t_early, t_late):
        expected = (
            fraction_exp_growth_two_phase(
                t, model.a_early, model.b_early, model.a_late, model.b_late,
                t_break_min=model.t_break_min,
            )
            * 0.02 * model.const.vol_rbc / model._vol_dv(t)
        )
        assert abs(model._uptake_dv(t) - expected) < 1e-18
    model.initial_values[model.HOST_KEY] = 0.04
    assert abs(
        model._uptake_dv(t_late)
        - 2.0 * expected
    ) < 1e-18
    model.initial_values[model.HOST_KEY] = 0.0
    assert model._uptake_dv(t_late) == 0.0
    # Late specific rate is larger than early at the same leftover host.
    model.initial_values[model.HOST_KEY] = 0.02
    assert model._uptake_dv(t_late) > model._uptake_dv(t_early)
    # f_exp is continuous at the Fig. 5B join (t just below vs at t_break).
    t_star = F_EXP_BREAK_T_MIN
    assert abs(model._uptake_dv(t_star - 1e-9) - model._uptake_dv(t_star)) < 1e-12


def test_model3_preserves_seed_fg_and_total_fe():
    model = Model3()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)


def test_model4a_preserves_seed_fg_and_total_fe():
    model = Model4a()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_hb_globin' in model.concentrations.columns
    assert float(model.concentrations['conc_hb_globin'].iloc[0]) < 0.05


def test_model4b_preserves_seed_fg_and_total_fe():
    model = Model4b()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)


def test_model5_preserves_seed_fg_and_total_fe():
    model = Model5()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_hb_htv' in model.concentrations.columns
    # Seed Hb is HTV cargo (~1 fg), not lumen native.
    assert abs(float(model.concentrations['conc_hb_htv'].iloc[0]) - 0.018 * 55.85) < 0.05
    assert float(model.concentrations['conc_hb_dv'].iloc[0]) < 0.05
    assert 'conc_hb_globin' not in model.concentrations.columns
    assay = model.concentrations['conc_hb_htv'] + model.concentrations['conc_hb_dv']
    assert 'conc_hb_assay' in model.concentrations.columns
    assert (model.concentrations['conc_hb_assay'] - assay).abs().max() < 1e-9


def test_model6_preserves_seed_fg_and_total_fe():
    model = Model6()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_hb_htv' in model.concentrations.columns
    assert 'conc_fe3pp' in model.concentrations.columns
    assert 'conc_fe3pp_aq' not in model.concentrations.columns
    assay = model.concentrations['conc_hb_htv'] + model.concentrations['conc_hb_dv']
    assert (model.concentrations['conc_hb_assay'] - assay).abs().max() < 1e-9


def test_model6_k_release_tracks_s_pm():
    model = Model6()
    k_plat = model.const.k_htv_release
    t20 = (20.0 - 16.0) * 60.0
    t40 = (40.0 - 16.0) * 60.0
    t44 = (44.0 - 16.0) * 60.0
    s20 = garnie_pm_amount_scale(t20)
    s40 = garnie_pm_amount_scale(t40)
    s44 = garnie_pm_amount_scale(t44)
    assert abs(model._k_release(t20) - k_plat * s20) < 1e-15
    assert abs(model._k_release(t40) - k_plat * s40) < 1e-15
    assert abs(model._k_release(t44) - k_plat * s44) < 1e-15
    assert s40 > 0.99
    assert abs(s44 - 1.0) < 0.01
    m5 = Model5()
    assert abs(m5.const.k_htv_release - k_plat) < 1e-15
    assert model._k_release(t20) < m5.const.k_htv_release


def test_model7_preserves_seed_fg_and_total_fe():
    model = Model7()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_fe3pp_aq' in model.concentrations.columns
    assert 'conc_fe3pp_lip' in model.concentrations.columns
    assert 'conc_fe3pp' not in model.concentrations.columns
    assay = model.concentrations['conc_hb_htv'] + model.concentrations['conc_hb_dv']
    assert (model.concentrations['conc_hb_assay'] - assay).abs().max() < 1e-9


def test_model7_scores_hm_as_aq_plus_lip():
    from haem_kinetics.components.fit_metrics import score_fractionation

    model = Model7()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    model.score_vs_experiment()
    df = model.concentrations
    lumped = df['conc_fe3pp_aq'] + df['conc_fe3pp_lip']
    alt = score_fractionation(
        df.drop(columns=['conc_fe3pp_aq', 'conc_fe3pp_lip']).assign(conc_fe3pp=lumped),
        model.exp_data,
        species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp', 'Hz': 'conc_hz'},
    )
    assert abs(model.fit_metrics['Hm']['rmse'] - alt['Hm']['rmse']) < 1e-12
    assert abs(model.fit_metrics['DV_Fe']['rmse'] - alt['DV_Fe']['rmse']) < 1e-12


def test_model99_whatif_overrides_and_conserves():
    import math
    from haem_kinetics.models.model99 import HTV_RELEASE_K_SCALE, K_XTAL_WHATIF

    model = Model99()
    assert model.const.k_htv_release < math.log(2) / 20.0
    assert abs(model.const.k_htv_release - (math.log(2) / 20.0) * HTV_RELEASE_K_SCALE) < 1e-12
    assert abs(model.const.K_xtal - K_XTAL_WHATIF) < 1e-15
    t20 = (20.0 - 16.0) * 60.0
    t44 = (44.0 - 16.0) * 60.0
    assert model._k_release(t20) < model._k_release(t44)
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)


def test_model9a_preserves_seed_fg_and_total_fe():
    model = Model9a()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_fe3pp_xtal' in model.concentrations.columns
    assay = model.concentrations['conc_hb_htv'] + model.concentrations['conc_hb_dv']
    assert (model.concentrations['conc_hb_assay'] - assay).abs().max() < 1e-9


def test_model9a_area_factor_is_one_at_seed():
    model = Model9a()
    t0 = 0.0
    y0 = model._prepare_y0([0.018, 0.0, 0.0, 0.36], t0=t0)
    model._set_initial_conc(init=y0)
    xtal = 0.01
    model.initial_values['conc_fe3pp_xtal'] = xtal
    assert abs(model._hz_area_factor(t0) - 1.0) < 1e-12
    assert abs(model._hz_rate(t0) - model.const.k_hz * xtal) < 1e-15


def test_model9b_preserves_seed_fg_and_total_fe():
    model = Model9b()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert model.AREA_EXPONENT == 0.5


def test_model9c_preserves_seed_fg_and_total_fe():
    model = Model9c()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert abs(model.AREA_EXPONENT - 1.0 / 3.0) < 1e-15


def test_model10_preserves_seed_fg_and_total_fe():
    model = Model10()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_fe3pp_xtal' in model.concentrations.columns


def test_model10_amount_species_include_all_fe():
    model = Model10()
    assert 'conc_hb_htv' in model.AMOUNT_SPECIES
    assert 'conc_fe2pp' in model.AMOUNT_SPECIES
    assert 'conc_fe3pp_aq' in model.AMOUNT_SPECIES
    assert 'conc_fe3pp_lip' in model.AMOUNT_SPECIES
    assert 'conc_fe3pp_xtal' in model.AMOUNT_SPECIES
    assert 'conc_hz' in model.AMOUNT_SPECIES
    # Only Hb_dv remains a lumen species
    assert 'conc_hb_dv' not in model.AMOUNT_SPECIES


def test_model8_preserves_seed_fg_and_total_fe():
    model = Model8()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_fe3pp_xtal' in model.concentrations.columns
    assert 'conc_fe3pp' not in model.concentrations.columns
    assay = model.concentrations['conc_hb_htv'] + model.concentrations['conc_hb_dv']
    assert (model.concentrations['conc_hb_assay'] - assay).abs().max() < 1e-9


def test_model8_scores_hm_as_aq_plus_lip_plus_xtal():
    from haem_kinetics.components.fit_metrics import score_fractionation

    model = Model8()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    model.score_vs_experiment()
    df = model.concentrations
    lumped = df['conc_fe3pp_aq'] + df['conc_fe3pp_lip'] + df['conc_fe3pp_xtal']
    alt = score_fractionation(
        df.drop(columns=['conc_fe3pp_aq', 'conc_fe3pp_lip', 'conc_fe3pp_xtal']).assign(
            conc_fe3pp=lumped
        ),
        model.exp_data,
        species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp', 'Hz': 'conc_hz'},
    )
    assert abs(model.fit_metrics['Hm']['rmse'] - alt['Hm']['rmse']) < 1e-12
    assert abs(model.fit_metrics['DV_Fe']['rmse'] - alt['DV_Fe']['rmse']) < 1e-12


def test_model8_fe3_seed_goes_to_aqueous():
    model = Model8()
    y0 = model._prepare_y0([0.018, 0.0, 0.01, 0.36], t0=0.0)
    v0 = model._vol_dv(0.0)
    scale = model.const.vol_dv / v0
    assert abs(y0[0] - 0.018) < 1e-15
    assert abs(y0[3] - 0.01 * scale) < 1e-15
    assert abs(y0[4]) < 1e-15
    assert abs(y0[5]) < 1e-15
    assert abs(y0[6] - 0.36 * scale) < 1e-15


def test_model7_fe3_seed_goes_to_aqueous():
    model = Model7()
    y0 = model._prepare_y0([0.018, 0.0, 0.01, 0.36], t0=0.0)
    v0 = model._vol_dv(0.0)
    scale = model.const.vol_dv / v0
    assert abs(y0[0] - 0.018) < 1e-15
    assert abs(y0[3] - 0.01 * scale) < 1e-15
    assert abs(y0[4]) < 1e-15
    assert abs(y0[5] - 0.36 * scale) < 1e-15


def test_htv_fg_uses_reference_volume_not_lumen():
    """AMOUNT_SPECIES fg is C·V_ref; lumen fg is C·V_DV(t)."""
    model = Model5()
    idx = pd.Index([16.0, 32.0])
    c_htv = 0.018
    c_lumen = 0.018
    df = pd.DataFrame({
        'conc_hb_htv': [c_htv, c_htv],
        'conc_hb_dv': [c_lumen, c_lumen],
        'conc_hb_rbc': [0.0, 0.0],
    }, index=idx)
    out = model._concentrations_to_fgcell(df)
    factor = (10 ** 15) * 55.85
    v_ref = model.const.vol_dv
    expected_htv = c_htv * v_ref * factor
    assert abs(float(out['conc_hb_htv'].iloc[0]) - expected_htv) < 1e-9
    assert abs(float(out['conc_hb_htv'].iloc[1]) - expected_htv) < 1e-9
    assert abs(float(out['conc_hb_dv'].iloc[1]) - float(out['conc_hb_dv'].iloc[0])) > 0.1


def test_htv_y0_not_scaled_by_lumen():
    model = Model5()
    y0 = model._prepare_y0([0.018, 0.0, 0.0, 0.36], t0=0.0)
    v0 = model._vol_dv(0.0)
    scale = model.const.vol_dv / v0
    assert abs(y0[0] - 0.018) < 1e-15
    assert abs(y0[4] - 0.36 * scale) < 1e-15


def test_pm_scale_holds_16_to_20h_and_plateau_is_one():
    s16 = garnie_pm_amount_scale(0.0)
    s20 = garnie_pm_amount_scale((20.0 - 16.0) * 60.0)
    assert abs(s16 - s20) < 1e-12
    knots = _garnie_pm_scale_knots()
    assert abs(0.5 * (knots[-2] + knots[-1]) - 1.0) < 1e-12
    # Lag is a large fraction of plateau (Fig. 3: PMs already present at 20 h)
    assert s20 > 0.3
    assert s20 < 0.55


def test_pm_scale_piecewise_linear_midpoint():
    knots = _garnie_pm_scale_knots()
    t28 = (28.0 - 16.0) * 60.0
    t32 = (32.0 - 16.0) * 60.0
    t_mid = 0.5 * (t28 + t32)
    s_mid = garnie_pm_amount_scale(t_mid)
    expected = 0.5 * (knots[_GARNIE_PM_AGES_H.index(28.0)] + knots[_GARNIE_PM_AGES_H.index(32.0)])
    assert abs(s_mid - expected) < 1e-12


def test_pm_percent_knots_match_figshare_xlsx():
    path = Path(__file__).resolve().parents[3] / 'data' / 'garnie' / 'PM_I_PM_IV_enzyme_analysis.xlsx'
    assert path.is_file()

    def _average_percent(sheet):
        df = pd.read_excel(path, sheet_name=sheet, header=None)
        start = None
        for i in range(len(df)):
            label = str(df.iloc[i, 0])
            if 'Average of Percent values' in label and 'BiP' not in label:
                start = i + 2  # skip the Hours header row
                break
        assert start is not None
        by_age = {}
        for i in range(start, start + len(_GARNIE_PM_AGES_H)):
            by_age[float(df.iloc[i, 0])] = float(df.iloc[i, 5])
        return tuple(by_age[h] for h in _GARNIE_PM_AGES_H)

    pmi = _average_percent('PM I')
    pmiv = _average_percent('PM IV')
    for a, b in zip(pmi, _GARNIE_PM_I_PERCENT):
        assert abs(a - b) < 1e-4
    for a, b in zip(pmiv, _GARNIE_PM_IV_PERCENT):
        assert abs(a - b) < 1e-4
