from pathlib import Path
import math

import pandas as pd

from haem_kinetics.components.constants import Constants
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
from haem_kinetics.models.model12a import Model12a
from haem_kinetics.models.model12b import Model12b
from haem_kinetics.models.model12c import Model12c
from haem_kinetics.models.model13 import Model13
from haem_kinetics.models.model14a import Model14a
from haem_kinetics.models.model14b import Model14b
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


def test_model12a_preserves_seed_fg_and_total_fe():
    model = Model12a()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_seed_fg_and_total_fe(model)
    assert 'conc_fe3pp_aq' in model.concentrations.columns
    assert 'conc_fe3pp_lip' in model.concentrations.columns


def test_model12a_hz_forms_from_aqueous_not_lipid():
    """Myburgh Model 3: crystallisation is first-order in aqueous hematin only."""
    model = Model12a()
    model.initial_values['conc_fe3pp_aq'] = 0.0
    model.initial_values['conc_fe3pp_lip'] = 1e-3
    assert model._hz_rate() == 0.0  # lipid pool is a buffer, not the substrate
    model.initial_values['conc_fe3pp_aq'] = 2e-6
    assert abs(model._hz_rate() - model.const.k_hz * 2e-6) < 1e-18
    # The lipid derivative carries no crystallisation sink (exchange only).
    model.initial_values['conc_fe3pp_lip'] = 0.0
    assert abs(model._d_fe3pp_lip(0.0) - model._exchange_rate()) < 1e-18


def test_model12b_uptake_is_myburgh_exponential():
    model = Model12b()
    t = 600.0
    t_abs = t + model.PARASITE_T0_MIN
    mole_rate = (
        model.UPTAKE_A_FG * 1e-15 * model.UPTAKE_B_PER_MIN
        * math.exp(model.UPTAKE_B_PER_MIN * t_abs) / model.MW_FE_G_PER_MOL
    )
    expected = mole_rate / model._vol_dv(t)
    # Rate magnitude is independent of the remaining host (not f_exp × host)...
    for host in (0.02, 0.001):
        model.initial_values[model.HOST_KEY] = host
        assert abs(model._uptake_dv(t) - expected) < 1e-24
    # ...but truncates when the finite host is spent (conserving; domain of the
    # physical uptake — no Hb left to internalise).
    model.initial_values[model.HOST_KEY] = 0.0
    assert model._uptake_dv(t) == 0.0


def _assert_total_fe_conserved(model):
    """Host+DV Fe is conserved even when the uptake drains the host."""
    df = model.concentrations
    tot = df[[c for c in df.columns
              if c.startswith('conc_') and c not in (
                  'conc_hb_assay', 'conc_hb_dv_obs', 'conc_fe3pp_free',
              )]].sum(axis=1)
    budget = model.const.total_fe_fg_cell
    assert abs(float(tot.iloc[0]) - budget) < 0.5
    assert abs(float(tot.iloc[-1]) - budget) < 0.5


def _assert_no_late_cliff(model):
    """Myburgh's constant [Hb_RBC]: monotonic late phase, no finite-host cliff."""
    df = model.concentrations.copy()
    df['hm'] = df['conc_fe3pp_aq'] + df['conc_fe3pp_lip'] + df.get('conc_fe2pp', 0.0)
    df['hb'] = df['conc_hb_htv'] + df['conc_hb_dv']
    times = df.index.to_numpy(dtype=float)
    late = df[times >= 40.0]
    # No sharp drop in the last 4 h: 44 h stays within 5% of the 40-44 h peak.
    for col in ('hb', 'hm', 'conc_hz'):
        peak = float(late[col].max())
        end = float(late[col].iloc[-1])
        assert end > 0.95 * peak, f'{col}: late cliff {end} vs peak {peak}'


def test_model12b_conserves_total_fe_and_exhausts_host():
    """12b keeps the conserving ladder: Myburgh's over-delivering uptake drains
    the finite host (~0 by 44 h), the honest late cliff. Contrast 12c."""
    model = Model12b()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_total_fe_conserved(model)
    host = model.concentrations['conc_hb_rbc']
    assert float(host.min()) > -0.01
    assert float(host.iloc[-1]) < 0.5


def test_model12c_constant_volume_and_enzyme():
    model = Model12c()
    assert model.variable_dv_volume is False
    assert model._vol_dv(0.0) == model.const.vol_dv
    assert model._vol_dv(1700.0) == model.const.vol_dv
    # Constant [E]: no s_PM schedule, so identical early vs late.
    e_early = model._enzyme_conc('plm_2', 100.0)
    e_late = model._enzyme_conc('plm_2', 1600.0)
    assert abs(e_early - e_late) < 1e-30
    assert abs(e_early - model.const.conc_enzymes['plm_2']) < 1e-30


def test_model12c_constant_host_no_cliff():
    model = Model12c()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    host = model.concentrations['conc_hb_rbc']
    assert abs(float(host.max()) - float(host.min())) < 1e-9
    _assert_no_late_cliff(model)


def test_model13_uses_upper_range_host_budget():
    """Model 13 = Model 12b uptake on an upper-range MCHC Fe budget.

    The 106 fg default is built from population-mean MCHC (34 g/dL); Model 13
    uses the clinical upper bound (36 g/dL), a larger but cited host pool.
    """
    model = Model13()
    assert model.MCHC_G_PER_DL == 36.0
    # Larger host pool than the mean-cell budget, but still same uptake law.
    assert model._full_hb_rbc_m() > Constants.compute_conc_hb_rcb()
    assert model.const.total_fe_fg_cell > 110.0
    # Uptake is inherited from 12b (Myburgh exponential): host-independent
    # magnitude, truncating only when the host is spent.
    model.initial_values[model.HOST_KEY] = 0.02
    t = 600.0
    rate_full = model._uptake_dv(t)
    model.initial_values[model.HOST_KEY] = 0.001
    assert abs(model._uptake_dv(t) - rate_full) < 1e-30
    model.initial_values[model.HOST_KEY] = 0.0
    assert model._uptake_dv(t) == 0.0


def test_model13_no_cliff_host_not_exhausted():
    """Upper-range budget keeps the host from emptying, so Myburgh's uptake
    never stops dead and the standing Hb/Hm pools do not collapse at 44 h."""
    model = Model13()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_total_fe_conserved(model)  # conserved at the model's own ~112 fg budget
    host = model.concentrations['conc_hb_rbc']
    # Contrast 12b (host → 0, cliff): here the host is drawn down but survives.
    assert float(host.min()) > 0.3
    _assert_no_late_cliff(model)


def test_model14a_release_decoupled_from_spm():
    """14a: inner-vesicle lysis is the constant Klemba rate, not gated by
    plasmepsin amount — so release no longer collapses early and the 24-26 h
    assay-Hb hump is removed."""
    model = Model14a()
    early = model._k_release(100.0)
    late = model._k_release(1600.0)
    assert early == late == model.const.k_htv_release
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_total_fe_conserved(model)
    _assert_no_late_cliff(model)
    # No spurious early hump: 26 h assay Hb no longer exceeds the 44 h value
    # (contrast Model 13, where low early s_PM piled cargo up to ~2.9 fg).
    df = model.concentrations
    hb26 = float(df.iloc[(df.index - 26.0).to_series().abs().values.argmin()]['conc_hb_assay'])
    hb44 = float(df.iloc[(df.index - 44.0).to_series().abs().values.argmin()]['conc_hb_assay'])
    assert hb26 < hb44


def test_model14b_slower_constant_release_than_14a():
    """14b: same decoupling as 14a but a slower provisional lysis t½ (30 min),
    a larger standing pool that fits early but overshoots late."""
    model = Model14b()
    assert model.HTV_LYSIS_T_HALF_MIN == 30.0
    r = model._k_release(500.0)
    assert r == model._k_release(1500.0)  # constant, decoupled from s_PM
    assert abs(r - math.log(2.0) / 30.0) < 1e-12
    # Slower than 14a's Klemba rate (t½ 20 min) → larger standing HTV pool.
    assert r < Model14a()._k_release(500.0)
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    _assert_total_fe_conserved(model)


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
