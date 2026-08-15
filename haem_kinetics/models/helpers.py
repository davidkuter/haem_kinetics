"""Shared kinetic helpers."""
import math
from typing import Tuple


# Garnie, Egan & Wicht, Commun. Biol. (2025) Fig. 2B — Dd2 aqueous DV lumen.
# https://doi.org/10.1038/s42003-025-08991-z
# Growth: Gompertz Y = A·exp(−exp(−k·(age−t_i))); A = 3.7 fL is the published
# Dd2 peak (plateau begins near 32 h). k and t_i are reconstructed by least
# squares to mean lumen volumes in data/garnie/phrodo_Dd2.xlsx at 20, 24, 32 h
# (the paper does not print fit coefficients).
# Collapse: linear to (46 h, 0) as published; slope from 36 h and 40 h means
# constrained through that point. 32–36 h is a linear join between the two
# published pieces (Garnie fitted growth and collapse separately).
_DD2_GOMPERTZ_A_FL = 3.7
_DD2_GOMPERTZ_K_PER_H = 0.31377421
_DD2_GOMPERTZ_TI_H = 21.69318223
_DD2_GROWTH_END_H = 32.0
_DD2_COLLAPSE_START_H = 36.0
_DD2_COLLAPSE_END_H = 46.0
_DD2_COLLAPSE_SLOPE_FL_PER_H = 0.07634190716418995


def _dd2_gompertz_fl(age_h: float) -> float:
    return _DD2_GOMPERTZ_A_FL * math.exp(
        -math.exp(-_DD2_GOMPERTZ_K_PER_H * (age_h - _DD2_GOMPERTZ_TI_H))
    )


def _dd2_gompertz_dfl_per_h(age_h: float) -> float:
    y = _dd2_gompertz_fl(age_h)
    return y * _DD2_GOMPERTZ_K_PER_H * math.exp(
        -_DD2_GOMPERTZ_K_PER_H * (age_h - _DD2_GOMPERTZ_TI_H)
    )


def garnie_dd2_vol_dv_fl(age_h: float) -> float:
    """Dd2 DV lumen volume (fL) vs hours post-invasion (Garnie Fig. 2B form)."""
    if age_h <= _DD2_GROWTH_END_H:
        return _dd2_gompertz_fl(age_h)
    v_growth_end = _dd2_gompertz_fl(_DD2_GROWTH_END_H)
    v_collapse_start = _DD2_COLLAPSE_SLOPE_FL_PER_H * (
        _DD2_COLLAPSE_END_H - _DD2_COLLAPSE_START_H
    )
    if age_h < _DD2_COLLAPSE_START_H:
        span = _DD2_COLLAPSE_START_H - _DD2_GROWTH_END_H
        frac = (age_h - _DD2_GROWTH_END_H) / span
        return v_growth_end + frac * (v_collapse_start - v_growth_end)
    return _DD2_COLLAPSE_SLOPE_FL_PER_H * (_DD2_COLLAPSE_END_H - age_h)


def garnie_dd2_dvol_dv_fl_per_h(age_h: float) -> float:
    """dV/d(age_h) in fL/h for garnie_dd2_vol_dv_fl."""
    if age_h <= _DD2_GROWTH_END_H:
        return _dd2_gompertz_dfl_per_h(age_h)
    v_growth_end = _dd2_gompertz_fl(_DD2_GROWTH_END_H)
    v_collapse_start = _DD2_COLLAPSE_SLOPE_FL_PER_H * (
        _DD2_COLLAPSE_END_H - _DD2_COLLAPSE_START_H
    )
    if age_h < _DD2_COLLAPSE_START_H:
        span = _DD2_COLLAPSE_START_H - _DD2_GROWTH_END_H
        return (v_collapse_start - v_growth_end) / span
    return -_DD2_COLLAPSE_SLOPE_FL_PER_H


def garnie_dd2_vol_dv_L(
    t_min: float,
    parasite_t0_h: float = 16.0,
) -> Tuple[float, float]:
    """
    Dd2 DV lumen volume (L) and dV/dt (L/min) at simulation time t_min.

    Simulation t is minutes from parasite_t0_h. Volume is the Garnie aqueous
    lumen (Hz excluded). Not an uptake law.
    """
    age_h = parasite_t0_h + t_min / 60.0
    vol_fl = garnie_dd2_vol_dv_fl(age_h)
    dvol_fl_per_h = garnie_dd2_dvol_dv_fl_per_h(age_h)
    vol_L = vol_fl * 1e-15
    dvol_L_per_min = dvol_fl_per_h * 1e-15 / 60.0
    return vol_L, dvol_L_per_min


def fraction_exp_growth(t: float, a: float = 0.1578, b: float = 0.001102) -> float:
    """Fractional exponential growth for Hb uptake (t in minutes from 16 h).

    Phenomenological accelerating rate of remaining host Hb — not V_DV(t).
    Empirically fit to cumulative DV Fe (see docs/models/model2.md).
    """
    return a * b * math.exp(b * t)


def enzyme_logistic_scale(
    t_min: float,
    parasite_t0_h: float = 16.0,
    t_mid_h: float = 26.0,
    steepness: float = 0.35,
) -> float:
    """
    Mature protease capacity vs parasite age (0→1), independent of uptake f_exp.

    Logistic in hours post-invasion; shaped after DV/protease build-up during
    the trophozoite window (Garnie-like), not matched to the uptake prefactor.
    """
    age_h = parasite_t0_h + t_min / 60.0
    return 1.0 / (1.0 + math.exp(-steepness * (age_h - t_mid_h)))


def lipid_aqueous_fraction(vol_fract_lip: float, k_partition: float) -> float:
    """
    Equilibrium aqueous fraction of total Fe(III) for DV-referenced amounts.

    Same formula as Constants.compute_lipid_seq_constant().
    """
    return (1.0 - vol_fract_lip) / (1.0 + vol_fract_lip + vol_fract_lip * k_partition)


def lipid_over_aq_ratio(vol_fract_lip: float, k_partition: float) -> float:
    """Equilibrium ratio [Fe3]_lip / [Fe3]_aq for DV-referenced concentrations."""
    phi = lipid_aqueous_fraction(vol_fract_lip, k_partition)
    if phi <= 0.0 or phi >= 1.0:
        raise ValueError(f'Invalid aqueous fraction phi={phi}')
    return (1.0 - phi) / phi
