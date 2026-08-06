"""Shared kinetic helpers for Models 5+."""
import math


def fraction_exp_growth(t: float, a: float = 0.1578, b: float = 0.001102) -> float:
    """Fractional exponential growth used in Model 3/5 (t in minutes from 16 h)."""
    return a * b * math.exp(b * t)


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


def vol_dv_fl(t_min: float, parasite_t0_h: float = 16.0) -> float:
    """
    Approximate Dd2 DV lumen volume (fL) vs simulation time (minutes from t0).

    Shape inspired by Combrink et al. 2025 (Gompertz growth then late collapse).
    Peak ~3.7 fL near ~32 h parasite age. Parameters are approximate placeholders
    for mechanistic modeling, not a digitised fit of the published curve.
    """
    age_h = parasite_t0_h + t_min / 60.0
    # Gompertz-like rise: V = Vmax * exp(-a * exp(-b * age))
    vmax = 3.7
    a = 4.0
    b = 0.18
    v_grow = vmax * math.exp(-a * math.exp(-b * age_h))
    # Linear collapse after ~36 h toward a small residual lumen
    if age_h <= 36.0:
        return max(v_grow, 0.05)
    # Collapse from V(36) toward ~0.3 fL by 46 h
    v36 = vmax * math.exp(-a * math.exp(-b * 36.0))
    frac = min(max((age_h - 36.0) / 10.0, 0.0), 1.0)
    return max(v36 * (1.0 - frac) + 0.3 * frac, 0.05)


def fg_to_molar(fg_per_cell: float, vol_fl: float) -> float:
    """Convert fg Fe/cell to M using current DV volume in fL."""
    # M = fg / (V_fL * MW_Fe); because fg = M * V_fL * 55.85
    return fg_per_cell / (vol_fl * 55.85)


def molar_to_fg(conc_m: float, vol_fl: float) -> float:
    """Convert M (DV basis) to fg Fe/cell."""
    return conc_m * vol_fl * 55.85


def sigmoidal_uptake_rate_per_min(
    t_min: float,
    remaining_host_fg: float,
    parasite_t0_h: float = 16.0,
    k_max: float = 0.0045,
    t_mid_h: float = 30.0,
    steepness: float = 0.35,
) -> float:
    """
    Rate of Fe uptake from host Hb into the DV (fg/cell/min).

    Uses a logistic schedule of a first-order drain on remaining host Fe:
        rate = k(t) * remaining_host_fg
    with k(t) rising through the trophozoite window (Combrink-like sigmoidal delivery).
    """
    if remaining_host_fg <= 0.0:
        return 0.0
    age_h = parasite_t0_h + t_min / 60.0
    # Suppress uptake during late lumen collapse
    if age_h >= 40.0:
        collapse = max(1.0 - (age_h - 40.0) / 6.0, 0.0)
    else:
        collapse = 1.0
    k_t = k_max / (1.0 + math.exp(-steepness * (age_h - t_mid_h)))
    return k_t * collapse * remaining_host_fg
