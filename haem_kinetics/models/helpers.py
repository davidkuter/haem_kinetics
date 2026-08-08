"""Shared kinetic helpers for Model 4+."""
import math


def fraction_exp_growth(t: float, a: float = 0.1578, b: float = 0.001102) -> float:
    """Fractional exponential growth used for Hb uptake (t in minutes from 16 h)."""
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
    the trophozoite window (Combrink-like), not matched to the uptake prefactor.
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
