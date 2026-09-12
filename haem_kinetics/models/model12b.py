import math

from haem_kinetics.models.model12a import Model12a


class Model12b(Model12a):
    """
    Model 12a Fe speciation with Myburgh (2023) Model 2 empirical uptake.

    Myburgh fitted total DV heme-Fe to an exponential T(t) = A·exp(B·t) with
    A = 13.1 fg and B = 8.3×10⁻⁴ min⁻¹ (t in minutes from invasion; his
    Table 4.5). Since the only source of DV Fe is Hb uptake, the uptake mole
    rate is dT/dt = A·B·exp(B·t) (fg Fe/min), converted to mol/min via the Fe
    molar mass and to M/min on the current lumen volume. Unlike our
    f_exp × host law this is monotonic and does not slow as host Hb depletes,
    which tests whether the residual late-phase dip in 12a is a host-depletion
    artifact.

    Boundary condition (ours, host-conserving): unlike Myburgh (who holds
    [Hb_RBC] constant, an infinite RBC reservoir), 12b keeps the conserving
    ladder's finite host — uptake is drawn from it and stops when it is
    exhausted. Because Myburgh's empirical exponential is calibrated to deliver
    more Fe than our 85 fg host holds, the host runs out near 43.5 h and assay
    Hb/Hm fall sharply after that. That late cliff is **not** a bug: it is the
    honest consequence of combining Myburgh's over-delivering uptake with mass
    conservation. Model 12c drops conservation (Myburgh's own constant-[Hb_RBC]
    boundary condition) and so has no cliff.
    """

    UPTAKE_A_FG = 13.1
    UPTAKE_B_PER_MIN = 8.3e-4
    MW_FE_G_PER_MOL = 55.845
    PARASITE_T0_MIN = 16.0 * 60.0

    def __init__(self, model_name: str = 'Model 12b'):
        super().__init__(model_name=model_name)

    def _uptake_dv(self, t):
        # Myburgh's empirical rate, but truncated when the finite host is spent
        # (domain of the physical uptake: no Hb left to internalise). Host
        # depletion is inherited from the conserving base (no override here).
        if self._nonneg(self.initial_values[self.HOST_KEY]) <= 0.0:
            return 0.0
        t_abs = t + self.PARASITE_T0_MIN
        mole_rate = (
            self.UPTAKE_A_FG * 1e-15 * self.UPTAKE_B_PER_MIN
            * math.exp(self.UPTAKE_B_PER_MIN * t_abs) / self.MW_FE_G_PER_MOL
        )
        return mole_rate / self._vol_dv(t)
