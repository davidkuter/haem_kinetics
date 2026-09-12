from haem_kinetics.models.helpers import elliott_cytostome_maturity
from haem_kinetics.models.model13 import Model13


class Model15b(Model13):
    """
    Model 13 with Elliott (2008) cytostomal-uptake maturity gating Myburgh v_up.

    The single change from Model 13 is the uptake *amount* before mid-trophozoite.
    Myburgh's exponential is a continuous cytostome-like feed from invasion.
    Elliott et al. (PNAS 2008) show that is the wrong early process: rings take
    up Hb mainly as a one-shot Big Gulp, and the cytostome's small-vesicle
    pathway only "increases its contribution to total hemoglobin uptake" once
    the parasite is a trophozoite (24–30 h post-invasion). Model 15b therefore
    multiplies Model 13's Myburgh rate by a Hermite smoothstep that is 0 before
    24 h and 1 after 30 h — Elliott's window, not Garnie Fig. 5B.

    Release stays k_release ∝ s_PM(t) (the Model 13 / 6 clock). Host
    conservation and the upper-range MCHC budget are unchanged. No free
    amplitude: the gate is 0 or 1 outside the cited window.

    Consequence: 16–24 h has no continuous Myburgh feed (ring Big Gulp is not
    that law). Early standing HTV should drop; whether the 29 h dip survives
    depends on how much cargo is left when the gate opens.
    """

    def __init__(self, model_name: str = 'Model 15b'):
        super().__init__(model_name=model_name)

    def _uptake_maturity(self, t):
        return elliott_cytostome_maturity(t)

    def _uptake_dv(self, t):
        return self._uptake_maturity(t) * super()._uptake_dv(t)
