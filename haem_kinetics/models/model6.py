from haem_kinetics.models.model5 import Model5
from haem_kinetics.models.helpers import garnie_pm_amount_scale


class Model6(Model5):
    """
    Model 5 + inner-vesicle lysis clocked by Garnie Fig. 3 s_PM(t).

    k_release(t) = k_htv_release · s_PM(t). Plateau (40–44 h) keeps
    Klemba t½ = 20 min; early s_PM ≈ 0.41 so early t½ is longer.
    Same states and scoring as Model 5. Not a Garnie Hb fit.
    """

    def __init__(self, model_name: str = 'Model 6'):
        super().__init__(model_name=model_name)

    def _k_release(self, t):
        return self.const.k_htv_release * garnie_pm_amount_scale(t)

    def _htv_release_lumen(self, t):
        """Cargo appearance in the lumen (M/min on current V_DV). Inner-vesicle lysis."""
        htv = self._nonneg(self.initial_values['conc_hb_htv'])
        if htv <= 0.0:
            return 0.0
        return self._k_release(t) * htv * self.const.vol_dv / self._vol_dv(t)

    def _d_hb_htv(self, t):
        htv = self._nonneg(self.initial_values['conc_hb_htv'])
        appear = self._uptake_dv(t) * self._vol_dv(t) / self.const.vol_dv
        return appear - self._k_release(t) * htv
