from haem_kinetics.models.model1 import Model1
from haem_kinetics.models.helpers import fraction_exp_growth


class Model2(Model1):
    """
    Model 1 + accelerating host→DV uptake via f_exp(t).

    Protease amount stays at PaxDB n_E; molarity is n_E / V_DV(t) (shared
    bookkeeping). Single Fe(III) pool; no lipid φ. f_exp is not V_DV(t).
    """

    def __init__(self, model_name: str = 'Model 2'):
        super().__init__(model_name=model_name)

    def _uptake_dv(self, t):
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        tot_hb_conc = host * self.const.vol_rbc / self._vol_dv(t)
        return fraction_exp_growth(t) * tot_hb_conc
