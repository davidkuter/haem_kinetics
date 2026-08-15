from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.native_hb import native_tetramer_rate


class Model4b(Model3):
    """
    Model 3 + native-Hb-competent enzymes only on a single DV Hb pool.

    v_dig uses PM I / PM II / FP-2 with the shared native rate law (same as 4a).
    HAP, PM IV, and FP-3 do not act on the native tetramer. Nick and haem
    release are lumped.
    """

    def __init__(self, model_name: str = 'Model 4b'):
        super().__init__(model_name=model_name)

    def _hb_removal(self, t):
        hb = self._nonneg(self.initial_values['conc_hb_dv'])
        return native_tetramer_rate(self, hb, t)
