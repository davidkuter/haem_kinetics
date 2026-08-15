"""Model 9a: Model 8 + crystal-area growth of haemozoin (sphere geometry).

v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{2/3}

n_Hz is amount (C · V_DV(t)), not lumen concentration, so dilution does not
fake a change in crystal size. n_Hz_start is the protocol seed: 0.36 M at
V_ref (~20 fg). The 2/3 exponent is sphere geometry, not a Garnie fit.
k_hz, K_xtal, and the exponent are not retuned.

See also Model 9b (1/2, rod) and Model 9c (1/3, needle).
"""
from haem_kinetics.models.model8 import Model8


class Model9a(Model8):
    """
    Model 8 chemistry with v_hz scaled by growing Hz surface area.

    At t = 0 the area factor is 1 (same rate as Model 8). Sphere geometry:
    area ∝ amount^{2/3}.
    """

    # Protocol init Hz at V_ref (same fg as [0.018, 0, 0, 0.36]).
    HZ_START_M = 0.36

    # Geometric exponent: sphere = 2/3, rod = 1/2, needle = 1/3
    AREA_EXPONENT = 2.0 / 3.0

    def __init__(self, model_name: str = 'Model 9a'):
        super().__init__(model_name=model_name)

    def _n_hz_start(self) -> float:
        return self.HZ_START_M * self.const.vol_dv

    def _hz_area_factor(self, t: float) -> float:
        n_start = self._n_hz_start()
        n_hz = self._nonneg(self.initial_values['conc_hz']) * self._vol_dv(t)
        if n_start <= 0.0 or n_hz <= 0.0:
            return 0.0
        return (n_hz / n_start) ** self.AREA_EXPONENT

    def _hz_rate(self, t: float):
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        if xtal <= 0.0:
            return 0.0
        return self.const.k_hz * xtal * self._hz_area_factor(t)

    def _d_fe3pp_xtal(self, t):
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        return self._lip_xtal_exchange() - self._hz_rate(t) + self._dilution(t, xtal)

    def _d_hz(self, t):
        hz = self._nonneg(self.initial_values['conc_hz'])
        return self._hz_rate(t) + self._dilution(t, hz)
