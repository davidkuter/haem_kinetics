"""Model 9b: Model 8 + crystal-area growth (rod/needle geometry).

v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{1/2}

Exponent 1/2 is for elongated crystals where area grows as length^1 × diameter
(fixed diameter, length ∝ volume). β-haematin is known to form elongated prisms.

See also Model 9a (2/3, sphere) and Model 9c (1/3, extreme elongation).
"""
from haem_kinetics.models.model9a import Model9a


class Model9b(Model9a):
    """
    Model 9a with rod/needle exponent 1/2.

    At t = 0 the area factor is 1. Rod geometry: area ∝ amount^{1/2}.
    """

    AREA_EXPONENT = 1.0 / 2.0

    def __init__(self, model_name: str = 'Model 9b'):
        super().__init__(model_name=model_name)
