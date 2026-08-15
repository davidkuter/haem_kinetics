"""Model 9c: Model 8 + crystal-area growth (extreme elongation).

v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{1/3}

Exponent 1/3 is the limiting case where area grows as length only (highly
elongated needle, area ∝ V^{1/3} when diameter is negligible).

See also Model 9a (2/3, sphere) and Model 9b (1/2, rod).
"""
from haem_kinetics.models.model9a import Model9a


class Model9c(Model9a):
    """
    Model 9a with extreme-elongation exponent 1/3.

    At t = 0 the area factor is 1. Needle geometry: area ∝ amount^{1/3}.
    """

    AREA_EXPONENT = 1.0 / 3.0

    def __init__(self, model_name: str = 'Model 9c'):
        super().__init__(model_name=model_name)
