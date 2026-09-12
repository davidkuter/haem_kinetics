import math

from haem_kinetics.models.model14a import Model14a


class Model14b(Model14a):
    """
    Model 14a with a slower, provisional inner-vesicle lysis time.

    Model 14a releases at the Klemba rate (t½ ≈ 20 min) but then runs low at
    20–29 h: a fast pool cannot hold the observed ~2 fg standing Hb on the low
    early uptake. Klemba's t½ < 20 min, however, bounds **plasmepsin
    biosynthesis / trafficking**, not necessarily inner-vesicle membrane lysis,
    which may be a slower, separate process. Model 14b tests that: release stays
    decoupled from s_PM (as 14a) but at a slower provisional t½ of 30 min
    (~1.5× the Klemba proxy).

    Result: the larger standing pool matches the early points well (20 h ≈ exp)
    but, because n_HTV ≈ v_up / k_release and late uptake is large, a *constant*
    slow rate then overshoots the late assay Hb (~4 fg at 44 h vs ~2 fg). The
    14a↔14b pair brackets the constant-rate hypothesis: 20 min fits late and
    misses early, 30 min fits early and misses late. No single constant rate
    fits both, so inner-vesicle lysis must **accelerate** through development —
    which is the direction Model 6's s_PM coupling encoded (only its early
    magnitude, from the low/noisy blot, was too extreme).

    Provisional: 30 min is not a measured vesicle-lysis time; it is a labelled
    exploration of the timescale, not a fitted knot.
    """

    # Provisional inner-vesicle lysis half-life (min); Klemba's <20 min bounds
    # PM trafficking, treated here as separate from vesicle-membrane lysis.
    HTV_LYSIS_T_HALF_MIN = 30.0

    def __init__(self, model_name: str = 'Model 14b'):
        super().__init__(model_name=model_name)

    def _k_release(self, t):
        return math.log(2.0) / self.HTV_LYSIS_T_HALF_MIN
