from haem_kinetics.models.helpers import garnie_pm_amount_scale
from haem_kinetics.models.model13 import Model13


class Model15a(Model13):
    """
    Model 13 with a two-component inner-vesicle lysis clock.

    The single change from Model 13 is how s_PM enters k_release. Model 6 set
    k_release(t) = k_Klemba · s_PM(t), treating the NF54 plasmepsin blot as if
    it were the lysis rate. Early s_PM ≈ 0.41 then makes t½ ≈ 49 min — slower
    than Klemba's t½ < 20 min bound — and piles HTV cargo so early assay Hb
    sits too high. The *shape* (rise, 29 h dip, rise) is the blot clock and
    is worth keeping; only the early *amount* is wrong.

    Model 15a splits lysis into two processes:
      * constitutive inner-membrane rupture / trafficking (Klemba), which
        does not wait on the blot;
      * a protease-modulated term that still follows s_PM(t), so the 29 h
        dip remains.

        k_release(t) = k_Klemba · [f + (1 − f) · s_PM(t)]

    Late (s_PM → 1) is unchanged from Model 13. Early s_PM ≈ 0.41 gives
    t½ ≈ 28 min instead of 49 min, so standing cargo falls without flattening
    the wiggle. Digestion still uses the full s_PM enzyme schedule.

    Provisional: f = 1/2 (equal constitutive and protease-modulated weights).
    Not fitted to Garnie Hb. The strict reading of Klemba (t½ never > 20 min)
    would force f = 1 and collapse to Model 14a's constant rate.
    """

    # Constitutive fraction of k_Klemba that does not scale with the PM blot.
    RELEASE_CONSTITUTIVE = 0.5

    def __init__(self, model_name: str = 'Model 15a'):
        super().__init__(model_name=model_name)

    def _k_release(self, t):
        s = garnie_pm_amount_scale(t)
        f = self.RELEASE_CONSTITUTIVE
        return self.const.k_htv_release * (f + (1.0 - f) * s)
