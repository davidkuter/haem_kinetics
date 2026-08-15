"""Model 10: what-if retune of Model 9 knobs vs Garnie Dd2.

Not a mechanistic ladder step. k_release and K_xtal are chosen so the
existing ODEs can be compared to the assay; they are not Klemba or NLB
geometry. f_exp stays at the Model 2b schedule (observed DV-inventory
delivery) — do not retune a, b here.

Inherits Model 9 crystal-area growth. Plateau k_htv_release is scaled;
Model 6 still multiplies by s_PM(t), so k_10(t) = 0.4774 · k_Klemba · s_PM(t).
Do not refit the scale or K_xtal.
"""
import math

from haem_kinetics.models.model9 import Model9


# Least-squares scale so (then) Model 7 Hb / s matches Garnie Dd2 (n_HTV ∝ 1/k).
# Chosen on the Model 2a inventory; not refit. Plateau t½ ≈ 41.9 min
# (Klemba bound is t½ < 20 min).
HTV_RELEASE_K_SCALE = 0.4774

# Grid of K_xtal with that k_release: 0.10 gives the lowest
# Hm χ²_red among {0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.24}.
K_XTAL_WHATIF = 0.10


class Model10(Model9):
    """
    Model 9 chemistry with Garnie-tuned k_release and K_xtal.

    Diagnostic only: does the interfacial + area topology have enough freedom
    to approach Dd2 Hb and Hm if those two knobs are freed? Uptake stays
    Model 2b f_exp. Lysis still follows s_PM(t) from Model 6. Area growth
    is inherited from Model 9.
    """

    def __init__(self, model_name: str = 'Model 10'):
        super().__init__(model_name=model_name)
        k_klemba = math.log(2) / 20.0
        self.const.k_htv_release = k_klemba * HTV_RELEASE_K_SCALE
        self.const.K_xtal = K_XTAL_WHATIF
