from haem_kinetics.models.model13 import Model13


class Model14a(Model13):
    """
    Model 13 with inner-vesicle release decoupled from plasmepsin amount.

    The single change from Model 13 is the HTV release clock. Model 6 set
    k_release(t) = k_Klemba · s_PM(t) — "a modelling convenience, not a claim"
    that vesicle lysis tracks protease expression. Because assay Hb is the HTV
    cargo at the quasi-steady value n_HTV ≈ v_up / (k_release · s_PM), and
    Garnie's early plasmepsin blot is low (~0.41 of plateau) and noisy (a
    dropped-outlier dip at 24 h), that coupling collapses the release rate
    exactly where the data is flat and piles the cargo into a spurious 24–26 h
    hump.

    Model 14a states the mechanism explicitly: inner-vesicle membrane lysis
    (delivering Hb to the DV lumen) is a fusion/rupture event, **not** gated by
    plasmepsin abundance — the plasmepsins digest Hb only *after* it reaches the
    lumen. Release is therefore the constant Klemba rate (t½ ≈ 20 min), not
    ∝ s_PM. Digestion still scales with enzyme amount (s_PM); only release is
    decoupled. Everything else is Model 13 (Myburgh conserving uptake,
    upper-range MCHC budget, aqueous-hematin crystallisation, HTV Hb).

    Result: the hump is removed and assay Hb rises smoothly. Because a fast pool
    (t½ ≤ 20 min) cannot retain the observed ~2 fg on the low early uptake
    (~2 fg/h), it runs slightly *low* at 20–29 h — an honest residual pointing
    at the inner-vesicle lysis timescale (Model 14b) rather than a fudge.
    """

    def __init__(self, model_name: str = 'Model 14a'):
        super().__init__(model_name=model_name)

    def _k_release(self, t):
        return self.const.k_htv_release
