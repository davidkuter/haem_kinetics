from haem_kinetics.components.constants import Constants
from haem_kinetics.models.model12b import Model12b


class Model13(Model12b):
    """
    Model 12b (Myburgh empirical uptake, host-conserving) on an upper-range
    red-cell Fe budget.

    The single change from Model 12b is the size of the host Fe pool. Model 12b
    inherits the ladder default budget of ~106 fg/cell, built from **population
    means**: MCHC = 34 g/dL and MCV = 90 fL, i.e. MCH ≈ 30.6 pg Hb/cell. On that
    average cell, Myburgh's empirical exponential (calibrated to the full Dd2 DV
    inventory) drains the ~85 fg host to exactly zero by 44 h; uptake then stops
    abruptly and the standing HTV cargo collapses within its ~20 min release
    half-life → the spurious 44 h cliff.

    Garnie's own Dd2 fractionation shows the average budget is already fully
    accounted for inside the DV by 44 h (Hz 97.7 + Hm 5.8 + Hb 2.0 ≈ 105.5 fg ≈
    106 fg), yet ~2 fg is still Hb-form and turning over — which is only possible
    if delivery is still running, i.e. if these cells carried **upper-range** Hb
    rather than the population mean. Model 13 therefore uses the upper end of the
    clinical MCHC reference range (36 g/dL; range 32–36) at the same 90 fL MCV,
    giving MCH ≈ 32.4 pg and a ~112 fg budget. The host is drawn down but not
    exhausted, so uptake never stops dead and the standing Hb/Hm pools are
    sustained through 44 h — no cliff, no added rate term.

    Provisional: MCHC is a per-donor / per-cell distribution; 36 g/dL is a cited
    upper bound, not a fitted number. Everything else (Myburgh exponential
    uptake, aqueous-hematin crystallisation, HTV Hb, host conservation) is
    exactly Model 12b.
    """

    # Upper end of the clinical MCHC reference range (g/dL); mean used elsewhere
    # is 34. MCV (vol_rbc = 90 fL) is unchanged, so MCH ≈ 32.4 pg (range 27–33).
    MCHC_G_PER_DL = 36.0
    _MEAN_MCHC_G_PER_DL = 34.0

    def __init__(self, model_name: str = 'Model 13'):
        super().__init__(model_name=model_name)
        scale = self.MCHC_G_PER_DL / self._MEAN_MCHC_G_PER_DL
        # Keep the plotted budget and the mass-balance diagnostics consistent
        # with the larger host pool this model integrates.
        self.const.conc_hb_rbc = Constants.compute_conc_hb_rcb() * scale
        self.const.total_fe_fg_cell = self.const.conc_hb_rbc * 90.0 * 55.85

    def _full_hb_rbc_m(self) -> float:
        """Uninfected-RBC haem-equivalent concentration at upper-range MCHC."""
        scale = self.MCHC_G_PER_DL / self._MEAN_MCHC_G_PER_DL
        return Constants.compute_conc_hb_rcb() * scale
