"""Model 7: Model 6 + falcipain-2/3 haemoglobin degradation."""
from haem_kinetics.models.model6 import Model6


class Model7(Model6):
    """
    Same transport, volume, and Fe(III)/Hz scheme as Model 6, with falcipains.

    Adds falcipain-2 and falcipain-3 to the Hb degradation sum. Provisional
    kcat/Km and PaxDB-style abundances are in Constants; refine against
    literature before quantitative claims. `fudge` remains 1.0 (no scaling).

    State vector (fg Fe/cell) — identical to Model 6:
      [conc_hb_dv, conc_fe2pp, conc_fe3pp_aq, conc_fe3pp_lip, conc_hz, conc_hb_host]
    """

    PROTEASES = ['plm_1', 'plm_2', 'hap', 'plm_4', 'fp_2', 'fp_3']

    def __init__(self, model_name: str = 'Model 7'):
        super().__init__(model_name=model_name)
        self.const.fudge = 1.0
