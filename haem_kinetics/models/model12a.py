from haem_kinetics.models.model7 import Model7


class Model12a(Model7):
    """
    Myburgh (2023) Hm/Hz speciation on our HTV Hb machinery.

    Same states as Model 7 (aqueous ⇄ lipid Fe(III) equilibrium partition,
    K_partition = 398, vol_fract_lip = 0.016). The single mechanistic change
    from Model 7 is the crystallisation substrate: Hz nucleates from *aqueous*
    hematin at the lipid–water interface (Myburgh Model 3, eq. 4.53), not from
    the bulk lipid pool. Bulk NLB lipid is a reservoir (Egan lipid-mediated
    β-haematin), so only the small aqueous fraction (~13 %) is crystal-competent
    and the lipid buffer keeps assayable free haem high.

    Hb uptake is ours (f_exp × remaining host, Model 2b). Assay Hb = HTV + lumen,
    assay Hm = Fe(II) + Fe(III)_aq + Fe(III)_lip. Literature constants only;
    nothing is retuned relative to Model 7.
    """

    def __init__(self, model_name: str = 'Model 12a'):
        super().__init__(model_name=model_name)

    def _hz_rate(self):
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        if aq <= 0.0:
            return 0.0
        return self.const.k_hz * aq

    def _d_fe3pp_aq(self, t):
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        remove = (
            self.const.k_fe3pp_red * aq * self.const.conc_supoxy
        ) + self._exchange_rate() + self._hz_rate()
        return self._ox_rate() - remove + self._dilution(t, aq)

    def _d_fe3pp_lip(self, t):
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        return self._exchange_rate() + self._dilution(t, lip)
