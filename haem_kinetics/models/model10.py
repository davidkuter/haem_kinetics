"""Model 10: Model 9a + full amount encoding for Fe species.

All haem-iron species downstream of Hb_dv are amount-encoded (not concentrated
by aqueous lumen collapse). The DV volume schedule affects only Hb_dv (the
protease-accessible lumen pool).

Mechanistic basis:
- NLBs are lipid droplets that don't shrink with aqueous lumen (Jackson 2004)
- Hz crystals are solid, excluded from pHrodo volume (Garnie 2025)
- Fe(II) and Fe(III) quickly equilibrate with NLBs / are bound to membranes,
  so effectively don't see the aqueous volume collapse

This tests whether the late-phase Hm/Hb dip in Model 9a is caused by the
concentration effect during DV collapse.
"""
from haem_kinetics.models.model9a import Model9a


class Model10(Model9a):
    """
    Model 9a + amount encoding for all Fe species.

    Fe2, all Fe(III) species, and Hz are amount-encoded (not concentrated
    by aqueous lumen collapse). Only Hb_dv remains a true lumen species.
    This removes the artificial late-phase concentration spike.
    """

    AMOUNT_SPECIES = [
        'conc_hb_htv',       # inner-vesicle cargo (already amount in Model 5+)
        'conc_fe2pp',        # Fe(II) — bound / quickly partitions
        'conc_fe3pp_aq',     # aqueous Fe(III)
        'conc_fe3pp_lip',    # bulk NLB lipid
        'conc_fe3pp_xtal',   # NLB-water interface
        'conc_hz',           # haemozoin crystals
    ]

    def __init__(self, model_name: str = 'Model 10'):
        super().__init__(model_name=model_name)

    def _vol_ratio(self, t: float) -> float:
        """Ratio V_DV(t) / V_ref for lumen→amount interface."""
        return self._vol_dv(t) / self.const.vol_dv

    def _n_hz_amount(self) -> float:
        """Hz amount (mol) from the amount-encoded species."""
        return self._nonneg(self.initial_values['conc_hz']) * self.const.vol_dv

    def _hz_area_factor(self, t: float) -> float:
        """Area factor for crystal-area growth.

        For Model 10, Hz is an amount species (M at V_ref), so the amount is
        conc_hz * V_ref, not conc_hz * V_DV(t).
        """
        n_start = self._n_hz_start()
        n_hz = self._n_hz_amount()
        if n_start <= 0.0 or n_hz <= 0.0:
            return 0.0
        return (n_hz / n_start) ** self.AREA_EXPONENT

    def _d_fe2pp(self, t):
        """Fe(II) derivative — no dilution (amount species).
        
        Fe2 is produced by digestion (from Hb_dv lumen pool). The digestion
        rate is scaled by V_DV(t)/V_ref to convert from lumen to amount basis,
        preserving mass balance.
        """
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        vol_ratio = self._vol_ratio(t)
        form = self._hb_removal(t) * vol_ratio + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp_aq'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate()  # No dilution

    def _d_fe3pp_aq(self, t):
        """Aqueous Fe(III) derivative — no dilution (amount species)."""
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        remove = (
            self.const.k_fe3pp_red * aq * self.const.conc_supoxy
        ) + self._exchange_rate()
        return self._ox_rate() - remove  # No dilution

    def _d_fe3pp_lip(self, t):
        """Lipid Fe(III) derivative — no dilution (amount species)."""
        return self._exchange_rate() - self._lip_xtal_exchange()

    def _d_fe3pp_xtal(self, t):
        """Interfacial Fe(III) derivative — no dilution (amount species)."""
        return self._lip_xtal_exchange() - self._hz_rate(t)

    def _d_hz(self, t):
        """Hz derivative — no dilution (amount species)."""
        return self._hz_rate(t)
