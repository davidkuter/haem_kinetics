"""Model 11: Model 9a + two-compartment volume model.

Separates the aqueous DV lumen (V_aq = V_DV(t), collapses 36-46h) from the
NLB compartment (V_nlb = constant). This tests whether separating the volumes
changes the late-phase dynamics.

Compartments:
- **Aqueous lumen** (V_aq = V_DV(t)): Hb_dv, Fe2, Fe3_aq — dilute during collapse
- **NLB compartment** (V_nlb = 1 fL constant): Fe3_lip, Fe3_xtal, Hz — no dilution

When Fe moves from aq → lip, the flux (mol/min) is the same, but the
concentration change is scaled by the volume ratio V_aq/V_nlb.

Mechanistic basis:
- NLBs are lipid droplets that don't shrink with aqueous lumen (Jackson 2004)
- Hz crystals are solid, excluded from pHrodo volume (Garnie 2025)
- The aq↔lip exchange occurs at the NLB-water interface
"""
from haem_kinetics.models.model9a import Model9a


class Model11(Model9a):
    """
    Model 9a + two-compartment volumes (V_aq shrinks, V_nlb constant).

    Aqueous species (Fe2, Fe3_aq) use V_DV(t) and have dilution terms.
    NLB species (Fe3_lip, Fe3_xtal, Hz) use V_nlb (constant) with no dilution.
    Fluxes between compartments are scaled by the volume ratio.
    """

    # NLB species use V_nlb instead of V_DV(t) - includes Fe3_aq since it
    # rapidly equilibrates with the NLB phase
    NLB_SPECIES = ['conc_fe3pp_aq', 'conc_fe3pp_lip', 'conc_fe3pp_xtal', 'conc_hz']

    # For fg conversion: HTV and NLB species use V_ref
    AMOUNT_SPECIES = [
        'conc_hb_htv',       # inner-vesicle cargo (already amount in Model 5+)
        'conc_fe3pp_aq',     # rapidly equilibrates with NLB
        'conc_fe3pp_lip',    # NLB compartment
        'conc_fe3pp_xtal',   # NLB compartment
        'conc_hz',           # NLB compartment (solid)
    ]

    def __init__(self, model_name: str = 'Model 11'):
        super().__init__(model_name=model_name)
        # NLB volume: constant, approximately 1 fL (same as V_ref)
        self._vol_nlb = self.const.vol_dv  # 1e-15 L
        # Minimum volume for enzyme concentration (cap how much [E] can rise)
        # Set to 100% of V_ref (1 fL) — enzymes don't concentrate at all during collapse
        # because they're membrane-associated, not freely soluble
        self._vol_e_min = self.const.vol_dv  # 1 fL

    def _vol_for_enzyme(self, t: float) -> float:
        """Effective volume for enzyme concentration, capped during collapse.
        
        Enzymes are membrane-associated; their effective concentration doesn't
        rise as sharply as if they were freely soluble in the shrinking lumen.
        """
        v_dv = self._vol_dv(t)
        return max(v_dv, self._vol_e_min)

    def _enzyme_conc(self, enzyme: str, t: float) -> float:
        """PaxDB enzyme amount with capped volume during collapse.
        
        Overrides base._enzyme_conc to use capped volume, preventing
        excessive enzyme concentration during late-phase DV collapse.
        """
        from haem_kinetics.models.helpers import garnie_pm_amount_scale
        v_eff = self._vol_for_enzyme(t)
        ppm = self.const.ppm_enzymes[enzyme]
        base_conc = self.const._dv_ppm_to_molar(ppm, vol_dv=v_eff)
        return garnie_pm_amount_scale(t) * base_conc

    def _vol_ratio_aq_nlb(self, t: float) -> float:
        """Ratio V_aq(t) / V_nlb for aq→NLB flux conversion."""
        return self._vol_dv(t) / self._vol_nlb

    def _n_hz_amount(self) -> float:
        """Hz amount (mol) from the NLB-encoded species."""
        return self._nonneg(self.initial_values['conc_hz']) * self._vol_nlb

    def _hz_area_factor(self, t: float) -> float:
        """Area factor for crystal-area growth.

        Hz is in the NLB compartment at V_nlb basis, so amount = conc * V_nlb.
        """
        n_start = self._n_hz_start()
        n_hz = self._n_hz_amount()
        if n_start <= 0.0 or n_hz <= 0.0:
            return 0.0
        return (n_hz / n_start) ** self.AREA_EXPONENT

    # --- Aqueous lumen species (V_aq = V_DV(t), with dilution) ---

    def _d_fe2pp(self, t):
        """Fe(II) derivative — aqueous lumen species with dilution."""
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        # Fe3 reduction term uses Fe3_aq (also in aqueous lumen)
        form = self._hb_removal(t) + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp_aq'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate() + self._dilution(t, fe2)

    def _d_fe3pp_aq(self, t):
        """Aqueous Fe(III) derivative — NLB compartment, no dilution.
        
        Fe3_aq rapidly equilibrates with NLB, so it effectively doesn't see
        the aqueous lumen collapse. Oxidation produces Fe3 from lumen Fe2.
        """
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        vol_ratio = self._vol_ratio_aq_nlb(t)
        # Oxidation from Fe2 (lumen) → Fe3_aq (NLB): scale by volume ratio
        ox_flux = self._ox_rate() * vol_ratio
        remove = (
            self.const.k_fe3pp_red * aq * self.const.conc_supoxy
        ) + self._exchange_rate()
        return ox_flux - remove  # No dilution

    # --- NLB compartment species (V_nlb = constant, no dilution) ---

    def _d_fe3pp_lip(self, t):
        """Lipid Fe(III) derivative — NLB compartment, no dilution.
        
        Both Fe3_aq and Fe3_lip are now in the NLB compartment (same V_nlb basis).
        """
        return self._exchange_rate() - self._lip_xtal_exchange()  # No dilution

    def _d_fe3pp_xtal(self, t):
        """Interfacial Fe(III) derivative — NLB compartment, no dilution."""
        return self._lip_xtal_exchange() - self._hz_rate(t)  # No dilution

    def _d_hz(self, t):
        """Hz derivative — NLB compartment (solid), no dilution."""
        return self._hz_rate(t)  # No dilution

    def _concentrations_to_fgcell(self, df):
        """Convert ODE outputs to fg Fe/cell.
        
        - Aqueous lumen species (Fe2, Fe3_aq, Hb_dv) use V_DV(t)
        - NLB species (Fe3_lip, Fe3_xtal, Hz) use V_nlb (= V_ref)
        - HTV uses V_ref (already in AMOUNT_SPECIES)
        - Host uses vol_rbc
        """
        t_min = (df.index.to_numpy(dtype=float) - 16.0) * 60.0
        vols = [self._vol_dv(float(t)) for t in t_min]
        vol_series = pd.Series(vols, index=df.index)
        factor = (10 ** 15) * 55.85
        amount = set(self.AMOUNT_SPECIES)
        out = pd.DataFrame(index=df.index)
        for col in df.columns:
            if col == self.HOST_KEY:
                out[col] = df[col] * self.const.vol_rbc * factor
            elif col in amount:
                out[col] = df[col] * self._vol_nlb * factor
            else:
                out[col] = df[col] * vol_series * factor
        return out


# Need pandas import for _concentrations_to_fgcell
import pandas as pd
