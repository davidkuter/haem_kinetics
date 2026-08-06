"""Model 6: Model 5 speciation + dynamic DV volume, depleting host Hb, sigmoidal uptake."""
import pandas as pd

from scipy.integrate import solve_ivp
from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.models.helpers import (
    fg_to_molar,
    lipid_over_aq_ratio,
    molar_to_fg,
    sigmoidal_uptake_rate_per_min,
    vol_dv_fl,
)
from haem_kinetics.components.experimental_data import ExperimentalData


class Model6(KineticsModel):
    """
    Extends Model 5 with mass-conserving transport and time-dependent DV volume.

    Assumptions:
     * States are fg Fe/cell (extensive)
     * Host Hb Fe is an explicit state and is depleted by uptake
     * Uptake rate follows a sigmoidal drain of remaining host Fe (Combrink-like)
     * DV lumen volume V(t) from helpers.vol_dv_fl (Gompertz then collapse)
     * Same aqueous/lipid Fe(III) pools; Hz from lipid pool at k_hz (no φ)
     * Plasmepsins only; [E] scales with lumen volume × fudge
     * Free haem = aqueous + lipid non-Hz Fe(III)

    State vector (fg Fe/cell):
      [conc_hb_dv, conc_fe2pp, conc_fe3pp_aq, conc_fe3pp_lip, conc_hz, conc_hb_host]
    """

    SPECIES = [
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp_aq',
        'conc_fe3pp_lip',
        'conc_hz',
        'conc_hb_host',
    ]
    PROTEASES = ['plm_1', 'plm_2', 'hap', 'plm_4']

    def __init__(self, model_name: str = 'Model 6'):
        super().__init__(model_name=model_name)
        budget = self.const.total_fe_fg_cell
        self._set_initial_conc(init=[1.0, 0.0, 0.0, 0.0, 23.0, budget - 24.0])
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()
        self._k_eff = lipid_over_aq_ratio(
            self.const.vol_fract_lip, self.const.K_partition
        )

    def _enzyme_conc(self, enzyme: str, t: float, vol_fl: float) -> float:
        """Effective enzyme molarity at current lumen volume."""
        # Convert fixed 1 fL reference M to current volume; scale with lumen growth
        growth = max(vol_fl / 3.7, 0.05)
        moles_per_ref = self.const.conc_enzymes[enzyme] * self.const.vol_dv
        return (moles_per_ref / (vol_fl * 1e-15)) * growth * self.const.fudge

    def _calc_enzyme_rate(self, enzyme, conc_hb_m, t, vol_fl):
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = self._enzyme_conc(enzyme, t, vol_fl)
        denom = Km + conc_hb_m
        if denom == 0:
            return 0.0
        return kcat * conc_enzyme / denom

    def _hb_removal_fg(self, t, vol_fl):
        hb_fg = self.initial_values['conc_hb_dv']
        if hb_fg <= 0.0 or vol_fl <= 0.0:
            return 0.0
        conc_hb_tet = fg_to_molar(hb_fg, vol_fl) / 4.0
        deg = 0.0
        for enzyme in self.PROTEASES:
            deg += self._calc_enzyme_rate(enzyme, conc_hb_tet, t, vol_fl)
        removal_m = 4.0 * deg * conc_hb_tet
        return molar_to_fg(removal_m, vol_fl)

    def _exchange_fg(self, vol_fl):
        aq_m = fg_to_molar(self.initial_values['conc_fe3pp_aq'], vol_fl)
        lip_m = fg_to_molar(self.initial_values['conc_fe3pp_lip'], vol_fl)
        flux_m = self.const.k_lipid_exchange * (aq_m - lip_m / self._k_eff)
        return molar_to_fg(flux_m, vol_fl)

    def _set_initial_conc(self, init: List[float]):
        if len(init) != 6:
            raise ValueError(
                'Model6 requires 6 initial values (fg/cell): '
                '[Hb_DV, Fe2, Fe3_aq, Fe3_lip, Hz, Hb_host]'
            )
        for key, val in zip(self.SPECIES, init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        self._set_initial_conc(init=init)
        vol_fl = vol_dv_fl(t)

        uptake = sigmoidal_uptake_rate_per_min(
            t, self.initial_values['conc_hb_host']
        )
        remove = self._hb_removal_fg(t, vol_fl)
        exchange = self._exchange_fg(vol_fl)

        fe2 = self.initial_values['conc_fe2pp']
        fe3_aq = self.initial_values['conc_fe3pp_aq']
        fe3_lip = self.initial_values['conc_fe3pp_lip']

        ox = self.const.k_fe2pp_ox * fe2 * self.const.conc_oxy
        red = self.const.k_fe3pp_red * fe3_aq * self.const.conc_supoxy
        hz_form = self.const.k_hz * fe3_lip

        d_hb = uptake - remove
        d_fe2 = remove + red - ox
        d_fe3_aq = ox - red - exchange
        d_fe3_lip = exchange - hz_form
        d_hz = hz_form
        d_host = -uptake

        return [d_hb, d_fe2, d_fe3_aq, d_fe3_lip, d_hz, d_host]

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            budget = self.const.total_fe_fg_cell
            init = [1.0, 0.0, 0.0, 0.0, 23.0, budget - 24.0]

        self._set_initial_conc(init)
        # Soft check on Fe budget
        tot = sum(init)
        if abs(tot - self.const.total_fe_fg_cell) > 5.0:
            # Renormalise host so total matches budget
            dv_fe = sum(init[:-1])
            init = list(init)
            init[-1] = max(self.const.total_fe_fg_cell - dv_fe, 0.0)
            self._set_initial_conc(init)

        self.solution = solve_ivp(self._integrate, t, init, method='BDF', **kwargs)
        self.time = 16 + self.solution.t / 60
        self.concentrations = pd.DataFrame(
            self.solution.y, columns=self.time, index=list(self.initial_values.keys())
        ).T
        # Already fg/cell — no molar conversion

        if plot:
            self._plot(
                save_file=plot,
                title=self.model_name,
                exp_data=self.exp_data,
                free_haem_cols=['conc_fe3pp_aq', 'conc_fe3pp_lip'],
                columns=['conc_hb_dv', 'conc_hz', 'conc_fe3pp_free'],
                plot_total_fe=True,
            )
