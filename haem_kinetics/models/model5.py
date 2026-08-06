"""Model 5: corrected aqueous/lipid Fe(III) pools; Hz from lipid pool without φ penalty."""
import pandas as pd

from scipy.integrate import solve_ivp
from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.models.helpers import fraction_exp_growth, lipid_over_aq_ratio
from haem_kinetics.components.experimental_data import ExperimentalData


class Model5(KineticsModel):
    """
    Successor to Model 3 with corrected lipid–haemozoin chemistry.

    Assumptions:
     * Exponential Hb transport and enzyme growth (same f_exp as Model 3)
     * Plasmepsin degradation with fudge multiplying [E]
     * Fe(III) split into aqueous and lipid-associated pools with kinetic exchange
     * Haemozoin forms from the lipid pool at k_hz (NOT multiplied by φ)
     * Free haem for comparison = aqueous + lipid non-Hz Fe(III)
     * Fixed DV volume; host Hb only corrected once at t0 (same limitation as Model 3)

    State vector (M, DV basis):
      [conc_hb_dv, conc_fe2pp, conc_fe3pp_aq, conc_fe3pp_lip, conc_hz]
    """

    SPECIES = [
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp_aq',
        'conc_fe3pp_lip',
        'conc_hz',
    ]

    def __init__(self, model_name: str = 'Model 5'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(init=[0.005, 0.0, 0.0, 0.0, 0.0])
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()
        self._k_eff = lipid_over_aq_ratio(
            self.const.vol_fract_lip, self.const.K_partition
        )

    def _calc_enzyme_rate(self, enzyme, conc_hb_dv, t):
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = (
            fraction_exp_growth(t)
            * self.const.conc_enzymes[enzyme]
            * self.const.fudge
        )
        denom = Km + conc_hb_dv
        if denom == 0:
            return 0.0
        return kcat * conc_enzyme / denom

    def _hb_removal(self, t):
        conc_hb = self.initial_values['conc_hb_dv'] / 4
        deg = 0.0
        for enzyme in ['plm_1', 'plm_2', 'hap', 'plm_4']:
            deg += self._calc_enzyme_rate(enzyme, conc_hb, t)
        return 4 * deg * conc_hb

    def _exchange_rate(self):
        """Net flux aqueous -> lipid (M/min), DV-referenced concentrations."""
        aq = self.initial_values['conc_fe3pp_aq']
        lip = self.initial_values['conc_fe3pp_lip']
        # At equilibrium lip = k_eff * aq
        return self.const.k_lipid_exchange * (aq - lip / self._k_eff)

    def _d_hb_dv(self, t):
        tot_hb_conc = self.const.conc_hb_rbc * self.const.vol_rbc / self.const.vol_dv
        form = fraction_exp_growth(t) * tot_hb_conc
        return form - self._hb_removal(t)

    def _d_fe2pp(self, t):
        form = self._hb_removal(t) + (
            self.const.k_fe3pp_red
            * self.initial_values['conc_fe3pp_aq']
            * self.const.conc_supoxy
        )
        remove = (
            self.const.k_fe2pp_ox
            * self.initial_values['conc_fe2pp']
            * self.const.conc_oxy
        )
        return form - remove

    def _d_fe3pp_aq(self):
        form = (
            self.const.k_fe2pp_ox
            * self.initial_values['conc_fe2pp']
            * self.const.conc_oxy
        )
        remove = (
            self.const.k_fe3pp_red
            * self.initial_values['conc_fe3pp_aq']
            * self.const.conc_supoxy
        ) + self._exchange_rate()
        return form - remove

    def _d_fe3pp_lip(self):
        form = self._exchange_rate()
        remove = self.const.k_hz * self.initial_values['conc_fe3pp_lip']
        return form - remove

    def _d_hz(self):
        return self.const.k_hz * self.initial_values['conc_fe3pp_lip']

    def _set_initial_conc(self, init: List[float]):
        if len(init) != 5:
            raise ValueError('Model5 requires 5 initial values: '
                             '[Hb_DV, Fe2, Fe3_aq, Fe3_lip, Hz] (M)')
        for key, val in zip(self.SPECIES, init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        self._set_initial_conc(init=init)
        return [
            self._d_hb_dv(t),
            self._d_fe2pp(t),
            self._d_fe3pp_aq(),
            self._d_fe3pp_lip(),
            self._d_hz(),
        ]

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0, 0.0]

        self._set_initial_conc(init)
        tot_init = sum(self.initial_values.values())
        self.const.conc_hb_rbc = self.const.conc_hb_rbc - (
            tot_init * self.const.vol_dv / self.const.vol_rbc
        )

        self.solution = solve_ivp(self._integrate, t, init, **kwargs)
        self.time = 16 + self.solution.t / 60
        self.concentrations = pd.DataFrame(
            self.solution.y, columns=self.time, index=list(self.initial_values.keys())
        ).T
        self.concentrations = self._molar_to_fgcell(self.concentrations)

        if plot:
            self._plot(
                save_file=plot,
                title=self.model_name,
                exp_data=self.exp_data,
                free_haem_cols=['conc_fe3pp_aq', 'conc_fe3pp_lip'],
                columns=['conc_hb_dv', 'conc_hz', 'conc_fe3pp_free'],
                plot_total_fe=True,
            )
