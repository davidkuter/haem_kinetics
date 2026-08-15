"""Model 6: Model 5 proteases + crystal-competent Fe(III) pool."""
import pandas as pd

from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.models.helpers import (
    enzyme_logistic_scale,
    fraction_exp_growth,
    lipid_over_aq_ratio,
)
from haem_kinetics.components.experimental_data import ExperimentalData


class Model6(KineticsModel):
    """
    One change vs Model 5: Fe(III) crystallization path.

    Keeps Model 5 proteases (PMs + FP2/3, logistic [E]).

    Adds crystal-competent pool:
      aq ⇄ lip ⇄ xtal → Hz at literature k_hz
    so most lipid-associated Fe can remain assayable free haem.
    """

    variable_dv_volume = False
    DV_SPECIES = [
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp_aq',
        'conc_fe3pp_lip',
        'conc_fe3pp_xtal',
        'conc_hz',
    ]
    PROTEASES = ['plm_1', 'plm_2', 'hap', 'plm_4', 'fp_2', 'fp_3']

    def __init__(self, model_name: str = 'Model 6'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(
            init=[0.005, 0.0, 0.0, 0.0, 0.0, 0.0, self._full_hb_rbc_m()]
        )
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()
        self._k_eff = lipid_over_aq_ratio(
            self.const.vol_fract_lip, self.const.K_partition
        )

    def _calc_enzyme_rate(self, enzyme, conc_hb_tetramer, t):
        if conc_hb_tetramer <= 0.0:
            return 0.0
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = (
            enzyme_logistic_scale(t)
            * self.const.conc_enzymes[enzyme]
        )
        denom = Km + conc_hb_tetramer
        if denom <= 0.0:
            return 0.0
        return kcat * conc_enzyme / denom

    def _hb_removal(self, t):
        conc_hb_tetramer = self._nonneg(self.initial_values['conc_hb_dv']) / 4.0
        if conc_hb_tetramer <= 0.0:
            return 0.0
        deg = 0.0
        for enzyme in self.PROTEASES:
            deg += self._calc_enzyme_rate(enzyme, conc_hb_tetramer, t)
        return 4.0 * deg * conc_hb_tetramer

    def _aq_lip_exchange(self):
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        return self.const.k_lipid_exchange * (aq - lip / self._k_eff)

    def _lip_xtal_exchange(self):
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        return self.const.k_xtal_exchange * (lip - xtal / self.const.K_xtal)

    def _uptake_dv(self, t):
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        tot_hb_conc = host * self.const.vol_rbc / self.const.vol_dv
        return fraction_exp_growth(t) * tot_hb_conc

    def _ox_rate(self):
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        if fe2 <= 0.0:
            return 0.0
        return self.const.k_fe2pp_ox * fe2 * self.const.conc_oxy

    def _hz_rate(self):
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        if xtal <= 0.0:
            return 0.0
        return self.const.k_hz * xtal

    def _d_hb_dv(self, t):
        return self._uptake_dv(t) - self._hb_removal(t)

    def _d_fe2pp(self, t):
        form = self._hb_removal(t) + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp_aq'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate()

    def _d_fe3pp_aq(self):
        remove = (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp_aq'])
            * self.const.conc_supoxy
        ) + self._aq_lip_exchange()
        return self._ox_rate() - remove

    def _d_fe3pp_lip(self):
        return self._aq_lip_exchange() - self._lip_xtal_exchange()

    def _d_fe3pp_xtal(self):
        return self._lip_xtal_exchange() - self._hz_rate()

    def _d_hz(self):
        return self._hz_rate()

    def _set_initial_conc(self, init: List[float]):
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model6 requires 6 DV values '
                '[Hb_DV, Fe2, Fe3_aq, Fe3_lip, Fe3_xtal, Hz], optionally + host'
            )
        for key, val in zip(self.DV_SPECIES + [self.HOST_KEY], init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        raw = [float(x) for x in init]
        self._set_initial_conc(init=raw)
        uptake = self._uptake_dv(t)
        dydt = [
            self._d_hb_dv(t),
            self._d_fe2pp(t),
            self._d_fe3pp_aq(),
            self._d_fe3pp_lip(),
            self._d_fe3pp_xtal(),
            self._d_hz(),
            self._d_host_from_dv_uptake(uptake),
        ]
        return dydt

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        y0 = self._pad_init_with_host(init) if len(init) == len(self.DV_SPECIES) else list(init)

        self.solution = self._solve_ivp(self._integrate, t, y0, **kwargs)
        self.time = 16 + self.solution.t / 60
        self.concentrations = pd.DataFrame(
            self.solution.y, columns=self.time, index=list(self.initial_values.keys())
        ).T
        self.concentrations = self._concentrations_to_fgcell(self.concentrations)

        if plot:
            self._plot(
                save_file=plot,
                title=self.model_name,
                exp_data=self.exp_data,
                free_haem_cols=['conc_fe3pp_aq', 'conc_fe3pp_lip'],
                columns=['conc_hb_dv', 'conc_hz', 'conc_fe3pp_free'],
                plot_total_fe=True,
            )
