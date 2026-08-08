"""Model 4: aqueous/lipid Fe(III) pools; Hz from lipid at k_hz (no φ)."""
import pandas as pd

from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.models.helpers import fraction_exp_growth, lipid_over_aq_ratio
from haem_kinetics.components.experimental_data import ExperimentalData


class Model4(KineticsModel):
    """
    Successor to Model 3 with corrected lipid–haemozoin chemistry.
    Concentrations kept non-negative in the RHS.
    """

    DV_SPECIES = [
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp_aq',
        'conc_fe3pp_lip',
        'conc_hz',
    ]

    def __init__(self, model_name: str = 'Model 4'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(init=[0.005, 0.0, 0.0, 0.0, 0.0, self._full_hb_rbc_m()])
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
            fraction_exp_growth(t)
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
        for enzyme in ['plm_1', 'plm_2', 'hap', 'plm_4', 'fp_2', 'fp_3']:
            deg += self._calc_enzyme_rate(enzyme, conc_hb_tetramer, t)
        return 4.0 * deg * conc_hb_tetramer

    def _exchange_rate(self):
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        return self.const.k_lipid_exchange * (aq - lip / self._k_eff)

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
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        if lip <= 0.0:
            return 0.0
        return self.const.k_hz * lip

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
        ) + self._exchange_rate()
        return self._ox_rate() - remove

    def _d_fe3pp_lip(self):
        return self._exchange_rate() - self._hz_rate()

    def _d_hz(self):
        return self._hz_rate()

    def _set_initial_conc(self, init: List[float]):
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model4 requires 5 DV values '
                '[Hb_DV, Fe2, Fe3_aq, Fe3_lip, Hz], optionally + host [Hb]_RBC'
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
            self._d_hz(),
            self._d_host_from_dv_uptake(uptake),
        ]
        return dydt

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0, 0.0]
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
