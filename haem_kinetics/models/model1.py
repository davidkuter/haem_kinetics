import pandas as pd

from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.components.experimental_data import ExperimentalData


class Model1(KineticsModel):
    """
    Baseline full speciation: linear uptake, PMs + FP2/3, Fe(II) oxidation, first-order Hz.
    Host [Hb]_RBC is depleted by uptake. Concentrations are kept non-negative in the RHS.
    """

    DV_SPECIES = ['conc_hb_dv', 'conc_fe2pp', 'conc_fe3pp', 'conc_hz']

    def __init__(self, model_name: str = 'Model 1'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(init=[0.005, 0.0, 0.0, 0.0, self._full_hb_rbc_m()])
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()

    def _calc_enzyme_rate(self, enzyme, conc_hb_tetramer):
        if conc_hb_tetramer <= 0.0:
            return 0.0
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = self.const.conc_enzymes[enzyme]
        denom = Km + conc_hb_tetramer
        if denom <= 0.0:
            return 0.0
        return kcat * conc_enzyme / denom

    def _hb_removal(self):
        conc_hb_tetramer = self._nonneg(self.initial_values['conc_hb_dv']) / 4.0
        if conc_hb_tetramer <= 0.0:
            return 0.0
        deg = 0.0
        for enzyme in ['plm_1', 'plm_2', 'hap', 'plm_4', 'fp_2', 'fp_3']:
            deg += self._calc_enzyme_rate(enzyme, conc_hb_tetramer)
        return 4.0 * deg * conc_hb_tetramer

    def _uptake_dv(self):
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        return self.const.k_hb_trans * host

    def _ox_rate(self):
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        if fe2 <= 0.0:
            return 0.0
        return self.const.k_fe2pp_ox * fe2 * self.const.conc_oxy

    def _hz_rate(self):
        fe3 = self._nonneg(self.initial_values['conc_fe3pp'])
        if fe3 <= 0.0:
            return 0.0
        return self.const.k_hz * fe3

    def _d_hb_dv(self):
        return self._uptake_dv() - self._hb_removal()

    def _d_fe2pp(self):
        form = self._hb_removal() + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate()

    def _d_fe3pp(self):
        fe3 = self._nonneg(self.initial_values['conc_fe3pp'])
        remove = (
            self.const.k_fe3pp_red * fe3 * self.const.conc_supoxy
        ) + self._hz_rate()
        return self._ox_rate() - remove

    def _d_hz(self):
        return self._hz_rate()

    def _set_initial_conc(self, init: List[float]):
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model1 requires 4 DV values [Hb_DV, Fe2, Fe3, Hz], '
                'optionally with host [Hb]_RBC appended'
            )
        for key, val in zip(self.DV_SPECIES + [self.HOST_KEY], init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        raw = [float(x) for x in init]
        self._set_initial_conc(init=raw)
        uptake = self._uptake_dv()
        dydt = [
            self._d_hb_dv(),
            self._d_fe2pp(),
            self._d_fe3pp(),
            self._d_hz(),
            self._d_host_from_dv_uptake(uptake),
        ]
        return dydt

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]
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
                columns=['conc_hb_dv', 'conc_hz', 'conc_fe3pp'],
                plot_total_fe=True,
            )
