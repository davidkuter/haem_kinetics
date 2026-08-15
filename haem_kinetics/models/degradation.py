import pandas as pd

from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.models.helpers import fraction_exp_growth
from haem_kinetics.components.experimental_data import ExperimentalData


class Degradation(KineticsModel):
    """
    Sandbox for Hb uptake / enzymatic release of Fe(II). No Fe(III)/Hz ODEs.
    Host depleted by uptake. Concentrations kept non-negative in the RHS.
    Uses shared variable_dv_volume for molar bookkeeping (same fg conversion
    as Models 1–3). The Hb ODE still omits digestion — diagnostic, not closed Fe.
    """

    DV_SPECIES = ['conc_hb_dv', 'conc_fe2pp']

    def __init__(self, model_name: str = 'Degradation'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(init=[0.0, 0.0, self._full_hb_rbc_m()])
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()

    def _calc_enzyme_rate(self, enzyme, conc_hb_tetramer, t):
        if conc_hb_tetramer <= 0.0:
            return 0.0
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = fraction_exp_growth(t) * self._enzyme_conc(enzyme, t)
        denom = Km + conc_hb_tetramer
        if denom <= 0.0:
            return 0.0
        return kcat * conc_enzyme / denom

    def _hb_removal(self, t):
        conc_hb_tetramer = self._nonneg(self.initial_values['conc_hb_dv']) / 4.0
        if conc_hb_tetramer <= 0.0:
            return 0.0
        deg = 0.0
        for enzyme in ['plm_1', 'plm_2', 'hap', 'plm_4']:
            deg += self._calc_enzyme_rate(enzyme, conc_hb_tetramer, t)
        return 4.0 * deg * conc_hb_tetramer

    def _uptake_dv(self, t):
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        tot_hb_conc = host * self.const.vol_rbc / self._vol_dv(t)
        return fraction_exp_growth(t) * tot_hb_conc

    def _d_hb_dv_kuter(self, t):
        hb = self._nonneg(self.initial_values['conc_hb_dv'])
        return self._uptake_dv(t) + self._dilution(t, hb)

    def _d_fe2pp(self, t):
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        return self._hb_removal(t=t) + self._dilution(t, fe2)

    def _set_initial_conc(self, init: List[float]):
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Degradation requires 2 DV values [Hb_DV, Fe2], optionally + host'
            )
        for key, val in zip(self.DV_SPECIES + [self.HOST_KEY], init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        raw = [float(x) for x in init]
        self._set_initial_conc(init=raw)
        uptake = self._uptake_dv(t)
        dydt = [
            self._d_hb_dv_kuter(t),
            self._d_fe2pp(t),
            self._d_host_from_dv_uptake(uptake, t),
        ]
        return dydt

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0]
        t0 = float(t[0]) if t is not None else 0.0
        y0 = self._prepare_y0(init, t0=t0)

        self.solution = self._solve_ivp(self._integrate, t, y0, **kwargs)
        self.time = 16 + self.solution.t / 60
        self.concentrations = pd.DataFrame(
            self.solution.y, columns=self.time, index=list(self.initial_values.keys())
        ).T
        self.concentrations = self._concentrations_to_fgcell(self.concentrations)
        self.concentrations['conc_hz'] = 0.0
        self.concentrations['conc_hb_dv_obs'] = (
            self.concentrations['conc_hb_dv'] - self.concentrations['conc_fe2pp']
        )

        if plot:
            self._plot(
                save_file=plot,
                title=self.model_name,
                exp_data=self.exp_data,
                columns=['conc_hb_dv_obs', 'conc_hz'],
                plot_total_fe=True,
            )
