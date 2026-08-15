import pandas as pd

from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.components.experimental_data import ExperimentalData


class Model1(KineticsModel):
    """
    Baseline full speciation: linear uptake, PMs + FP2/3, Fe(II) oxidation, first-order Hz.
    Host [Hb]_RBC is depleted by uptake. Concentrations are kept non-negative in the RHS.

    variable_dv_volume is shared bookkeeping (dilution, [E] = n_E/V(t),
    fg = C·V), not this model's mechanistic change.
    """

    DV_SPECIES = ['conc_hb_dv', 'conc_fe2pp', 'conc_fe3pp', 'conc_hz']

    def __init__(self, model_name: str = 'Model 1'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(init=[0.005, 0.0, 0.0, 0.0, self._full_hb_rbc_m()])
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()

    def _calc_enzyme_rate(self, enzyme, conc_hb_tetramer, t):
        if conc_hb_tetramer <= 0.0:
            return 0.0
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = self._enzyme_conc(enzyme, t)
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

    def _uptake_dv(self, t):
        """Linear host→DV appearance (M/min on V_DV(t)).

        `k_hb_trans` was defined as a concentration rate at the 1 fL reference;
        scaling by V_ref/V(t) keeps the mole delivery independent of lumen size.
        """
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        return self.const.k_hb_trans * host * self.const.vol_dv / self._vol_dv(t)

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

    def _d_hb_dv(self, t):
        hb = self._nonneg(self.initial_values['conc_hb_dv'])
        return self._uptake_dv(t) - self._hb_removal(t) + self._dilution(t, hb)

    def _d_fe2pp(self, t):
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        form = self._hb_removal(t) + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate() + self._dilution(t, fe2)

    def _d_fe3pp(self, t):
        fe3 = self._nonneg(self.initial_values['conc_fe3pp'])
        remove = (
            self.const.k_fe3pp_red * fe3 * self.const.conc_supoxy
        ) + self._hz_rate()
        return self._ox_rate() - remove + self._dilution(t, fe3)

    def _d_hz(self, t):
        hz = self._nonneg(self.initial_values['conc_hz'])
        return self._hz_rate() + self._dilution(t, hz)

    def _set_initial_conc(self, init: List[float]):
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                f'{type(self).__name__} requires 4 DV values [Hb_DV, Fe2, Fe3, Hz], '
                'optionally with host [Hb]_RBC appended'
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
            self._d_fe3pp(t),
            self._d_hz(t),
            self._d_host_from_dv_uptake(uptake, t),
        ]
        return dydt

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]
        t0 = float(t[0]) if t is not None else 0.0
        y0 = self._prepare_y0(init, t0=t0)

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
