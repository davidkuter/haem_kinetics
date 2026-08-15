import pandas as pd

from typing import List, Optional

from haem_kinetics.models.model3 import Model3
from haem_kinetics.models.native_hb import globin_haem_release_rate, native_tetramer_rate
from haem_kinetics.components.fit_metrics import score_fractionation


class Model4a(Model3):
    """
    Model 3 + Goldberg ordered pathway (two DV Hb pools).

    Native tetramer is nicked by PM I / PM II / FP-2; nicked globin still
    scores as assay Hb until peptide-MM haem release by all six proteases.
    """

    DV_SPECIES = [
        'conc_hb_dv',
        'conc_hb_globin',
        'conc_fe2pp',
        'conc_fe3pp',
        'conc_hz',
    ]

    def __init__(self, model_name: str = 'Model 4a'):
        super().__init__(model_name=model_name)
        # Model1.__init__ passes a 4-DV vector; remap onto the five-pool state.
        self._set_initial_conc(
            init=[0.005, 0.0, 0.0, 0.0, 0.0, self._full_hb_rbc_m()]
        )

    def _nick_rate(self, t):
        native = self._nonneg(self.initial_values['conc_hb_dv'])
        return native_tetramer_rate(self, native, t)

    def _globin_release(self, t):
        globin = self._nonneg(self.initial_values['conc_hb_globin'])
        return globin_haem_release_rate(self, globin, t)

    def _d_hb_dv(self, t):
        hb = self._nonneg(self.initial_values['conc_hb_dv'])
        return self._uptake_dv(t) - self._nick_rate(t) + self._dilution(t, hb)

    def _d_hb_globin(self, t):
        globin = self._nonneg(self.initial_values['conc_hb_globin'])
        return self._nick_rate(t) - self._globin_release(t) + self._dilution(t, globin)

    def _d_fe2pp(self, t):
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        form = self._globin_release(t) + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate() + self._dilution(t, fe2)

    def _set_initial_conc(self, init: List[float]):
        init = self._expand_init(list(init))
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model4a requires 5 DV values [Hb_native, globin, Fe2, Fe3, Hz] '
                'or the 4-value Model 3 init [Hb, Fe2, Fe3, Hz] (globin=0)'
            )
        for key, val in zip(self.DV_SPECIES + [self.HOST_KEY], init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        raw = [float(x) for x in init]
        self._set_initial_conc(init=raw)
        uptake = self._uptake_dv(t)
        return [
            self._d_hb_dv(t),
            self._d_hb_globin(t),
            self._d_fe2pp(t),
            self._d_fe3pp(t),
            self._d_hz(t),
            self._d_host_from_dv_uptake(uptake, t),
        ]

    def _expand_init(self, init: Optional[List[float]]) -> List[float]:
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]
        if len(init) == 4:
            return [init[0], 0.0, init[1], init[2], init[3]]
        return list(init)

    def _prepare_y0(self, init: List[float], t0: float = 0.0) -> List[float]:
        return super()._prepare_y0(self._expand_init(init), t0=t0)

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        init = self._expand_init(init)
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
                columns=['conc_hb_dv', 'conc_hb_globin', 'conc_hz', 'conc_fe3pp'],
                plot_total_fe=True,
            )

    def score_vs_experiment(self, exp_data=None, free_haem_cols=None):
        data = exp_data if exp_data is not None else getattr(self, 'exp_data', None)
        df = self.concentrations.copy()
        df['conc_hb_assay'] = df['conc_hb_dv'] + df['conc_hb_globin']
        self.fit_metrics = score_fractionation(
            df,
            data,
            species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp', 'Hz': 'conc_hz'},
            free_haem_cols=free_haem_cols,
        )
        return self.fit_metrics
