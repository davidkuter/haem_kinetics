import pandas as pd

from typing import List, Optional

from haem_kinetics.models.model4b import Model4b
from haem_kinetics.components.fit_metrics import score_fractionation


class Model5(Model4b):
    """
    Model 4b + inaccessible HTV / inner-vesicle Hb cargo.

    f_exp delivers into conc_hb_htv (not the protease-accessible lumen).
    First-order release feeds Model 4b lumen chemistry (native-competent
    enzymes, lumped haem release). Assay Hb is HTV + lumen Hb.
    """

    AMOUNT_SPECIES = ['conc_hb_htv']
    DV_SPECIES = [
        'conc_hb_htv',
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp',
        'conc_hz',
    ]

    def __init__(self, model_name: str = 'Model 5'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(
            init=[0.005, 0.0, 0.0, 0.0, 0.0, self._full_hb_rbc_m()]
        )

    def _htv_release_lumen(self, t):
        """Cargo appearance in the lumen (M/min on current V_DV)."""
        htv = self._nonneg(self.initial_values['conc_hb_htv'])
        if htv <= 0.0:
            return 0.0
        return self.const.k_htv_release * htv * self.const.vol_dv / self._vol_dv(t)

    def _d_hb_htv(self, t):
        htv = self._nonneg(self.initial_values['conc_hb_htv'])
        appear = self._uptake_dv(t) * self._vol_dv(t) / self.const.vol_dv
        return appear - self.const.k_htv_release * htv

    def _d_hb_dv(self, t):
        hb = self._nonneg(self.initial_values['conc_hb_dv'])
        return self._htv_release_lumen(t) - self._hb_removal(t) + self._dilution(t, hb)

    def _set_initial_conc(self, init: List[float]):
        init = self._expand_init(list(init))
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model5 requires 5 DV values [HTV, Hb_lumen, Fe2, Fe3, Hz] '
                'or the 4-value Model 3 init [Hb, Fe2, Fe3, Hz] (Hb seeds HTV)'
            )
        for key, val in zip(self.DV_SPECIES + [self.HOST_KEY], init):
            self.initial_values[key] = val

    def _integrate(self, t, init):
        raw = [float(x) for x in init]
        self._set_initial_conc(init=raw)
        uptake = self._uptake_dv(t)
        return [
            self._d_hb_htv(t),
            self._d_hb_dv(t),
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
        self.concentrations['conc_hb_assay'] = (
            self.concentrations['conc_hb_htv'] + self.concentrations['conc_hb_dv']
        )

        if plot:
            self._plot(
                save_file=plot,
                title=self.model_name,
                exp_data=self.exp_data,
                columns=['conc_hb_assay', 'conc_hz', 'conc_fe3pp'],
                plot_total_fe=True,
            )

    def score_vs_experiment(self, exp_data=None, free_haem_cols=None):
        data = exp_data if exp_data is not None else getattr(self, 'exp_data', None)
        df = self.concentrations.copy()
        if 'conc_hb_assay' not in df.columns:
            df['conc_hb_assay'] = df['conc_hb_htv'] + df['conc_hb_dv']
        self.fit_metrics = score_fractionation(
            df,
            data,
            species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp', 'Hz': 'conc_hz'},
            free_haem_cols=free_haem_cols,
        )
        return self.fit_metrics
