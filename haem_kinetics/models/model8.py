import pandas as pd

from typing import List, Optional

from haem_kinetics.models.model7 import Model7
from haem_kinetics.components.fit_metrics import score_fractionation


class Model8(Model7):
    """
    Model 7 + interfacial (crystal-competent) Fe(III).

    Bulk NLB Fe(III) exchanges with an interfacial pool; Hz forms from the
    interface at literature k_hz. Assay Hm is aq + lip + xtal (all still
    non-crystalline Fe(III)PPIX). Assay Hb remains HTV + lumen.
    """

    AMOUNT_SPECIES = ['conc_hb_htv']
    DV_SPECIES = [
        'conc_hb_htv',
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp_aq',
        'conc_fe3pp_lip',
        'conc_fe3pp_xtal',
        'conc_hz',
    ]
    FREE_HAEM_COLS = ['conc_fe3pp_aq', 'conc_fe3pp_lip', 'conc_fe3pp_xtal']

    def __init__(self, model_name: str = 'Model 8'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(
            init=[0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, self._full_hb_rbc_m()]
        )

    def _lip_xtal_exchange(self):
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        return self.const.k_xtal_exchange * (lip - xtal / self.const.K_xtal)

    def _hz_rate(self):
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        if xtal <= 0.0:
            return 0.0
        return self.const.k_hz * xtal

    def _d_fe3pp_lip(self, t):
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        return self._exchange_rate() - self._lip_xtal_exchange() + self._dilution(t, lip)

    def _d_fe3pp_xtal(self, t):
        xtal = self._nonneg(self.initial_values['conc_fe3pp_xtal'])
        return self._lip_xtal_exchange() - self._hz_rate() + self._dilution(t, xtal)

    def _set_initial_conc(self, init: List[float]):
        init = self._expand_init(list(init))
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model8 requires 7 DV values '
                '[HTV, Hb_lumen, Fe2, Fe3_aq, Fe3_lip, Fe3_xtal, Hz] '
                'or the 4-value Model 3 init [Hb, Fe2, Fe3, Hz] '
                '(Hb seeds HTV; Fe3 seeds aq; lip and xtal = 0)'
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
            self._d_fe3pp_aq(t),
            self._d_fe3pp_lip(t),
            self._d_fe3pp_xtal(t),
            self._d_hz(t),
            self._d_host_from_dv_uptake(uptake, t),
        ]

    def _expand_init(self, init: Optional[List[float]]) -> List[float]:
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]
        if len(init) == 4:
            return [init[0], 0.0, init[1], init[2], 0.0, 0.0, init[3]]
        if len(init) == 5:
            return [init[0], init[1], init[2], init[3], 0.0, 0.0, init[4]]
        if len(init) == 6:
            return [init[0], init[1], init[2], init[3], init[4], 0.0, init[5]]
        return list(init)

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
                free_haem_cols=self.FREE_HAEM_COLS,
                columns=['conc_hb_assay', 'conc_hz', 'conc_fe3pp_free'],
                plot_total_fe=True,
            )

    def score_vs_experiment(self, exp_data=None, free_haem_cols=None):
        data = exp_data if exp_data is not None else getattr(self, 'exp_data', None)
        df = self.concentrations.copy()
        if 'conc_hb_assay' not in df.columns:
            df['conc_hb_assay'] = df['conc_hb_htv'] + df['conc_hb_dv']
        cols = free_haem_cols or self.FREE_HAEM_COLS
        self.fit_metrics = score_fractionation(
            df,
            data,
            species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp_free', 'Hz': 'conc_hz'},
            free_haem_cols=cols,
        )
        return self.fit_metrics
