import pandas as pd

from typing import List, Optional

from haem_kinetics.models.model6 import Model6
from haem_kinetics.models.helpers import lipid_over_aq_ratio
from haem_kinetics.components.fit_metrics import score_fractionation


class Model7(Model6):
    """
    Model 6 + aqueous ⇄ lipid Fe(III) partition.

    Hz forms from the lipid pool at literature k_hz (no φ). Assay Hm is
    aqueous + lipid-associated Fe(III). Assay Hb remains HTV + lumen.
    """

    AMOUNT_SPECIES = ['conc_hb_htv']
    DV_SPECIES = [
        'conc_hb_htv',
        'conc_hb_dv',
        'conc_fe2pp',
        'conc_fe3pp_aq',
        'conc_fe3pp_lip',
        'conc_hz',
    ]

    def __init__(self, model_name: str = 'Model 7'):
        super().__init__(model_name=model_name)
        self._k_eff = lipid_over_aq_ratio(
            self.const.vol_fract_lip, self.const.K_partition
        )
        self._set_initial_conc(
            init=[0.005, 0.0, 0.0, 0.0, 0.0, 0.0, self._full_hb_rbc_m()]
        )

    def _exchange_rate(self):
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        return self.const.k_lipid_exchange * (aq - lip / self._k_eff)

    def _hz_rate(self):
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        if lip <= 0.0:
            return 0.0
        return self.const.k_hz * lip

    def _d_fe2pp(self, t):
        fe2 = self._nonneg(self.initial_values['conc_fe2pp'])
        form = self._hb_removal(t) + (
            self.const.k_fe3pp_red
            * self._nonneg(self.initial_values['conc_fe3pp_aq'])
            * self.const.conc_supoxy
        )
        return form - self._ox_rate() + self._dilution(t, fe2)

    def _d_fe3pp_aq(self, t):
        aq = self._nonneg(self.initial_values['conc_fe3pp_aq'])
        remove = (
            self.const.k_fe3pp_red * aq * self.const.conc_supoxy
        ) + self._exchange_rate()
        return self._ox_rate() - remove + self._dilution(t, aq)

    def _d_fe3pp_lip(self, t):
        lip = self._nonneg(self.initial_values['conc_fe3pp_lip'])
        return self._exchange_rate() - self._hz_rate() + self._dilution(t, lip)

    def _set_initial_conc(self, init: List[float]):
        init = self._expand_init(list(init))
        if len(init) == len(self.DV_SPECIES):
            init = self._pad_init_with_host(init)
        if len(init) != len(self.DV_SPECIES) + 1:
            raise ValueError(
                'Model7 requires 6 DV values [HTV, Hb_lumen, Fe2, Fe3_aq, Fe3_lip, Hz] '
                'or the 4-value Model 3 init [Hb, Fe2, Fe3, Hz] (Hb seeds HTV; Fe3 seeds aq)'
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
            self._d_hz(t),
            self._d_host_from_dv_uptake(uptake, t),
        ]

    def _expand_init(self, init: Optional[List[float]]) -> List[float]:
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]
        if len(init) == 4:
            return [init[0], 0.0, init[1], init[2], 0.0, init[3]]
        if len(init) == 5:
            return [init[0], init[1], init[2], init[3], 0.0, init[4]]
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
                free_haem_cols=['conc_fe3pp_aq', 'conc_fe3pp_lip'],
                columns=['conc_hb_assay', 'conc_hz', 'conc_fe3pp_free'],
                plot_total_fe=True,
            )

    def score_vs_experiment(self, exp_data=None, free_haem_cols=None):
        data = exp_data if exp_data is not None else getattr(self, 'exp_data', None)
        df = self.concentrations.copy()
        if 'conc_hb_assay' not in df.columns:
            df['conc_hb_assay'] = df['conc_hb_htv'] + df['conc_hb_dv']
        cols = free_haem_cols or ['conc_fe3pp_aq', 'conc_fe3pp_lip']
        self.fit_metrics = score_fractionation(
            df,
            data,
            species_map={'Hb': 'conc_hb_assay', 'Hm': 'conc_fe3pp_free', 'Hz': 'conc_hz'},
            free_haem_cols=cols,
        )
        return self.fit_metrics
