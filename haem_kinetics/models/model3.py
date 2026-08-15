import pandas as pd

from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.models.helpers import fraction_exp_growth, garnie_dd2_vol_dv_L
from haem_kinetics.components.experimental_data import ExperimentalData


class Model3(KineticsModel):
    """
    Model 2 chemistry with Garnie Dd2 V_DV(t) for molar bookkeeping.

    Uptake stays f_exp (mole delivery unchanged). Enzyme amount is the PaxDB
    count; [E](t) = n_E / V_DV(t). Concentration ODEs include the dilution
    term −C·(dV/dt)/V so C·V is conserved when volume changes.
    """

    DV_SPECIES = ['conc_hb_dv', 'conc_fe2pp', 'conc_fe3pp', 'conc_hz']

    def __init__(self, model_name: str = 'Model 3'):
        super().__init__(model_name=model_name)
        self._set_initial_conc(init=[0.005, 0.0, 0.0, 0.0, self._full_hb_rbc_m()])
        self.exp_data = ExperimentalData()
        self.exp_data.no_drug_dd2()

    def _vol_dv(self, t):
        vol, _dvol = garnie_dd2_vol_dv_L(t)
        if vol <= 0.0:
            raise ValueError(
                f'DV lumen volume is {vol} L at t={t} min; '
                'concentration ODEs are undefined at V≤0 (Garnie collapse reaches 0 at 46 h)'
            )
        return vol

    def _dvol_dv(self, t):
        _vol, dvol = garnie_dd2_vol_dv_L(t)
        return dvol

    def _dilution(self, t, conc):
        """Chain rule for C = n/V: dC/dt includes −C·(dV/dt)/V."""
        vol = self._vol_dv(t)
        return -conc * self._dvol_dv(t) / vol

    def _enzyme_conc(self, enzyme, t):
        ppm = self.const.ppm_enzymes[enzyme]
        return self.const._dv_ppm_to_molar(ppm, vol_dv=self._vol_dv(t))

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
        host = self._nonneg(self.initial_values[self.HOST_KEY])
        if host <= 0.0:
            return 0.0
        tot_hb_conc = host * self.const.vol_rbc / self._vol_dv(t)
        return fraction_exp_growth(t) * tot_hb_conc

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
                'Model3 requires 4 DV values [Hb_DV, Fe2, Fe3, Hz], '
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
            -uptake * self._vol_dv(t) / self.const.vol_rbc,
        ]
        return dydt

    def _reference_init_to_true_m(self, dv_init: List[float], t0: float = 0.0) -> List[float]:
        """Map 1 fL-reference molarities to M at V_DV(t0), preserving fg inventory."""
        v0 = self._vol_dv(t0)
        scale = self.const.vol_dv / v0
        return [c * scale for c in dv_init]

    def _concentrations_to_fgcell(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert ODE M to fg Fe/cell using V_DV(t) for DV species."""
        t_min = (df.index.to_numpy(dtype=float) - 16.0) * 60.0
        vols = [garnie_dd2_vol_dv_L(t)[0] for t in t_min]
        vol_series = pd.Series(vols, index=df.index)
        out = pd.DataFrame(index=df.index)
        for col in df.columns:
            if col == self.HOST_KEY:
                out[col] = df[col] * self.const.vol_rbc * (10 ** 15) * 55.85
            else:
                out[col] = df[col] * vol_series * (10 ** 15) * 55.85
        return out

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]
        if len(init) == len(self.DV_SPECIES):
            host = self._initial_host_rbc_m(init)
            y0 = self._reference_init_to_true_m(init) + [host]
        else:
            dv = list(init[:len(self.DV_SPECIES)])
            host = float(init[len(self.DV_SPECIES)])
            y0 = self._reference_init_to_true_m(dv) + [host]

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
