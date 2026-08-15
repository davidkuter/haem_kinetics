import matplotlib.pyplot as plt
import pandas as pd

from scipy.integrate import solve_ivp
from typing import List, Optional

from haem_kinetics.components.constants import Constants
from haem_kinetics.components.experimental_data import ExperimentalData
from haem_kinetics.components.fit_metrics import (
    SKIP_TOTAL,
    format_fit_metrics,
    score_fractionation,
)
from haem_kinetics.models.helpers import variable_dv_volume_L


class KineticsModel:
    """Base class for haem speciation ODEs.

    Host RBC haemoglobin (`conc_hb_rbc`) is part of the ODE state in all models so
    uptake depletes the finite Fe budget (~106 fg/cell). Callers may still pass
    DV-only initial conditions; `run()` appends the remaining host concentration.

    Active models use variable_dv_volume for molar bookkeeping (dilution,
    [E] = n_E/V(t), fg = C·V(t)). `constants.vol_dv` (1 fL) is the reference
    volume for the init API and PaxDB amount. Legacy models set
    variable_dv_volume = False and keep a fixed 1 fL lumen.
    """

    HOST_KEY = 'conc_hb_rbc'
    variable_dv_volume = True
    DV_SPECIES: List[str] = []
    # Species stored as moles/cell encoded as M at V_ref (1 fL). Not diluted by
    # dV_lumen/dt and converted to fg with V_ref, not Garnie lumen volume.
    AMOUNT_SPECIES: List[str] = []

    def __init__(self, model_name):

        # General
        self.model_name = model_name

        # Grab the constants required to integrate equations
        self.const = Constants()

        # To be solved
        self.initial_values = {}    # Stores initial concentrations of haem species using in integration
        self.differential_eqs = []  # Stores the differential eqs to be integrated

        # Solved results
        self.concentrations = pd.DataFrame()    # Stores the solved time-course concentrations of haem species
        self.time = []              # Stores the time-series for the solution
        self.solution = None        # Stores the entire integrated solution (gives access to additional info if needed)
        self.fit_metrics = None     # Dd2 fractionation scores after score_vs_experiment()

    def _set_initial_conc(self, init):
        """
        Function to set the initial concentration of haem species. Since each model could have different haem species
        or ordering of haem species, this class function must be overwritten for each model
        :param init: Initial values for haem concentrations (Order matters!)
        :return:
        """
        raise NotImplementedError('_set_initial_conc must be overwritten by the model class')

    def _nonneg(self, x: float) -> float:
        """Physical domain helper: rate laws are undefined for negative concentrations."""
        x = float(x)
        return x if x > 0.0 else 0.0

    def _solve_ivp(self, fun, t_span, y0, **kwargs):
        """Integrate the stated ODEs (BDF). Tight tol = accurate solve, not a model change."""
        opts = {'method': 'BDF', 'rtol': 1e-8, 'atol': 1e-12, **kwargs}
        return solve_ivp(fun, t_span, y0, **opts)

    def _full_hb_rbc_m(self) -> float:
        """Uninfected-RBC haem-equivalent concentration (M), before any parasite uptake."""
        return Constants.compute_conc_hb_rcb()

    def _vol_dv(self, t: float) -> float:
        """Current DV lumen volume (L). Variable V_DV(t) unless a subclass opts out."""
        if not self.variable_dv_volume:
            return self.const.vol_dv
        vol, _dvol = variable_dv_volume_L(t)
        if vol <= 0.0:
            raise ValueError(
                f'DV lumen volume is {vol} L at t={t} min; '
                'concentration ODEs are undefined at V≤0 (collapse reaches 0 at 46 h)'
            )
        return vol

    def _dvol_dv(self, t: float) -> float:
        """dV_DV/dt (L/min)."""
        if not self.variable_dv_volume:
            return 0.0
        _vol, dvol = variable_dv_volume_L(t)
        return dvol

    def _dilution(self, t: float, conc: float) -> float:
        """Chain rule for C = n/V: dC/dt includes −C·(dV/dt)/V."""
        vol = self._vol_dv(t)
        return -conc * self._dvol_dv(t) / vol

    def _enzyme_conc(self, enzyme: str, t: float) -> float:
        """PaxDB amount in the current lumen: [E] = n_E / V_DV(t)."""
        ppm = self.const.ppm_enzymes[enzyme]
        return self.const._dv_ppm_to_molar(ppm, vol_dv=self._vol_dv(t))

    def _initial_host_rbc_m(self, dv_init: List[float]) -> float:
        """
        Remaining host [Hb]_RBC (M) after accounting for Fe already in DV species.

        `dv_init` is in 1 fL-reference molarities so the fg seed is
        tot_dv · vol_dv_ref; convert that amount to an RBC-basis concentration.
        """
        tot_dv = float(sum(dv_init))
        host = self._full_hb_rbc_m() - tot_dv * self.const.vol_dv / self.const.vol_rbc
        return max(host, 0.0)

    def _pad_init_with_host(self, dv_init: List[float]) -> List[float]:
        """Append remaining host [Hb]_RBC to a DV-only initial state vector."""
        return list(dv_init) + [self._initial_host_rbc_m(dv_init)]

    def _reference_init_to_true_m(self, dv_init: List[float], t0: float = 0.0) -> List[float]:
        """Map 1 fL-reference molarities to M at V_DV(t0), preserving fg inventory.

        Amount-encoded species stay at V_ref molarity (they are not in the lumen).
        """
        v0 = self._vol_dv(t0)
        scale = self.const.vol_dv / v0
        amount = set(self.AMOUNT_SPECIES)
        out = []
        for name, c in zip(self.DV_SPECIES, dv_init):
            out.append(c if name in amount else c * scale)
        return out

    def _prepare_y0(self, init: List[float], t0: float = 0.0) -> List[float]:
        """Build the integrator y0: true-M DV species plus host (RBC basis)."""
        n_dv = len(self.DV_SPECIES)
        if len(init) == n_dv:
            host = self._initial_host_rbc_m(init)
            return self._reference_init_to_true_m(init, t0) + [host]
        dv = list(init[:n_dv])
        host = float(init[n_dv])
        return self._reference_init_to_true_m(dv, t0) + [host]

    def _d_host_from_dv_uptake(self, uptake_dv_m_per_min: float, t: Optional[float] = None) -> float:
        """d[Hb]_RBC/dt given Hb appearance rate in the DV (M/min on current V_DV).

        Legacy callers omit `t` and use the 1 fL reference volume.
        """
        vol = self.const.vol_dv if t is None else self._vol_dv(t)
        return -uptake_dv_m_per_min * vol / self.const.vol_rbc

    def _concentrations_to_fgcell(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert ODE outputs to fg Fe/cell.

        Lumen DV species use V_DV(t) from the age index (16 + t/60 h); host RBC
        Hb uses vol_rbc; amount-encoded species (AMOUNT_SPECIES) use V_ref.
        No clipping: the fg time course must reflect the ODE solution.
        """
        t_min = (df.index.to_numpy(dtype=float) - 16.0) * 60.0
        vols = [self._vol_dv(float(t)) for t in t_min]
        vol_series = pd.Series(vols, index=df.index)
        factor = (10 ** 15) * 55.85
        amount = set(self.AMOUNT_SPECIES)
        out = pd.DataFrame(index=df.index)
        for col in df.columns:
            if col == self.HOST_KEY:
                out[col] = df[col] * self.const.vol_rbc * factor
            elif col in amount:
                out[col] = df[col] * self.const.vol_dv * factor
            else:
                out[col] = df[col] * vol_series * factor
        return out

    def _plot(self, save_file: str, title: str, columns: Optional[List[str]] = None,
              exp_data: Optional[ExperimentalData] = None,
              free_haem_cols: Optional[List[str]] = None,
              plot_total_fe: bool = False):
        # Plot the converted solution as-is (no clipping). Negatives mean the
        # integrator failed the stated ODEs — do not hide that in the figure.
        df_plot = self.concentrations.copy()
        if free_haem_cols:
            df_plot['conc_fe3pp_free'] = df_plot[free_haem_cols].sum(axis=1)
        if plot_total_fe:
            skip = {self.HOST_KEY, *SKIP_TOTAL}
            species = [c for c in df_plot.columns
                       if c.startswith('conc_') and c not in skip]
            if self.HOST_KEY in df_plot.columns:
                df_plot['conc_fe_total'] = df_plot[species].sum(axis=1) + df_plot[self.HOST_KEY]
            else:
                df_plot['conc_fe_total'] = df_plot[species].sum(axis=1)

        if columns is None:
            plot_df = df_plot
        else:
            plot_df = df_plot[columns]

        # Set up plot
        n_axes = 3 if plot_total_fe else 2
        fig, axes = plt.subplots(1, n_axes, figsize=(10 * n_axes, 10))
        if n_axes == 2:
            axes = list(axes)
        font_size = 16
        plt.rcParams.update({'font.size': font_size})
        fig.suptitle(title, fontsize=font_size + 8)
        x_range = range(int(self.time[-1]), int(self.time[0]), -5)
        for ax in axes:
            ax.xaxis.label.set_fontsize(font_size)
            ax.yaxis.label.set_fontsize(font_size)
            ax.tick_params(axis='x', labelsize=font_size)
            ax.tick_params(axis='y', labelsize=font_size)

        plt.setp(axes, xticks=x_range, xlabel='Time (hrs)', ylabel='Fe (fg/cell)')

        # Plot Hz
        axes[1].plot(plot_df.index, plot_df['conc_hz'], 'r', label='conc_hz')
        if exp_data is not None and not exp_data.data.empty:
            axes[1].errorbar(exp_data.data.index, exp_data.data['Hz'], yerr=exp_data.data['Hz:SEM'].values,
                             label='Exp Hz', fmt="o", mfc='white', ecolor='r', color='r')
        axes[1].legend(loc='upper left')

        # Plot remaining Haem species
        cols = [col for col in plot_df.columns if col not in ('conc_hz', 'conc_fe_total')]
        axes[0].plot(plot_df.index, plot_df[cols], label=cols)
        if exp_data is not None and not exp_data.data.empty:
            axes[0].errorbar(exp_data.data.index, exp_data.data['Hm'], yerr=exp_data.data['Hm:SEM'].values,
                             label='Exp Haem', fmt="o", mfc='white', ecolor='orange', color='orange')
            axes[0].errorbar(exp_data.data.index, exp_data.data['Hb'], yerr=exp_data.data['Hb:SEM'].values,
                             label='Exp Hb', fmt="o", mfc='white', ecolor='b', color='b')
        axes[0].legend(loc='upper left')

        if plot_total_fe and 'conc_fe_total' in df_plot.columns:
            budget = self.const.total_fe_fg_cell
            axes[2].plot(df_plot.index, df_plot['conc_fe_total'], 'k', label='Total Fe (model)')
            axes[2].axhline(budget, color='gray', linestyle='--',
                            label=f'Budget ({budget:.0f} fg)')
            # Avoid matplotlib offset notation (e.g. +1.06e2 with ±0.015 ticks),
            # which makes conserved ~106 fg look like values near zero.
            axes[2].ticklabel_format(axis='y', style='plain', useOffset=False)
            ymin = min(0.0, float(df_plot['conc_fe_total'].min()) * 0.95)
            ymax = max(budget * 1.15, float(df_plot['conc_fe_total'].max()) * 1.05)
            axes[2].set_ylim(ymin, ymax)
            axes[2].legend(loc='upper left')

        plt.savefig(save_file)

    def score_vs_experiment(
        self,
        exp_data: Optional[ExperimentalData] = None,
        free_haem_cols: Optional[List[str]] = None,
    ):
        """Score fg trajectories against heme fractionation ages (see fit_metrics.py)."""
        data = exp_data if exp_data is not None else getattr(self, 'exp_data', None)
        self.fit_metrics = score_fractionation(
            self.concentrations, data, free_haem_cols=free_haem_cols
        )
        return self.fit_metrics

    def format_fit_metrics(self) -> str:
        if not self.fit_metrics:
            return ''
        return format_fit_metrics(self.fit_metrics)

    def _molar_to_fgcell(self, df: pd.DataFrame) -> pd.DataFrame:
        """DV M → fg/cell using V_DV(t) from the age index (host column excluded)."""
        return self._concentrations_to_fgcell(df)

    def run(self, t, init, kwargs):
        """
        API that solves differential equations and saves output.

        :param t:
        :param init:
        :param kwargs:
        :return:
        """

        raise NotImplementedError('"run" must be overwritten by the model class')
