import matplotlib.pyplot as plt
import pandas as pd

from typing import List, Optional

from haem_kinetics.components.constants import Constants
from haem_kinetics.components.experimental_data import ExperimentalData


class KineticsModel:
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

    def _set_initial_conc(self, init):
        """
        Function to set the initial concentration of haem species. Since each model could have different haem species
        or ordering of haem species, this class function must be overwritten for each model
        :param init: Initial values for haem concentrations (Order matters!)
        :return:
        """
        raise NotImplementedError('_set_initial_conc must be overwritten by the model class')

    # ToDo: Implement later
    def _plot(self, save_file: str, title: str, columns: Optional[List[str]] = None,
              exp_data: Optional[ExperimentalData] = None,
              free_haem_cols: Optional[List[str]] = None,
              plot_total_fe: bool = False):
        # Set which data will be plotted
        df_plot = self.concentrations.copy()
        if free_haem_cols:
            df_plot['conc_fe3pp_free'] = df_plot[free_haem_cols].sum(axis=1)
        if plot_total_fe:
            fe_cols = [c for c in df_plot.columns
                       if c.startswith('conc_') and c not in ('conc_hb_host',)]
            # Prefer explicit DV species if host Fe is tracked separately
            species = [c for c in fe_cols if c != 'conc_fe3pp_free']
            if 'conc_hb_host' in df_plot.columns:
                df_plot['conc_fe_total'] = df_plot[species].sum(axis=1) + df_plot['conc_hb_host']
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
            axes[2].plot(df_plot.index, df_plot['conc_fe_total'], 'k', label='Total Fe (model)')
            axes[2].axhline(self.const.total_fe_fg_cell, color='gray', linestyle='--',
                            label=f'Budget ({self.const.total_fe_fg_cell:.0f} fg)')
            axes[2].legend(loc='upper left')

        plt.savefig(save_file)

    def _molar_to_fgcell(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Converts mol/L (molar) concentrations to fg/cell. We only need to consider mol/L -> fg because we start off
        the calculation on a per cell basis. When we set the rate of Hb transport, we use 106 fg Fe which is the max
        mass of iron PER CELL. This the per cell actually inherent through the entire calculation.

        I.e. concentration M is actually M/cell, etc.

        mol/cell -> fg/cell: mol * 10^15 * mw of Fe

        :param df: Dataframe of concentrations in molar
        """
        return df * self.const.vol_dv * (10 ** 15) * 55.85

    def run(self, t, init, kwargs):
        """
        API that solves differential equations and saves output.

        :param t:
        :param init:
        :param kwargs:
        :return:
        """

        raise NotImplementedError('"run" must be overwritten by the model class')
