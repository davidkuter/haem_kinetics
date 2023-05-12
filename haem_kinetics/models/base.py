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
        self.const.compute_values()

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
              exp_data: Optional[ExperimentalData] = None):
        # Set which data will be plotted
        df_plot = self.concentrations if columns is None else self.concentrations[columns]

        # Set up plot
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
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
        axes[1].plot(df_plot.index, df_plot['conc_hz'], 'r', label='conc_hz')
        if exp_data:
            axes[1].errorbar(exp_data.data.index, exp_data.data['Hz'], yerr=exp_data.data['Hz:SEM'].values,
                             label='Exp Hz', fmt="o", mfc='white', ecolor='r', color='r')
        axes[1].legend(loc='upper right')

        # Plot remaining Haem species
        cols = [col for col in df_plot.columns if col != 'conc_hz']
        axes[0].plot(df_plot.index, df_plot[cols], label=cols)
        if exp_data:
            axes[0].errorbar(exp_data.data.index, exp_data.data['Hm'], yerr=exp_data.data['Hm:SEM'].values,
                             label='Exp Haem', fmt="o", mfc='white', ecolor='orange', color='orange')
            axes[0].errorbar(exp_data.data.index, exp_data.data['Hb'], yerr=exp_data.data['Hb:SEM'].values,
                             label='Exp Hb', fmt="o", mfc='white', ecolor='b', color='b')
        axes[0].legend(loc='upper left')

        # Format

        # Save output
        axes[1].legend(loc='upper left')
        plt.savefig(save_file)

    def _molar_to_fgcell(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Converts mol/L (molar) concentrations to fg/cell.

        mol -> fg: mol * 10^15 / mw of Fe = fg/L
        fg/L -> fg/cell: fg/cell * volume of cell


        :param df: Dataframe of concentrations in molar
        """

        return df * self.const.vol_rbc * (10 ** 15) / 55.85

    def run(self, t, init, kwargs):
        """
        API that solves differential equations and saves output.

        :param t:
        :param init:
        :param kwargs:
        :return:
        """

        raise NotImplementedError('"run" must be overwritten by the model class')
