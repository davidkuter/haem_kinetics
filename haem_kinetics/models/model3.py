import math
import pandas as pd

from scipy.integrate import solve_ivp
from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.components.experimental_data import ExperimentalData


class Model3(KineticsModel):
    """
    This is the simplest model to simulate haemoglobin catabolism in the malaria parasite.
    In this case, we have assumed:
     * Exponential transport of Hb into the DV
     * Hb enzymatic degradation is rate-limited by HAP
     * Exponential increase in enzyme concentration that follows Hb transport rate
     * There is fudge factor to lower enzyme concentration or kcat.
     * O2- is effectively 0 M given the presence of SOD, thus the reduction of Fe(III)PP is ignored.
     * A portion of Fe3PPIX is sequestered in a lipid droplet
    """
    def __init__(self, model_name: str = 'Model 3'):
        super().__init__(model_name=model_name)

        # Initialise concentrations
        self._set_initial_conc(init=[0.005, 0.0, 0.0, 0.0])

        # Set Experimental data
        self.exp_data = ExperimentalData()
        # self.exp_data.no_drug_nf54()
        self.exp_data.no_drug_dd2()

    @staticmethod
    def _fraction_exp_growth(t):
        """
        Fractional exponential growth
        :param t:
        :return:
        """
        a = 0.1578
        b = 0.001102
        # return a * b * (math.e ** (b * t))
        return a * (math.e ** (b * t))

    def _calc_enzyme_rate(self, enzyme, conc_hb_dv, t):
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60  # Converts s-1 to min-1
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = self._fraction_exp_growth(t) * self.const.conc_enzymes[enzyme] / self.const.fudge
        denom = Km + conc_hb_dv

        if denom == 0:
            return 0
        else:
            return kcat * conc_enzyme / denom

    def _hb_removal(self, t):
        """
        Haemoglobin degradation by enzymes
        :return:
        """
        conc_hb_dv = self.initial_values['conc_hb_dv'] / 4
        # plm_1_deg = self._calc_enzyme_rate(enzyme='plm_1',
        #                                    conc_hb_dv=conc_hb_dv)
        # plm_2_deg = self._calc_enzyme_rate(enzyme='plm_2',
        #                                    conc_hb_dv=conc_hb_dv)
        hap_deg = self._calc_enzyme_rate(enzyme='hap',
                                         conc_hb_dv=conc_hb_dv, t=t)
        # plm_4_deg = self._calc_enzyme_rate(enzyme='plm_4',
        #                                    conc_hb_dv=conc_hb_dv)
        conc_hb_dv = self.initial_values['conc_hb_dv']

        # removal = (plm_1_deg + plm_2_deg + hap_deg + plm_4_deg) * conc_hb_dv
        # removal = (plm_1_deg + plm_2_deg) * conc_hb_dv
        removal = hap_deg * conc_hb_dv
        return removal

    def _d_hb_dv(self, t):

        # Formation
        form = self._fraction_exp_growth(t) * self.const.conc_hb_rbc / 40  # Need to investigate this (40)

        # Removal
        remove = self._hb_removal(t)

        return form - remove

    def _d_fe2pp(self, t):

        # Formation
        form = self._hb_removal(t) + \
               (self.const.k_fe3pp_red * self.const.compute_lipid_seq_constant()
                * self.initial_values['conc_fe3pp'] * self.const.conc_supoxy)

        # Removal
        remove = self.const.k_fe2pp_ox * self.initial_values['conc_fe2pp'] * self.const.conc_oxy

        return form - remove

    def _d_fe3pp(self):

        # Formation
        form = self.const.k_fe2pp_ox * self.initial_values['conc_fe2pp'] * self.const.conc_oxy

        # Removal
        remove = (self.const.k_fe3pp_red * self.const.compute_lipid_seq_constant()
                  * self.initial_values['conc_fe3pp'] * self.const.conc_supoxy) + \
                 (self.const.k_hz * self.const.compute_lipid_seq_constant() * self.initial_values['conc_fe3pp'])

        return form - remove

    def _d_hz(self):

        # Formation
        return self.const.k_hz * self.const.compute_lipid_seq_constant() * self.initial_values['conc_fe3pp']

    def _set_initial_conc(self, init: List[float]):
        """
        Sets the initial concentrations of haem species to be integrated

        :param init: List of concentrations (in M) - order matters
        :return:
        """

        if len(init) != 4:
            raise ValueError(f'The number of initial values needed for this model')

        self.initial_values['conc_hb_dv'] = init[0]  # Concentration of Haemoglobin in the digestive vacuole
        self.initial_values['conc_fe2pp'] = init[1]  # Concentration of Free Fe(II) haem
        self.initial_values['conc_fe3pp'] = init[2]  # Concentration of Free Fe(III) haem
        self.initial_values['conc_hz'] = init[3]     # Concentration of haemozoin

    def _integrate(self, t, init):
        """
        Function that will integrate differential equations for the model. The differential equations must be set
        by within the model but adding them to the self.differential_eqs list.

        :param t: Time range that will be integrated over
        :param init: Initial values for haem concentrations (Order matters!)
        :return: returns concentrations for each haem species in the model at each time point. E.g.:
                   t0      t1      t2     t3     t4
                 [
                  [0.    0.001    0.01   0.02   0.03]     # Haem species 1
                  [0.    0.       0.001  0.002  0.003]    # Haem species 2
                 ]
        """
        # Set initial concentration values
        self._set_initial_conc(init=init)

        return [self._d_hb_dv(t), self._d_fe2pp(t), self._d_fe3pp(), self._d_hz()]

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        """

        :param t:
        :param init:
        :param plot: [Optional] Name of file to save plot to. If None, no plot is generated.
        :param kwargs:
        :return:
        """
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]

        if kwargs is False:
            kwargs = {}

        # Solve the differential equations
        self.solution = solve_ivp(self._integrate, t, init, **kwargs)
        self.time = 16 + self.solution.t / 60  # In hours, offset by 16 for parasite life-cycle
        self.concentrations = pd.DataFrame(self.solution.y, columns=self.time, index=list(self.initial_values.keys())).T
        self.concentrations = self.concentrations * 1000 * 0.2232  # convert to fg/cell

        # Plot graph
        if plot:
            self._plot(save_file=plot, title=self.model_name, exp_data=self.exp_data,
                       columns=['conc_hb_dv', 'conc_hz', 'conc_fe3pp'])
