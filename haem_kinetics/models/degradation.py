import math
import pandas as pd

from scipy.integrate import solve_ivp
from typing import List, Optional

from haem_kinetics.models.base import KineticsModel
from haem_kinetics.components.experimental_data import ExperimentalData


class Degradation(KineticsModel):
    """
    This is the simplest model to simulate haemoglobin catabolism in the malaria parasite.
    In this case, we have assumed O2- is effectively 0 M given the presence of SOD, thus the reduction of
    Fe(III)PP is ignored. While this model correctly predicts an increase in Hz formation over a time period
    that is relevant to the life cycle of the troph, it is unable to account for the basal “free haem”
    levels as measured by Combrink et al. Consequently, an alteration to the model was necessary which is
    described in Model 2.
    """
    def __init__(self, model_name: str = 'Degradation'):
        super().__init__(model_name=model_name)

        # Initialise concentrations
        self._set_initial_conc(init=[0.0, 0.0])

        # Set Experimental data
        self.exp_data = ExperimentalData()
        # self.exp_data.no_drug_nf54()
        self.exp_data.no_drug_dd2()

    def _calc_hb_trans_linear(self):
        """
        conc_hb_rbc (mol/L) * volRBC (L) = mol  (in RBC, i.e per cell)
        mol / volDV (L) = conc_hb_dv (mol/L)  (i.e. max conc. Fe in DV at tmax)

        rate = conc_hb_dv / (hrs * 60)  # M/min
        k = rate / conc_hb_rbc
        :return:
        """
        conc_hb_dv = self.const.conc_hb_rbc * self.const.vol_rbc / self.const.vol_dv  # Max conc Fe in DV
        rate = conc_hb_dv / (44 * 60)
        return rate / self.const.conc_hb_rbc

    def _calc_enzyme_rate(self, enzyme, conc_hb_dv, t):
        """
        Michaelis-Menton equation

        rate = kcat * conc_enzyme / (Km + conc_substrate)

        :param enzyme:
        :param conc_hb_dv:
        :param t:
        :return:
        """
        kcat = self.const.k_enzymes[enzyme]['kcat'] * 60  # Converts s-1 to min-1
        Km = self.const.k_enzymes[enzyme]['Km']
        conc_enzyme = self._fraction_exp_growth(t) * self.const.conc_enzymes[enzyme] * self.const.fudge
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

        deg = 0
        for enzyme in ['plm_1', 'plm_2', 'hap', 'plm_4']:
        # for enzyme in ['hap']:
            deg += self._calc_enzyme_rate(enzyme=enzyme,
                                          conc_hb_dv=conc_hb_dv,
                                          t=t)
        removal = 4 * deg * conc_hb_dv  # 4 * rate of Hb deg because we release 4 haems per Hb
        return removal

    # def _d_hb_dv_linear(self):
    #
    #     form = self._calc_hb_trans_linear() * self.const.conc_hb_rbc
    #     # remove = self._hb_removal()
    #     remove = 0
    #
    #     return form - remove

    # def _d_hb_dv_dirkie(self, t):
    #
    #     # Formation
    #     # form = self.const.k_hb_trans * self.const.conc_hb_rbc
    #     a = 13.1
    #     b = 8.3e-4
    #     form = a * b * (math.e ** (b * t)) / 55.845
    #
    #     # Removal
    #     # remove = self._hb_removal()
    #     remove = 0
    #
    #     return form - remove

    @staticmethod
    def _fraction_exp_growth(t):
        """
        Fractional exponential growth
        :param t:
        :return:
        """
        a = 0.1578    # t0 = 16hrs
        b = 0.001102  # t0 = 16hrs
        # a = 0.06427   # t0 = 0hrs
        # b = 0.001036  # t0 = 0hrs

        return a * b * (math.e ** (b * t))

    @staticmethod
    def _fraction_lin_growth():
        """
        Fractional linear growth
        :param t:
        :return:
        """
        return 0.0003788

    def _d_hb_dv_kuter(self, t):
        """
        abe^(bt)
        :return:
        """
        tot_hb_conc = (self.const.conc_hb_rbc * self.const.vol_rbc / self.const.vol_dv)
        form = self._fraction_exp_growth(t) * tot_hb_conc

        # Removal
        # remove = self._hb_removal(t=t)
        remove = 0

        return form - remove

    # def _d_hb_dv_sigmoid(self, t):
    #     """
    #
    #     :return:
    #     """
    #     ln_t = math.log(t) if t !=0 else t
    #     top = 1.924
    #     bot = 0.5633
    #     k = 7.121
    #     h = 2.493
    #     denom = (1 + (math.e ** (h * (k - ln_t))))**2
    #     form = h * (top - bot) * (math.e ** (h * (k - ln_t))) / denom
    #
    #     # Removal
    #     remove = self._hb_removal(t=t)
    #
    #     return form - remove

    def _d_fe2pp(self, t):

        # Formation
        form = self._hb_removal(t=t)

        return form

    def _set_initial_conc(self, init: List[float]):
        """
        Sets the initial concentrations of haem species to be integrated

        :param init: List of concentrations (in M) - order matters
        :return:
        """

        if len(init) != 2:
            raise ValueError(f'The number of initial values needed for this model')

        self.initial_values['conc_hb_dv'] = init[0]  # Concentration of Haemoglobin in the digestive vacuole
        self.initial_values['conc_fe2pp'] = init[1]  # Concentration of Free Fe(II) haem

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

        # return [self._d_hb_dv_dirkie(t), self._d_fe2pp()]
        return [self._d_hb_dv_kuter(t), self._d_fe2pp(t)]
        # return [self._d_hb_dv_linear(), self._d_fe2pp(t)]

    def run(self, t, init: Optional[List[float]] = None, plot: Optional[str] = None, **kwargs):
        """

        :param t:
        :param init:
        :param plot: [Optional] Name of file to save plot to. If None, no plot is generated.
        :param kwargs:
        :return:
        """
        if init is None:
            init = [0.0, 0.0]

        if kwargs is False:
            kwargs = {}

        # Reset the conc of Hb in RBC based on initial values supplied
        self._set_initial_conc(init)
        tot_init = 0
        for _, v in self.initial_values.items():
            tot_init += v
        self.const.conc_hb_rbc = self.const.conc_hb_rbc - (tot_init * self.const.vol_dv / self.const.vol_rbc)

        # Solve the differential equations
        self.solution = solve_ivp(self._integrate, t, init, method='BDF', **kwargs)
        self.time = 16 + self.solution.t / 60  # In hours, offset by 16 for parasite life-cycle
        # self.time = self.solution.t / 60  # In hours
        self.concentrations = pd.DataFrame(self.solution.y, columns=self.time, index=list(self.initial_values.keys())).T
        self.concentrations = self._molar_to_fgcell(df=self.concentrations)
        self.concentrations['conc_hz'] = 0.0
        self.concentrations['conc_hb_dv_obs'] = self.concentrations['conc_hb_dv'] - self.concentrations['conc_fe2pp']

        # Plot graph
        if plot:
            self._plot(save_file=plot, title=self.model_name, exp_data=self.exp_data,
                       # columns=['conc_hb_dv', 'conc_hz', 'conc_fe2pp'])
                       # columns=['conc_hb_dv', 'conc_fe2pp', 'conc_hb_dv_obs', 'conc_hz'])
                       columns=['conc_hb_dv_obs', 'conc_hz'])
