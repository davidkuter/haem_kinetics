from typing import List, Optional

from haem_kinetics.models.base import KineticsModel


class Model1(KineticsModel):
    """
    This is the simplest model to simulate haemoglobin catabolism in the malaria parasite.
    In this case, we have assumed O2- is effectively 0 M given the presence of SOD, thus the reduction of
    Fe(III)PP is ignored. While this model correctly predicts an increase in Hz formation over a time period
    that is relevant to the life cycle of the troph, it is unable to account for the basal “free haem”
    levels as measured by Combrink et al. Consequently, an alteration to the model was necessary which is
    described in Model 2.
    """
    def __int__(self, model_name='Model 1'):
        super().__init__(self, model_name)

        # Initialise concentrations
        self._set_initial_conc(init=[0.0, 0.0, 0.0, 0.0])

    def _d_hb_dv(self):

        # Formation
        form = self.const.k_hb_trans * self.const.conc_hb_rbc

        # Removal
        remove = self.const.k_hb_deg * self.initial_values['conc_hb_dv']

        return form - remove

    def _d_fe2pp(self):

        # Formation
        form = (self.const.k_hb_deg * self.initial_values['conc_hb_dv']) + \
               (self.const.k_fe3pp_red * self.initial_values['conc_fe3pp'] * self.const.conc_supoxy)

        # Removal
        remove = self.const.k_fe2pp_ox * self.initial_values['conc_fe2pp'] * self.const.conc_oxy

        return form - remove

    def _d_fe3pp(self):

        # Formation
        form = self.const.k_fe2pp_ox * self.initial_values['conc_fe2pp'] * self.const.conc_oxy

        # Removal
        remove = (self.const.k_fe3pp_red * self.initial_values['conc_fe3pp'] * self.const.conc_supoxy) +\
                 (self.const.k_hz * self.initial_values['conc_fe3pp'])

        return form - remove

    def _d_hz(self):

        # Formation
        return self.const.k_hz * self.initial_values['conc_fe3pp']

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

    def _set_diff_eqs(self):
        """
        Adds the differential equations to the list that will be integrated.
        :return:
        """
        self.differential_eqs.extend = [self._d_hb_dv(), self._d_fe2pp(), self._d_fe3pp(), self._d_hz()]

    def run(self, t, init: Optional[List[float]] = None, **kwargs):
        """

        :param t:
        :param init:
        :param kwargs:
        :return:
        """
        if init is None:
            init = [0.0, 0.0, 0.0, 0.0]

        if kwargs is False:
            kwargs = {}

        # Solve the differential equations
        self._solve(t, init, **kwargs)

        # Plot graph
        self._plot()
