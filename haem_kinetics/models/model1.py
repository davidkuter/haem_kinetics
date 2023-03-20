from haem_kinetics.components.constants import Constants


class Model1:
    """
    This is the simplest model to simulate haemoglobin catabolism in the malaria parasite.
    In this case, we have assumed O2- is effectively 0 M given the presence of SOD, thus the reduction of
    Fe(III)PP is ignored. While this model correctly predicts an increase in Hz formation over a time period
    that is relevant to the life cycle of the troph, it is unable to account for the basal “free haem”
    levels as measured by Combrink et al. Consequently, an alteration to the model was necessary which is
    described in Model 2.
    """
    def __int__(self):
        # Grab the constants required to integrate equations
        self.const = Constants().compute_values()

        # Initialise variables to store output
        self.conc_hb_dv = 0.0  # Concentration of Haemoglobin in the digestive vacuole
        self.conc_fe2pp = 0.0  # Concentration of Free Fe(II) haem
        self.conc_fe3pp = 0.0  # Concentration of Free Fe(III) haem
        self.conc_hz = 0.0     # Concentratin of haemozoin

    def _d_hb_dv(self):

        # Formation
        form = self.const.k_hb_trans * self.const.conc_hb_rbc

        # Removal
        remove = self.const.k_hb_deg * self.conc_hb_dv

        return form - remove

    def _d_fe2pp(self):

        # Formation
        form = (self.const.k_hb_deg * self.conc_hb_dv) + \
               (self.const.k_fe3pp_red * self.conc_fe3pp * self.const.conc_supoxy)

        # Removal
        remove = self.const.k_fe2pp_ox * self.conc_fe2pp * self.const.conc_oxy

        return form - remove

    def _d_fe3pp(self):

        # Formation
        form = self.const.k_fe2pp_ox * self.conc_fe2pp * self.const.conc_oxy

        # Removal
        remove = (self.const.k_fe3pp_red * self.conc_fe3pp * self.const.conc_supoxy) +\
                 (self.const.k_hz * self.conc_fe3pp)

        return form - remove

    def _d_hz(self):

        # Formation
        return self.const.k_hz * self.conc_fe3pp

    def compute(self):
        """

        :return:
        """
        self.conc_hb_dv = self._d_hb_dv()
        self.conc_fe2pp = self._d_fe2pp()
        self.conc_fe3pp = self._d_fe3pp()
        self.conc_hz = self._d_hz()
