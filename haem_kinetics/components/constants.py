
class Constants:
    def __init__(self):

        # -------------------------------------------------------------------------------------
        # Miscellaneous Constants
        # -------------------------------------------------------------------------------------
        self.avogadro = 6.022e23
        self.fudge = 10  # Fudge factor

        # - Volumes
        # self.vol_rbc = 90e-15  # Volume of RBC is 90 fL, reported here in L
        self.vol_dv = 1e-15  # Volume of digestive vacuole is 4 fL, reported here in L
        self.vol_fract_lip = 0.016  # Fractional volume of a lipid nanosphere relative to the DV volume

        # - Other
        self.num_prots = 1.9e8  # Average number of proteins in a P.falciparium

        # -------------------------------------------------------------------------------------
        # Concentrations
        # -------------------------------------------------------------------------------------
        # - Concentration of haemoglobin in the red blood cell (RBC)
        self.conc_hb_rbc = 0.0
        # - Concentration of oxygen [O2]
        #   From Prof. Egan: "Based on Hb saturation curve (30% at 3% O2)"
        self.conc_oxy = 1e-3  # Molar
        # - Concentration of superoxide [O2-]
        #   We assume this to be 0 because of superoxide dismutase
        self.conc_supoxy = 0  # Molar
        # - Concentration of enzymes
        #   Values of enzyme concentrations (ppm) were obtained from paxDB
        self.conc_enzymes = {'plm_1': self._dv_ppm_to_molar(ppm=754),
                             'plm_2': self._dv_ppm_to_molar(ppm=585),
                             'hap': self._dv_ppm_to_molar(ppm=1_373),
                             'plm_4': self._dv_ppm_to_molar(ppm=1_377)}

        # -------------------------------------------------------------------------------------
        # Rate Constants
        # -------------------------------------------------------------------------------------
        # - Estimated Rate constant for the transport of haemoglobin in the digestive vacuole
        #   DKuter (2023-03): I think this is an estimation from a Roepe paper?
        self.k_hb_trans = 0.0011  # M.min-1

        # - Observed rate constant for the degradation of haemoglobin
        self.k_hb_deg = 0.0
        # - Rate of Fe(II) haem oxidation by O2
        #   https://pubs.acs.org/doi/pdf/10.1021/bi00878a025
        self.k_fe2pp_ox = 193_800  # min-1
        # - Rate of Fe(III) haem reduction by O2-
        #   This is an unknown value but doesn't matter when we assume [O2-] is 0 (see above)
        self.k_fe3pp_red = 180e-9
        # - Rate of haemozoin formation
        #   https://link.springer.com/article/10.1186/1475-2875-11-337
        self.k_hz = 0.12  # min-1
        # - Enzyme rate constants
        #   1. https://www.sciencedirect.com/science/article/pii/0166685196026515
        #   2. https://www.sciencedirect.com/science/article/abs/pii/S0731708503005661
        #   3. https://pubs.acs.org/doi/abs/10.1021/bi048252q
        self.k_enzymes = {'plm_1': {'kcat': 2.3,     # s-1 (ref 1)
                                    'Km': 0.49e-6},  # (ref 1)
                          'plm_2': {'kcat': 11,      # s-1 (ref 2)
                                    'Km': 2.6e-6},   # (ref 2)
                          'hap': {'kcat': 0.1,       # s-1 (ref 1)
                                  'Km': 0},         # Km is unknown and is set to 0
                          'plm_4': {'kcat': 1.05,    # s-1
                                    'Km': 0.33e-6}}

        # -------------------------------------------------------------------------------------
        # Equilibrium constants
        # -------------------------------------------------------------------------------------
        self.K_partition = 398  # Fe(III)PPIX lipid partitioning coefficient



    def _dv_ppm_to_molar(self, ppm) -> float:
        """
        Converts the concentration of an enzyme in the digestive vacuole from ppm to Molar.
        Formula:

        mol = ppm * (number of proteins in cell) / (avogadro's constant)
        Molar = mol / (volume of the digestive vacuole)

        :return: Enzyme concentration in Molar
        """
        # Convert from parts per million to parts per 1
        conc = ppm * 10**-6

        # Convert to mol
        conc = conc * self.num_prots / self.avogadro

        # Convert to Molar
        return conc / self.vol_dv

    def compute_conc_hb_rcb(self):
        """
        Computes the concentration of Haemoglobin (in M) in the red blood cell (RBC)
        """

        # Average haemoglobin concentration in the RBC from:
        # https://medlineplus.gov/ency/article/003648.htm
        hb_conc = 34  # g.dL-1
        hb_conc = hb_conc / 0.1  # There are 0.1 dL in 1 L

        # Convert to Molar by dividing by MW of Hb
        hb_conc = hb_conc / 64_500

        # Convert from per haemoglobin molecule to per haem
        self.conc_hb_rbc = hb_conc * 4 * self.fudge

    @staticmethod
    def _compute_kobs(kcat: float, Km: float, enzyme_conc: float) -> float:
        """
        Computes the observed rate constant from the Michaelis–Menten values using the
        formula below. Note the catalytic rate constant input must be in s-1 but
        the returned observed rate constant is in min-1.

        k_obs = kcat * (enzyme concentration) / (Km * (enzyme concentration))

        :param kcat: Catalytic rate constant (in s-1)
        :param Km: Equilibrium constant
        :param enzyme_conc: Enzyme concentration (in M)
        :return: observed rate constant (in min-1)
        """
        k_obs = kcat * enzyme_conc / (Km + enzyme_conc)  # in s-1

        # Convert to min-1
        return k_obs * 60

    def compute_rate_hb_deg_old(self):
        """
        This is the method used to compute a constant kobs for haemoglobin degradation. I suspect there is a mistake
        in the derivation of the _compute_kobs which Prof. and I made. We use the conc of the enzyme instead of the
        conc of the substrate (Hb_dv) - according to Michaelis-Menten, this is wrong.
        :return:
        """
        # Convert from Michaelis–Menten to k observed
        k_obs_plm_1 = self._compute_kobs(kcat=self.k_enzymes['plm_1']['kcat'],
                                         Km=self.k_enzymes['plm_1']['Km'],
                                         enzyme_conc=self.conc_enzymes['plm_1'])
        k_obs_plm_2 = self._compute_kobs(kcat=self.k_enzymes['plm_2']['kcat'],
                                         Km=self.k_enzymes['plm_2']['Km'],
                                         enzyme_conc=self.conc_enzymes['plm_2'])
        k_obs_hap = self._compute_kobs(kcat=self.k_enzymes['hap']['kcat'],
                                         Km=self.k_enzymes['hap']['Km'],
                                         enzyme_conc=self.conc_enzymes['hap'])
        k_obs_plm_4 = self._compute_kobs(kcat=self.k_enzymes['hap']['kcat'],
                                         Km=self.k_enzymes['hap']['Km'],
                                         enzyme_conc=self.conc_enzymes['hap'])

        # Compute the overall observed rate constant for haemoglobin degradation
        # self.k_hb_deg = k_obs_plm_1*conc_plm_1 + k_obs_plm_2*conc_plm_2 + k_obs_hap*conc_hap
        self.k_hb_deg = k_obs_plm_1 * self.conc_enzymes['plm_1'] + k_obs_plm_2 * self.conc_enzymes['plm_2'] + \
                        k_obs_hap * self.conc_enzymes['hap'] + k_obs_plm_4 * self.conc_enzymes['plm_4']

    def compute_lipid_seq_constant(self):
        # return 1 / (1 + (self.vol_fract_lip * self.K_partition))
        return (1 - self.vol_fract_lip) / (1 + self.vol_fract_lip + (self.vol_fract_lip * self.K_partition))

    def compute_values(self):

        # Compute Haemoglobin concentration in the red blood cell
        self.compute_conc_hb_rcb()

        # Compute rate constant for haemoglobin degradation
        self.compute_rate_hb_deg_old()
