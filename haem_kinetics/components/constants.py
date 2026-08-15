import math


class Constants:
    def __init__(self):

        # -------------------------------------------------------------------------------------
        # Miscellaneous Constants
        # -------------------------------------------------------------------------------------
        self.avogadro = 6.022e23

        # - Volumes
        self.vol_rbc = 90e-15  # Volume of RBC is 90 fL, reported here in L
        # Reference DV volume (1 fL) for the init API and PaxDB amount
        # n_E = [E]_1fL · 1 fL. Active models convert to molarity with variable_dv_volume.
        self.vol_dv = 1e-15  # L (1 fL) — reference, not the integrated lumen volume
        self.vol_fract_lip = 0.016  # Fractional volume of a lipid nanosphere relative to the DV volume

        # - Other
        self.num_prots = 1.9e8  # Average number of proteins in a P.falciparum
        # Total Fe budget per infected RBC (~106 fg) for mass-balance diagnostics
        self.total_fe_fg_cell = self.compute_conc_hb_rcb() * 90.0 * 55.85

        # -------------------------------------------------------------------------------------
        # Concentrations
        # -------------------------------------------------------------------------------------
        # - Concentration of haemoglobin in the red blood cell (RBC)
        self.conc_hb_rbc = self.compute_conc_hb_rcb()
        # - Concentration of oxygen [O2]
        #   From Prof. Egan: "Based on Hb saturation curve (30% at 3% O2)"
        self.conc_oxy = 1e-3  # Molar
        # - Concentration of superoxide [O2-]
        #   We assume this to be 0 because of superoxide dismutase
        self.conc_supoxy = 0  # Molar
        # - Enzyme abundances (ppm) -> DV M via _dv_ppm_to_molar
        #   ([E] = f(ppm, N_prot, V_DV); not an independent constant).
        #   Source: PaxDB "P.falciparum 3D7 - Whole organism, Dd2, SC (Tao,MCP,2014)"
        #   https://pax-db.org/ — gene IDs in comments below.
        self.ppm_enzymes = {
            'plm_1': 752,     # PF3D7_1407900 / Q7KQM4
            'plm_2': 1204,    # PF3D7_1408000 / Q8I6V3
            'hap': 1373,      # PF3D7_1408100 / Q8IM15 (HAP; PaxDB: "PM III")
            'plm_4': 3139,    # PF3D7_1407800 / Q8IM16
            'fp_2': 20.2,     # PF3D7_1115700 / Q8I6U4 (falcipain-2a)
            'fp_3': 23.5,     # PF3D7_1115400 / Q8IIL0
        }
        self.conc_enzymes = {
            name: self._dv_ppm_to_molar(ppm) for name, ppm in self.ppm_enzymes.items()
        }

        # -------------------------------------------------------------------------------------
        # Rate Constants
        # -------------------------------------------------------------------------------------
        # This is the value needed to obtain ~100 fg/cell Hb in the DV at t~45hrs
        # self.k_hb_trans = 0.011 / self.conc_hb_rbc  # 0.011 M.min-1 / haemoglobin conc (M)
        self.k_hb_trans = 0.000007986 / self.conc_hb_rbc
        # - Observed rate constant for the degradation of haemoglobin
        self.k_hb_deg = 0.0
        # - Rate of Fe(II) haem oxidation by O2
        #   https://pubs.acs.org/doi/pdf/10.1021/bi00878a025
        self.k_fe2pp_ox = 193_800  # min-1
        # - Rate of Fe(III) haem reduction by O2-
        #   This is an unknown value but doesn't matter when we assume [O2-] is 0 (see above)
        self.k_fe3pp_red = 180e-9
        # - Rate of haemozoin formation (lipid-mediated β-haematin)
        #   https://link.springer.com/article/10.1186/1475-2875-11-337
        #   Active Model 6: apply to lipid-associated Fe(III), not bulk aqueous haem.
        #   Do not multiply by φ. Crystal-competent (xtal) is a later step if needed.
        self.k_hz = 0.12  # min-1
        # - Aqueous <-> lipid Fe(III) exchange; large => near-equilibrium partition
        self.k_lipid_exchange = 50.0  # min-1
        # - Lipid-associated <-> crystal-competent Fe(III) exchange
        #   Not used in active Model 6 (aq ⇄ lip only). Kept for a later numbered
        #   model / legacy Model 6 if assay Hm still drains.
        self.k_xtal_exchange = 5.0  # min-1
        self.K_xtal = 0.08  # [Fe3]_xtal / [Fe3]_lip at equilibrium
        # - Enzyme rate constants (peptide-substrate MM; see docs/enzyme_kinetics.md)
        #   Plasmepsins / HAP: Banerjee et al. PNAS 2002 Table 1 (PM I/II from Luker 1996)
        #     https://doi.org/10.1073/pnas.022630099
        #     https://doi.org/10.1016/0166-6851(96)02651-5
        #   Falcipains: Ramjee et al. Biochem J 2006 (FRET peptides)
        #     https://doi.org/10.1042/BJ20060422
        self.k_enzymes = {
            'plm_1': {'kcat': 2.3, 'Km': 0.49e-6},   # s-1, M — Luker/Banerjee α33–34
            'plm_2': {'kcat': 11, 'Km': 2.6e-6},     # s-1, M — Luker/Banerjee α33–34
            'hap': {'kcat': 0.05, 'Km': 0.29e-6},    # s-1, M — Banerjee native HAP
            'plm_4': {'kcat': 1.05, 'Km': 0.33e-6},  # s-1, M — Banerjee recombinant PM IV
            'fp_2': {'kcat': 0.79, 'Km': 0.9e-6},    # s-1, M — Ramjee FP-2 best FRET
            'fp_3': {'kcat': 0.204, 'Km': 4.0e-6},   # s-1, M — Ramjee FP-3 Leu-Arg FRET
        }
        # Native tetramer (Models 4a/4b). No paper reports kcat on native Hb in s-1.
        # Gluzman et al. JCI 1994: equal globin-degrading units, PM I is considerably
        # more active on native Hb than PM II — using PM II's peptide kcat (11 s-1)
        # would invert that ranking. PM I peptide kcat is used as a provisional
        # native kcat for PM I; PM II native kcat is set equal to that value (not
        # 11). FP-2 keeps its peptide kcat (Shenai 2000: FP-2 does cut native Hb;
        # Vmax is already ~4 fg/h). Not fit to Garnie standing Hb.
        # https://doi.org/10.1172/jci117140
        self.k_enzymes_native = {
            'plm_1': {'kcat': 2.3, 'Km': 0.49e-6},
            'plm_2': {'kcat': 2.3, 'Km': 2.6e-6},
            'fp_2': {'kcat': 0.79, 'Km': 0.9e-6},
        }
        # HTV / inner-vesicle cargo → lumen (Model 5). First-order; lumps traffic,
        # outer fusion, and inner-membrane lysis. Klemba JCB 2004: if PM II is
        # trafficked exclusively through the cytostome, delivery t½ is < 20 min
        # (PM biosynthesis/maturation t½ ≈ 20 min; Francis 1997; Banerjee 2003).
        # Provisional: use that 20 min bound, not a Garnie standing-Hb fit.
        # https://doi.org/10.1083/jcb.200307147
        self.k_htv_release = math.log(2) / 20.0  # min-1

        # -------------------------------------------------------------------------------------
        # Equilibrium constants
        # -------------------------------------------------------------------------------------
        self.K_partition = 398  # Fe(III)PPIX lipid partitioning coefficient

    def _dv_ppm_to_molar(self, ppm, vol_dv=None) -> float:
        """
        Converts the concentration of an enzyme in the digestive vacuole from ppm to Molar.
        Formula:

        mol = ppm * (number of proteins in cell) / (avogadro's constant)
        Molar = mol / (volume of the digestive vacuole)

        Enzyme *amount* is independent of V_DV; molarity scales as 1/V_DV.
        :return: Enzyme concentration in Molar
        """
        vol = self.vol_dv if vol_dv is None else vol_dv
        if vol <= 0.0:
            raise ValueError(f'vol_dv must be positive, got {vol}')
        # Convert from parts per million to parts per 1
        conc = ppm * 10**-6

        # Convert to mol
        conc = conc * self.num_prots / self.avogadro

        # Convert to Molar
        return conc / vol

    @staticmethod
    def compute_conc_hb_rcb() -> float:
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
        return hb_conc * 4

    def compute_lipid_seq_constant(self):
        return (1 - self.vol_fract_lip) / (1 + self.vol_fract_lip + (self.vol_fract_lip * self.K_partition))
