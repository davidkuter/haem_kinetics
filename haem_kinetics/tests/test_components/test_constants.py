from haem_kinetics.components.constants import Constants


def test_compute_conc_hb_rcb():
    """
    Test computing the concentration of haemoglobin in the red blood cell
    Expected value is 80 mM
    :return:
    """

    constants = Constants()
    constants.compute_conc_hb_rcb()
    conc_rbc = round(constants.conc_hb_rbc, 3)

    assert conc_rbc == 0.080


def test_dv_ppm_to_molar():
    """
    Test computing the concentration of haem degrading enzymes.
    Note: Values will change when num_proteins value is changed
    Expected values are:
    - Plasmepsin I: 6.3E-5 M
    - Plasmepsin II: 4.9E-5 M
    - HAP: 1.1E-4 M
    :return:
    """

    # Concentrations in ppm
    conc_plm_1 = 754  # Plasmepsin I
    conc_plm_2 = 585  # Plasmepsin II
    conc_hap = 1_373  # Histo-aspartic protease

    # Initialise
    constants = Constants()

    # Convert ppm to Molar
    conc_plm_1 = round(constants._dv_ppm_to_molar(ppm=conc_plm_1), 6)
    conc_plm_2 = round(constants._dv_ppm_to_molar(ppm=conc_plm_2), 6)
    conc_hap = round(constants._dv_ppm_to_molar(ppm=conc_hap), 5)

    assert conc_plm_1 == 6.3e-5
    assert conc_plm_2 == 4.9e-5
    assert conc_hap == 1.1e-4


def test_compute_kobs():
    """
    Test computing the observed rate constant of haemoglobin degrading enzymes.
    Expected values are:
    - Plasmepsin I: 137 min-1
    - Plasmepsin II: 927 min-1
    - HAP: 6 min-1
    :return:
    """
    # Rate constants
    kcat_plm_1 = 2.3  # s-1 (ref 1, above)
    Km_plm_1 = 0.49e-6  # (ref 1, above)
    kcat_plm_2 = 11  # # s-1 (ref 2, above)
    Km_plm_2 = 2.6e-6  # (ref 2, above)
    kcat_hap = 0.1  # s-1 (ref 1, above).
    Km_hap = 0.0  # Km is unknown and is set to 0

    # Concentrations (from dv_ppm_to_molar)
    conc_plm_1 = 6.3e-5
    conc_plm_2 = 4.9e-5
    conc_hap = 1.1e-4

    # Initialise
    constants = Constants()

    # Convert to observed rate constant
    k_obs_plm_1 = constants._compute_kobs(kcat=kcat_plm_1,
                                          Km=Km_plm_1,
                                          enzyme_conc=conc_plm_1)
    k_obs_plm_2 = constants._compute_kobs(kcat=kcat_plm_2,
                                          Km=Km_plm_2,
                                          enzyme_conc=conc_plm_2)
    k_obs_hap = constants._compute_kobs(kcat=kcat_hap,
                                        Km=Km_hap,
                                        enzyme_conc=conc_hap)

    assert round(k_obs_plm_1) == 137
    assert round(k_obs_plm_2) == 627
    assert round(k_obs_hap) == 6


def test_compute_rate_hb_deg():
    """
    Test computing the observed rate constant for haemoglobin degradation
    Note: Values will change when num_proteins value is changed
    Expected value is 0.04 M.min-1
    :return:
    """

    constants = Constants()
    constants.compute_rate_hb_deg()

    assert round(constants.k_hb_deg, 2) == 0.04
