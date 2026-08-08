from haem_kinetics.components.constants import Constants


def test_compute_conc_hb_rcb():
    """Uninfected RBC Hb as haem-equivalents (~21 mM)."""
    c = Constants()
    assert round(c.conc_hb_rbc, 3) == 0.021


def test_dv_ppm_to_molar_tao_dd2():
    """PaxDB Tao 2014 Dd2 ppm -> DV M with current N_prot and V_DV."""
    c = Constants()

    def expected(ppm: float) -> float:
        return ppm * 1e-6 * c.num_prots / (c.avogadro * c.vol_dv)

    assert abs(c.conc_enzymes["plm_1"] - expected(752)) < 1e-15
    assert abs(c.conc_enzymes["plm_2"] - expected(1204)) < 1e-15
    assert abs(c.conc_enzymes["hap"] - expected(1373)) < 1e-15
    assert abs(c.conc_enzymes["plm_4"] - expected(3139)) < 1e-15
    assert abs(c.conc_enzymes["fp_2"] - expected(20.2)) < 1e-15
    assert abs(c.conc_enzymes["fp_3"] - expected(23.5)) < 1e-15
