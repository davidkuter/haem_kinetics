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


def test_native_kcat_gluzman_ranking_not_garnie_fit():
    """PM II native kcat is not the peptide 11 s-1 (would invert Gluzman)."""
    c = Constants()
    assert c.k_enzymes['plm_2']['kcat'] == 11
    assert c.k_enzymes_native['plm_1']['kcat'] == c.k_enzymes['plm_1']['kcat']
    assert c.k_enzymes_native['plm_2']['kcat'] == c.k_enzymes_native['plm_1']['kcat']
    assert c.k_enzymes_native['plm_2']['kcat'] != c.k_enzymes['plm_2']['kcat']
    assert c.k_enzymes_native['fp_2']['kcat'] == c.k_enzymes['fp_2']['kcat']
    assert 'hap' not in c.k_enzymes_native
    assert 'plm_4' not in c.k_enzymes_native
    assert 'fp_3' not in c.k_enzymes_native


def test_k_htv_release_is_klemba_half_life_bound():
    """ln(2)/20 min — Klemba delivery t½ bound, not Garnie-fitted."""
    import math
    c = Constants()
    assert abs(c.k_htv_release - math.log(2) / 20.0) < 1e-15
