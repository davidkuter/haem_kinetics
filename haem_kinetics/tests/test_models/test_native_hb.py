from haem_kinetics.models.model4a import Model4a
from haem_kinetics.models.model4b import Model4b
from haem_kinetics.models.native_hb import (
    GLOBIN_ENZYMES,
    NATIVE_HB_ENZYMES,
    native_tetramer_rate,
)


def test_native_competent_set_is_pmi_pmii_fp2():
    assert NATIVE_HB_ENZYMES == ('plm_1', 'plm_2', 'fp_2')
    assert GLOBIN_ENZYMES == ('plm_1', 'plm_2', 'hap', 'plm_4', 'fp_2', 'fp_3')
    assert 'hap' not in NATIVE_HB_ENZYMES
    assert 'plm_4' not in NATIVE_HB_ENZYMES
    assert 'fp_3' not in NATIVE_HB_ENZYMES


def test_model4a_expands_four_value_init_with_globin_zero():
    assert Model4a()._expand_init([0.018, 0.0, 0.0, 0.36]) == [
        0.018, 0.0, 0.0, 0.0, 0.36
    ]


def test_4a_and_4b_share_native_rate():
    hb = 1.0e-3
    t = 240.0
    m4a = Model4a()
    m4b = Model4b()
    m4a.initial_values['conc_hb_dv'] = hb
    m4b.initial_values['conc_hb_dv'] = hb
    r_a = native_tetramer_rate(m4a, hb, t)
    r_b = native_tetramer_rate(m4b, hb, t)
    assert r_a > 0.0
    assert abs(r_a - r_b) < 1e-18
    assert abs(m4b._hb_removal(t) - r_a) < 1e-18
    assert abs(m4a._nick_rate(t) - r_a) < 1e-18


def test_model4a_scores_assay_hb_as_native_plus_globin():
    model = Model4a()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    model.score_vs_experiment()
    df = model.concentrations
    assay = df['conc_hb_dv'] + df['conc_hb_globin']
    assert abs(float(assay.iloc[0]) - float(df['conc_hb_dv'].iloc[0])) < 0.05
    assert model.fit_metrics['Hb']['n'] == 9
    assert 'conc_hb_assay' not in df.columns
