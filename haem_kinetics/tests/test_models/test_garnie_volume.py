from haem_kinetics.models.helpers import garnie_dd2_vol_dv_fl, garnie_dd2_vol_dv_L
from haem_kinetics.models.model3 import Model3


def test_garnie_dd2_vol_matches_reconstructed_means():
    """Gompertz+collapse reconstruction vs xlsx/paper landmarks."""
    assert abs(garnie_dd2_vol_dv_fl(20.0) - 0.675) < 0.01
    assert abs(garnie_dd2_vol_dv_fl(24.0) - 2.278) < 0.01
    assert abs(garnie_dd2_vol_dv_fl(32.0) - 3.557) < 0.02
    assert abs(garnie_dd2_vol_dv_fl(36.0) - 0.763) < 0.01
    assert abs(garnie_dd2_vol_dv_fl(46.0)) < 1e-12


def test_garnie_dd2_vol_positive_in_trophozoite_window():
    vol, _dvol = garnie_dd2_vol_dv_L(0.0)
    assert vol > 0.0
    vol_end, _dvol = garnie_dd2_vol_dv_L(1700.0)  # ~44.3 h
    assert vol_end > 0.0


def test_model3_preserves_seed_fg_and_total_fe():
    model = Model3()
    model.run(t=[0, 1700], init=[0.018, 0.0, 0.0, 0.36], t_eval=range(0, 1700, 20))
    df = model.concentrations
    start = float(df.iloc[0][['conc_hb_dv', 'conc_fe2pp', 'conc_fe3pp', 'conc_hz']].sum())
    assert abs(start - (0.018 + 0.36) * 55.85) < 0.05
    tot = df[[c for c in df.columns if c.startswith('conc_')]].sum(axis=1)
    budget = model.const.total_fe_fg_cell
    assert abs(float(tot.iloc[0]) - budget) < 0.5
    assert abs(float(tot.iloc[-1]) - budget) < 0.5
