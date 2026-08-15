"""Native-Hb vs globin enzyme assignment (Goldberg ordered pathway).

Used by Models 4a and 4b so the two encodings cannot drift.
"""

NATIVE_HB_ENZYMES = ('plm_1', 'plm_2', 'fp_2')
GLOBIN_ENZYMES = ('plm_1', 'plm_2', 'hap', 'plm_4', 'fp_2', 'fp_3')


def mm_haem_eq_rate(model, enzymes, k_table, conc_hb_haem, t) -> float:
    """Haem-equivalent MM rate (M/min): 4 · Σ kcat[E][tet] / (Km + [tet])."""
    conc_tet = conc_hb_haem / 4.0
    if conc_tet <= 0.0:
        return 0.0
    deg = 0.0
    for enzyme in enzymes:
        kcat = k_table[enzyme]['kcat'] * 60.0
        Km = k_table[enzyme]['Km']
        conc_enzyme = model._enzyme_conc(enzyme, t)
        denom = Km + conc_tet
        if denom <= 0.0:
            continue
        deg += kcat * conc_enzyme / denom
    return 4.0 * deg * conc_tet


def native_tetramer_rate(model, conc_hb_haem, t) -> float:
    """Nick / haem-release from native tetramer (PM I, PM II, FP-2)."""
    return mm_haem_eq_rate(
        model, NATIVE_HB_ENZYMES, model.const.k_enzymes_native, conc_hb_haem, t
    )


def globin_haem_release_rate(model, conc_globin_haem, t) -> float:
    """Haem release from nicked globin (all six; peptide MM table)."""
    return mm_haem_eq_rate(
        model, GLOBIN_ENZYMES, model.const.k_enzymes, conc_globin_haem, t
    )
