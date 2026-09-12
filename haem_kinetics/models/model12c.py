from haem_kinetics.models.model12b import Model12b
from haem_kinetics.models.native_hb import globin_haem_release_rate


class Model12c(Model12b):
    """
    Full Myburgh (2023) replication with our HTV layer only to hold Hb.

    Reproduces Myburgh's working configuration as closely as our code allows:
      * constant 1 fL lumen (variable_dv_volume = False) — his Models 1–3
        assumption; his Model 4 showed variable volume leaves the mass profile
        unchanged (our Models 10/11 confirm this for first-order kinetics);
      * empirical exponential uptake (Model 12b);
      * six proteases in parallel with constant enzyme concentration and the
        peptide-substrate MM table (Banerjee/Luker plasmepsins, Ramjee
        falcipains — the same sources as his Table 4.3), no s_PM schedule;
      * aqueous-hematin crystallisation with the lipid partition buffer (12a).

    Myburgh could not hold assay Hb above ~10⁻⁶ fg because he had no protected
    Hb pool. The only element retained from our ladder is the HTV inner-vesicle
    cargo (released at a constant first-order rate, not s_PM-clocked) so assay
    Hb = HTV + lumen is non-negligible. This isolates whether that single
    addition is what his pathway was missing for Hb.

    Host boundary condition (Myburgh's, overriding 12b): [Hb_RBC] is held
    **constant** — "any decrease in Hb in the RBC due to uptake is compensated
    for by the decrease in RBC cytoplasm volume" (§4.4.1.1) — an infinite RBC
    reservoir. The empirical uptake therefore runs smoothly through 44 h with
    no finite-host cliff. host+DV Fe is **not** conserved (unlike 12b and the
    rest of the ladder); the DV total is the experimental target.
    """

    variable_dv_volume = False

    def __init__(self, model_name: str = 'Model 12c'):
        super().__init__(model_name=model_name)

    def _d_host_from_dv_uptake(self, uptake_dv_m_per_min, t=None):
        # Myburgh constant-[Hb_RBC] reservoir: host does not deplete, so the
        # empirical uptake is never truncated (no late cliff). Not conserving.
        return 0.0

    def _enzyme_conc(self, enzyme, t):
        ppm = self.const.ppm_enzymes[enzyme]
        return self.const._dv_ppm_to_molar(ppm, vol_dv=self._vol_dv(t))

    def _hb_removal(self, t):
        hb = self._nonneg(self.initial_values['conc_hb_dv'])
        return globin_haem_release_rate(self, hb, t)

    def _k_release(self, t):
        return self.const.k_htv_release
