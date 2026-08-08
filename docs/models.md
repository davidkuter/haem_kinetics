# Haem kinetics models

Simulates haemoglobin (Hb) uptake into the *Plasmodium falciparum* digestive vacuole (DV), enzymatic release of haem, Fe(II)→Fe(III) oxidation, and detoxification to haemozoin (Hz). Experimental targets are Combrink-style heme fractionation time courses (fg Fe/cell), primarily **Dd2**.

**Time base:** simulation `t` is in **minutes** from a trophozoite offset of **16 h** post-invasion (`parasite age = 16 + t/60` hours).

## Modeling principles

Every addition to the ladder must be **mechanistically accountable** (chemistry / cell biology / cited measurement). We do **not** add numerical patches, efficiency fudges, or ad hoc rate caps to improve plots or stabilize the solver. Solver tolerances exist only to integrate the **stated** ODEs accurately. If the solution goes negative or total Fe drifts, that is a signal to fix the model or the integration settings — not to clip or reshape the RHS.

---

## Model pages

| Model | File | One-line summary |
|-------|------|------------------|
| [Degradation](models/degradation.md) | `degradation.py` | Uptake / Fe(II) release sandbox |
| [Model 1](models/model1.md) | `model1.py` | Linear uptake; PMs + FP2/3; full Fe speciation; no lipid |
| [Model 2](models/model2.md) | `model2.py` | Model 1 + lipid factor φ on Fe(III) rates |
| [Model 3](models/model3.md) | `model3.py` | Exponential uptake; enzymes track `f_exp`; φ on Hz |
| [Model 4](models/model4.md) | `model4.py` | Aqueous ⇄ lipid Fe(III); Hz from lipid at `k_hz` (no φ) |
| [Model 5](models/model5.md) | `model5.py` | Model 4 + logistic enzyme clock |
| [Model 6](models/model6.md) | `model6.py` | Model 5 + crystal-competent Fe(III) pool |

---

## Quick comparison

| Model | Hb uptake | Proteases | Enzyme schedule | Lipid / Fe(III) | Hz rate |
|-------|-----------|-----------|-----------------|-----------------|---------|
| Degradation | Exponential | PMs | `f_exp` | — | — |
| 1 | Linear | **PMs + FP2/3** | constant | None | `k_hz × [Fe3]` |
| 2 | Linear | PMs + FP2/3 | constant | φ on Fe3 rates | `k_hz × φ × [Fe3]` |
| 3 | `f_exp(t)` | PMs + FP2/3 | `f_exp` | φ on Fe3 rates | same as 2 |
| 4 | `f_exp(t)` | PMs + FP2/3 | `f_exp` | aq ⇄ lip | `k_hz × [Fe3]_lip` |
| 5 | `f_exp(t)` | PMs + FP2/3 | **logistic** | aq ⇄ lip | `k_hz × [Fe3]_lip` |
| 6 | `f_exp(t)` | PMs + FP2/3 | logistic | aq ⇄ lip ⇄ **xtal** | `k_hz × [Fe3]_xtal` |

**Incremental ladder:** one mechanistic change per step.

| Step | Problem in previous model | What this model changes |
|------|---------------------------|-------------------------|
| 1 | Need a minimal closed Fe path | Linear uptake + PMs + FP2/3 + ox + Hz |
| 2 | Model 1 drains free Fe³⁺; no basal Hm | Multiply Fe³⁺ rates by φ |
| 3 | Linear uptake / constant [E] mistimes trophozoite | `f_exp` uptake; proteases track `f_exp` |
| 4 | φ×`k_hz` wrong for lipid-mediated Hz | Explicit aq ⇄ lip; Hz from lip at full `k_hz` |
| 5 | Digestion still locked to uptake `f_exp` | **Logistic** `[E]` vs age |
| 6 | Full `k_hz` on lipid drains free haem | Crystal-competent pool; Hz from xtal only |

`v_dig` includes only proteases that liberate haem from Hb / haem-bearing globin (PMs + falcipains). Downstream peptidases are omitted. Falcipains are present from Model 1 onward (Degradation remains PMs-only as a sandbox).

---

## Shared framework

```mermaid
flowchart LR
  Host["Hb_RBC host pool"] -->|uptake| HbDV["Hb_DV"]
  HbDV -->|proteases| Fe2["Fe2PPIX"]
  Fe2 -->|"k_ox x O2"| Fe3["Fe3PPIX speciation"]
  Fe3 -->|crystallization| Hz["Haemozoin"]
```

**Common rules (all models unless a page says otherwise):**

- Host RBC haemoglobin (`conc_hb_rbc`) is an ODE state; uptake depletes it so total Fe ≈ **106 fg/cell**.
- Consumption rates are zero when their substrate is ≤ 0 (domain of the physical rate law). BDF defaults to `rtol=1e-8`, `atol=1e-12` so the **stated** ODEs are integrated accurately — not a change to the chemistry.
- Total Fe must stay ≈ **106 fg/cell** without clipping. Negatives or drift mean investigate the model or the solve.
- Integration defaults to SciPy **`BDF`**.
- DV species are integrated in **M** (fixed `vol_dv`) and converted to fg/cell for plotting.
- `[O2−]` = 0 (SOD) → Fe(III) reduction is off.
- Parameters: [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py)
- Experiment tables: [`experimental_data.py`](../haem_kinetics/components/experimental_data.py)

### Shared notation

| Symbol | Code / meaning |
|--------|----------------|
| `t` | Time in **minutes** from 16 h post-invasion |
| `[Hb]_DV` | `conc_hb_dv` — DV Hb as haem-equivalents (M) |
| `[Hb]_tet` | `[Hb]_DV` / 4 — tetramer basis for MM |
| `[Hb]_RBC` | `conc_hb_rbc` — remaining host Hb (M, RBC basis) |
| `v_up` | Uptake rate into DV (M haem-eq / min on `V_DV`) |
| `v_dig` | Digestion / haem-release rate (M haem-eq / min) |
| `v_ox` | Fe(II)→Fe(III) oxidation rate |
| `v_hz` | Haemozoin formation rate |
| `f_exp(t)` | `a · b · exp(b · t)` with `a` = 0.1578, `b` = 0.001102 |
| `s_log(t)` | Logistic enzyme scale vs parasite age (Models 5–6) |

Host mass balance:

```text
# Remaining host RBC Hb
d[Hb]_RBC / dt = − v_up · V_DV / V_RBC
```

### Shared physical / volume constants

| Constant | Value | Units | Description |
|----------|------:|-------|-------------|
| `V_RBC` | 90×10⁻¹⁵ | L | Volume of the host red blood cell |
| `V_DV` | 1×10⁻¹⁵ | L | Fixed digestive-vacuole volume used for M ↔ fg conversion |
| `f_lip` | 0.016 | — | Fractional volume of lipid nanospheres relative to the DV |
| `N_A` | 6.022×10²³ | mol⁻¹ | Avogadro's number |
| `N_prot` | 1.9×10⁸ | — | Average number of proteins per *P. falciparum* cell |
| `[Hb]_RBC,0` | ≈ 0.0211 | M | Uninfected RBC haemoglobin concentration (haem-equivalents) |
| Total Fe budget | ≈ 106 | fg/cell | Total iron inventory per infected RBC (`[Hb]_RBC,0 · V_RBC × 55.85`) |

### Shared enzyme inputs (PaxDB ppm → derived DV `[E]`)

**Do not treat `[E]` as an independent constant.** PaxDB reports abundance in ppm; molar DV concentration is derived:

`[E] = ppm × 10⁻⁶ × N_prot / (N_A · V_DV)`

With fixed `V_DV` this number is constant in time for a run, but it still depends on `ppm` and `V_DV`. Default ppm source: PaxDB **P. falciparum 3D7 — Whole organism, Dd2, SC (Tao, MCP, 2014)** (no DV-specific tissue). `fp_2` is falcipain-2a. Table `kcat` is s⁻¹; code uses `kcat [min⁻¹] = 60 × kcat [s⁻¹]`.

**kcat / Km citations:** see [`docs/enzyme_kinetics.md`](enzyme_kinetics.md). Defaults use Banerjee et al. *PNAS* 2002 Table 1 (PM I/II from Luker 1996; HAP & PM IV measured there) and Ramjee et al. *Biochem. J.* 2006 for falcipains. All are **peptide** assays, not native-Hb turnover.

| Constant | Value | Units | Description |
|----------|------:|-------|-------------|
| `ppm_plm_1` | 752 | — | PaxDB Tao 2014 Dd2 whole-organism abundance of plasmepsin-1 (PF3D7_1407900; input; `[E]` derived) |
| `kcat_plm_1` | 2.3 | s⁻¹ | Luker/Banerjee α33–34 peptide (native PM I) |
| `Km_plm_1` | 0.49×10⁻⁶ | M | Luker/Banerjee α33–34 peptide (native PM I) |
| `ppm_plm_2` | 1204 | — | PaxDB Tao 2014 Dd2 whole-organism abundance of plasmepsin-2 (PF3D7_1408000; input; `[E]` derived) |
| `kcat_plm_2` | 11 | s⁻¹ | Luker/Banerjee α33–34 peptide (native PM II) |
| `Km_plm_2` | 2.6×10⁻⁶ | M | Luker/Banerjee α33–34 peptide (native PM II) |
| `ppm_hap` | 1373 | — | PaxDB Tao 2014 Dd2 whole-organism abundance of HAP / PM3 (PF3D7_1408100; input; `[E]` derived) |
| `kcat_hap` | 0.05 | s⁻¹ | Banerjee 2002 Table 1 (native HAP, α33–34) |
| `Km_hap` | 0.29×10⁻⁶ | M | Banerjee 2002 Table 1 (native HAP, α33–34) |
| `ppm_plm_4` | 3139 | — | PaxDB Tao 2014 Dd2 whole-organism abundance of plasmepsin-4 (PF3D7_1407800; input; `[E]` derived) |
| `kcat_plm_4` | 1.05 | s⁻¹ | Banerjee 2002 Table 1 (recombinant PM IV, α33–34) |
| `Km_plm_4` | 0.33×10⁻⁶ | M | Banerjee 2002 Table 1 (recombinant PM IV, α33–34) |
| `ppm_fp_2` | 20.2 | — | PaxDB Tao 2014 Dd2 whole-organism abundance of falcipain-2a (PF3D7_1115700; input; `[E]` derived) |
| `kcat_fp_2` | 0.79 | s⁻¹ | Ramjee 2006 best FP-2 FRET peptide |
| `Km_fp_2` | 0.9×10⁻⁶ | M | Ramjee 2006 best FP-2 FRET peptide |
| `ppm_fp_3` | 23.5 | — | PaxDB Tao 2014 Dd2 whole-organism abundance of falcipain-3 (PF3D7_1115400; input; `[E]` derived) |
| `kcat_fp_3` | 0.204 | s⁻¹ | Ramjee 2006 FP-3 Leu-Arg FRET peptide |
| `Km_fp_3` | 4.0×10⁻⁶ | M | Ramjee 2006 FP-3 Leu-Arg FRET peptide |

Falcipains (fp_2/3) appear from Model 1 onward. `[E]_i` in equations is always the derived molarity from `ppm_i`. Generic MM haem-release rate:

```text
# Haem release from Hb (MM sum; 4 haem-eq per tetramer)
v_dig = 4 · Σ_i  (60 · kcat_i) · [E]_i,eff · [Hb]_tet / (Km_i + [Hb]_tet)
```

where `[E]_i,eff` is constant (Models 1–2), `f_exp(t)·[E]_i` (Models 3–4 / Degradation), or `s_log(t)·[E]_i` (Models 5–6).

---

## How to run

```bash
pip install -e .
python examples/run.py   # writes examples/model1.png … model6.png
```

```python
from haem_kinetics.models.model6 import Model6

model = Model6()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model6.png',
)
```

---

## Cross-cutting issues

1. **φ on `k_hz` (Models 2–3):** chemically inconsistent with lipid-mediated β-haematin; corrected in Model 4.
2. **Fixed `vol_dv = 1 fL`:** Combrink 2025 reports dynamic Dd2 lumen (~3.7 fL peak).
3. **Models 1–2:** full PaxDB [E] without an `f_exp` gate digests DV Hb almost instantly.
4. **NF54 digits** in `experimental_data.py` may not match Combrink 2025; prefer Dd2.
5. Peptide `kcat`/`Km` (Banerjee/Luker; Ramjee) applied to DV Hb is an approximation — see [`enzyme_kinetics.md`](enzyme_kinetics.md). Model 6 `K_xtal` is provisional. If peptide kinetics over-digest vs Dd2, introduce a dedicated Hb-efficiency step with protein-vs-peptide evidence rather than an uncited rescale.
6. ppm are Tao 2014 Dd2 whole-organism (not DV-specific).

---

## References

- Combrink et al., *Communications Biology* (2025): DV volume, uptake, basal Hb/Hm/Hz.
- Egan et al., *Malaria Journal* 11:337 (2012): lipid-mediated β-haematin kinetics.
- Enzyme kcat/Km: [`enzyme_kinetics.md`](enzyme_kinetics.md) (Banerjee 2002; Luker 1996; Ramjee 2006).
