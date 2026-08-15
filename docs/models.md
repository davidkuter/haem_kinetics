# Haem kinetics models

Simulates haemoglobin (Hb) uptake into the *Plasmodium falciparum* digestive vacuole (DV), enzymatic release of haem, Fe(II)→Fe(III) oxidation, and detoxification to haemozoin (Hz). Experimental targets are Garnie et al. (*Commun. Biol.* 2025) **heme fractionation** time courses (fg Fe/cell), primarily **Dd2**. Those numbers come from saponin-isolated trophozoites, not isolated DVs; Garnie still treats Hb, Hm, and Hz as DV-localized. Assay vs pHrodo: [garnie_fractionation.md](garnie_fractionation.md).

**Time base:** simulation `t` is in **minutes** from a trophozoite offset of **16 h** post-invasion (`parasite age = 16 + t/60` hours).

## Modeling principles

Every addition to the ladder must be **mechanistically accountable** (chemistry / cell biology / cited measurement). We do **not** add numerical patches, efficiency fudges, or ad hoc rate caps to improve plots or stabilize the solver. Solver tolerances exist only to integrate the **stated** ODEs accurately. If the solution goes negative or total Fe drifts, that is a signal to fix the model or the integration settings — not to clip or reshape the RHS.

---

## Model pages (active ladder)

| Model | File | One-line summary |
|-------|------|------------------|
| [Degradation](models/degradation.md) | `degradation.py` | Uptake / Fe(II) release sandbox |
| [Model 1](models/model1.md) | `model1.py` | Linear uptake; PMs + FP2/3; full Fe speciation; no lipid |
| [Model 2a / 2b](models/model2.md) | `model2.py`, `model2b.py` | Empirical `f_exp` × remaining host (2b: two-phase at Fig. 5B 29 h, continuous) |
| [Model 3](models/model3.md) | `model3.py` | Model 2b + Garnie Fig. 3 `s_PM(t)` amount schedule |
| [Model 4a / 4b](models/model4.md) | `model4a.py`, `model4b.py` | Native-Hb-competent enzymes (ordered pathway vs single pool) |
| [Model 5](models/model5.md) | `model5.py` | Model 4b + inaccessible inner-vesicle Hb cargo (already in the DV) |
| [Model 6](models/model6.md) | `model6.py` | Model 5 + `k_release(t) ∝ s_PM(t)` (Garnie Fig. 3 lysis clock) |
| [Model 7](models/model7.md) | `model7.py` | Model 6 + aqueous ⇄ lipid Fe(III); Hz from lipid at `k_hz` |
| [Model 8](models/model8.md) | `model8.py` | Model 7 + interfacial (crystal-competent) Fe(III); Hz from xtal |
| [Model 9](models/model9.md) | `model9.py` | Model 8 + crystal-area growth: `v_hz ∝ [Fe3]_xtal · (n_Hz / n_Hz_start)^{2/3}` |
| [Model 10](models/model10.md) | `model10.py` | **What-if** (not mechanistic): Model 9 with Garnie-tuned `k_release` and `K_xtal` (`f_exp` unchanged) |

**Archived prior ladder** (φ → f_exp+[E] → lipid → logistic → xtal): [`docs/models/legacy/`](models/legacy/README.md) and `haem_kinetics/models/legacy/`.

---

## Quick comparison

| Model | Hb uptake | Proteases | Enzyme schedule | Lipid / Fe(III) | Hz rate |
|-------|-----------|-----------|-----------------|-----------------|---------|
| Degradation | Exponential | PMs | `f_exp` | — | — |
| 1 | Linear | **PMs + FP2/3** | `n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 2a | empirical `f_exp(t)` × remaining host | PMs + FP2/3 | `n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 2b | two-phase `f_exp` × remaining host (29 h break) | PMs + FP2/3 | `n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 3 | Model 2b two-phase `f_exp` | PMs + FP2/3 | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 4a | empirical `f_exp(t)` | **PM I/II/FP-2 on native**; all six peptide MM on globin | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 4b | empirical `f_exp(t)` | **PM I/II/FP-2 native rate only** (lumped) | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 5 | `f_exp` into **inner-vesicle cargo in the DV**; first-order lysis to lumen | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 6 | same as 5; **`k_release(t) ∝ s_PM(t)`** | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 7 | same as 6 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | **aq ⇄ lip** | `k_hz × [Fe3]_lip` |
| 8 | same as 6 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip ⇄ **xtal** | `k_hz × [Fe3]_xtal` |
| 9 | same as 8 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | same | `k_hz × [Fe3]_xtal × (n_Hz / n_Hz_start)^{2/3}` |
| 10 | same as 9; **slower plateau `k_release` (fit)** | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | same; **`K_xtal` fit** | same area law as 9 |

**Incremental ladder:** one mechanistic change per step.

| Step | Problem in previous model | What this model changes |
|------|---------------------------|-------------------------|
| 1 | Need a minimal closed Fe path | Linear uptake + PMs + FP2/3 + ox + Hz |
| 2a / 2b | Linear uptake leaves most Fe in host | Empirical `f_exp(t)` × remaining host; enzyme **amount** unchanged. **2a** = one phase (comparison). **2b** = two phases, break at Garnie Fig. 5B 29 h, `f_exp` continuous at the join (ladder default). Models 3–10 inherit Model 2b |
| 3 | Full PaxDB amount from `t` = 0 collapses DV Hb | Garnie Fig. 3 `s_PM(t)` on that amount; `f_exp` unchanged |
| 4a / 4b | Peptide `kcat` on native tetramer; PM II dominates `Vmax` | Goldberg ordered pathway: native-competent set vs peptide MM on globin |
| 5 | Lumen `Vmax` still ≫ uptake; standing Hb ~0 | Inaccessible **inner-vesicle** cargo in the DV before 4b lumen proteases; `k_release` from Klemba lysis bound |
| 6 | Constant `k_release` makes assay Hb track `v_up(t)` | `k_release(t) = k_Klemba · s_PM(t)` (Fig. 3 clock; plateau keeps Klemba `t½`) |
| 7 | Single Fe(III) pool drains Hm at lipid-assay `k_hz` | Aqueous ⇄ lipid Fe(III); `v_hz` on lipid pool at full `k_hz` (no `φ`) |
| 8 | Lipid pool still drains Hm at full `k_hz` | Interfacial / crystal-competent Fe(III); `v_hz` on xtal at full `k_hz` |
| 9 | First-order `k_hz` does not grow Egan’s accelerating structure | `v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{2/3}` (amount, sphere geometry) |
| 10 | What-if: can Hb/Hm match if two knobs are freed? | **Not chemistry** — slower plateau `k_release` and larger `K_xtal` vs Garnie Dd2; inherits Model 9 area; `f_exp` unchanged |

Model 2a’s `f_exp` is a provisional first-order remaining-host clock with `a`, `b` fit to cumulative **fractionation** Fe (`Hb+Hm+Hz`; see [model2.md](models/model2.md) and [garnie_fractionation.md](garnie_fractionation.md)) — not lumen volume, not pHrodo Fig. 2C, not cytostome kinetics, and not “delivery to the parasite” with a later DV step. On Dd2, ODE DV Fe **R² = 0.612**, RMSE 15.89 fg (~36 fg still in host at 44 h). Model 2b is the same constitutive class with a faster late phase after Garnie Fig. 5B’s 29 h break, `f_exp` continuous at the join (DV Fe RMSE 2.56, R² = 0.990, host ~5 fg at 44 h); **Models 3–10 inherit 2b**. The **form** is general; `a` and `b` can be refit to another strain. Variable `V_DV(t)` (`variable_dv_volume`) is **shared bookkeeping** from Model 1 onward (not a numbered ladder step). The present numerical schedule is Garnie Dd2 **aqueous lumen**. Model 3’s one change is the blot-derived amount clock (see [model3.md](models/model3.md)). Model 4’s one change is native vs peptide substrate (two encodings: [model4.md](models/model4.md)). Model 5’s one change is inaccessible **inner-vesicle** cargo already in the DV on the **4b** pathway (see [model5.md](models/model5.md)); 4a’s globin pool is not carried because it never accumulated. A native-Hb `kcat_app` was attempted and still lacks moles of enzyme ([enzyme_kinetics.md](enzyme_kinetics.md)); that gap is not filled with a Garnie-fitted scalar. Model 6’s one change is `k_release(t) ∝ s_PM(t)` so lysis follows the same Fig. 3 clock (see [model6.md](models/model6.md)). Model 7’s one change is Egan lipid partition of Fe(III) (see [model7.md](models/model7.md)), not `φ` on `k_hz`. Model 8’s one change is interfacial (crystal-competent) Fe(III) so bulk NLB haem is not the `k_hz` substrate (see [model8.md](models/model8.md)). Model 9’s one change is crystal-area growth of that `v_hz` (see [model9.md](models/model9.md)): amount, not lumen concentration; `n_Hz_start` is the 20 fg seed; `2/3` is sphere geometry. Model 10 is a **what-if** on that topology ([model10.md](models/model10.md)), not a cited mechanism: slower plateau `k_release` and slightly larger `K_xtal`, chosen on the **2a** inventory and not refit. It does **not** retune `f_exp`.

`v_dig` includes only proteases that liberate haem from Hb / haem-bearing globin (PMs + falcipains). Downstream peptidases are omitted. Falcipains are present from Model 1 onward (Degradation remains PMs-only as a sandbox).

---

## Fit vs Garnie Dd2 (tracking)

Scores are **diagnostics**, not an objective to minimize with extra terms. A ladder step counts as an improvement only if the change is mechanistic **and** the relevant series move toward the assay.

**Protocol:** `t = [0, 1700]` min from 16 h, `t_eval` step 20 min, init `[0.018, 0, 0, 0.36]`. Model interpolated onto Dd2 ages 20–44 h (`n` = 9). Hm is scored as `conc_fe3pp` (Models 1–6), aq + lip (Model 7), or aq + lip + xtal (Models 8–10). Model 4a scores assay Hb as native + globin; Models 5–10 as inner-vesicle + lumen Hb. Init Hb seeds the inner-vesicle pool in Models 5–10. What the tubes measure: [garnie_fractionation.md](garnie_fractionation.md).

| Symbol | Definition |
|--------|------------|
| RMSE | √ mean((pred − obs)²) (fg/cell) |
| MAE | mean(\|pred − obs\|) (fg/cell) |
| mean signed error | mean(pred − obs) (fg/cell). **Not** mean squared error. Negative = model low. Equals −MAE only if every residual has the same sign. |
| χ²_red | mean(((pred − obs)/SEM)²); ~1 means residuals match assay scatter |

Mean signed error is listed on each model page. Recompute with `python examples/run.py`.

| Model | RMSE Hb | χ²_red Hb | RMSE Hm | χ²_red Hm | RMSE Hz | χ²_red Hz | RMSE DV Fe | χ²_red DV Fe |
|-------|--------:|----------:|--------:|----------:|--------:|----------:|-----------:|-------------:|
| 1 | 1.91 | 37.13 | 3.63 | 284 | 37.66 | 9.84 | 42.76 | 13.70 |
| 2a | 1.91 | 37.13 | 3.38 | 231 | 11.63 | 0.57 | 15.89 | 1.17 |
| 2b | 1.91 | 37.13 | 3.21 | 215 | 5.32 | 0.48 | 2.56 | 0.09 |
| 3 | 1.91 | 37.13 | 3.21 | 215 | 5.32 | 0.48 | 2.56 | 0.09 |
| 4a | 1.91 | 37.13 | 3.21 | 215 | 5.32 | 0.48 | 2.56 | 0.09 |
| 4b | 1.91 | 37.13 | 3.21 | 215 | 5.32 | 0.48 | 2.56 | 0.09 |
| 5 | 0.65 | 4.40 | 3.20 | 216 | 3.97 | 0.31 | 2.56 | 0.09 |
| 6 | 0.33 | 0.74 | 3.20 | 217 | 3.59 | 0.22 | 2.56 | 0.09 |
| 7 | 0.33 | 0.74 | 3.13 | 207 | 3.54 | 0.22 | 2.56 | 0.09 |
| 8 | 0.33 | 0.74 | 3.31 | 102 | 4.10 | 0.17 | 2.56 | 0.09 |
| 9 | 0.33 | 0.74 | 1.14 | 53.3 | 2.30 | 0.08 | 2.56 | 0.09 |
| 10 | 2.33 | 52.7 | 1.24 | 1.29 | 2.82 | 0.11 | 2.56 | 0.09 |

**Reading the ladder:** Model 2a improved Hz and internalized Fe relative to Model 1, but still leaves ~36 fg in the host at 44 h (inventory **R² = 0.612**). Model 2b (ladder default; Models 3–10 inherit it) uses a faster late `f_exp` after Fig. 5B’s 29 h break, still × remaining host, continuous at the join: DV Fe RMSE 15.89 → 2.56 (R² = 0.990), host ~5 fg at 44 h, Hz RMSE 11.63 → 5.32. The last gulp is still a little short. Hb stayed collapsed through Model 4 (lumen `Vmax` ≫ uptake). Shared `variable_dv_volume` did not move M1/M2 fg scores (amount-linear rates). Model 3 matches Model 2b: lag `s_PM(20 h) ≈ 0.41` of plateau is already enough enzyme to collapse DV Hb. Models 4a and 4b also match Model 3. Model 5 is the first step that moves Hb (χ²_red 37 → 4.4): assay Hb tracks inner-vesicle cargo while lumen Hb stays ~0. DV Fe is unchanged (same Model 2b internalized inventory). Model 6 clocks lysis with `s_PM(t)` so standing cargo is `v_up / (k · s_PM)` rather than `v_up / k` ([model6.md](models/model6.md)): Hb χ²_red 4.4 → 0.74. Model 7 splits Fe(III) aq ⇄ lip and crystallizes from the lipid pool at full literature `k_hz` (no `φ`). Hm barely moves (χ²_red 217 → 207). Model 8 puts `k_hz` on interfacial Fe only: Hm signed error flips from −2.88 to +2.60 fg (high, not drained). `K_xtal = 3δ/R = 0.08` from cited NLB radius and film thickness ([model8.md](models/model8.md)) — do not retune it. Model 9 scales that `v_hz` by growing crystal area ([model9.md](models/model9.md)): Hb and DV Fe stay at Model 8; Hm RMSE 3.31 → 1.14 (signed +2.60 → 0.00), Hz RMSE 4.10 → 2.30. Do not fit the `2/3` exponent. Model 10 is a **what-if** on the same ODEs ([model10.md](models/model10.md)): plateau `k_release` and `K_xtal` were chosen on **2a** and are not refit; on 2b plus the Model 6 clock they overshoot Hb (χ²_red 0.74 → 53). Ranked next *ladder* step if 2b’s last ~5 fg is still short: Garnie Fig. 5B Dd2 phases as `v_up` (0.9 then 4.8 fg/h), not a new `f_exp` `b`.

Code: [`haem_kinetics/components/fit_metrics.py`](../haem_kinetics/components/fit_metrics.py); `model.score_vs_experiment()`.

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
- Lumen DV species are integrated in **M at variable `V_DV(t)`** (`variable_dv_volume`) and converted to fg/cell as `C·V_DV(t)`. `constants.vol_dv = 1 fL` is the **reference** volume for the init API and PaxDB amount `n_E = [E]_1fL · 1 fL`. Init `[0.018, 0, 0, 0.36]` stays 1 fL-reference M so the fg seed is unchanged; `run()` rescales lumen species to true M at `V(t=0)`. Model 5–10 inner-vesicle cargo (`conc_hb_htv`) is an **amount** (`AMOUNT_SPECIES`): inside the DV but not pHrodo aqueous lumen, encoded as M at `V_ref`, not diluted by `dV_DV/dt`, converted as `C·V_ref`. The current lumen schedule is Garnie Dd2 (`variable_dv_volume_L`).
- `[O2−]` = 0 (SOD) → Fe(III) reduction is off.
- Parameters: [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py)
- Experiment tables: [`experimental_data.py`](../haem_kinetics/components/experimental_data.py) ([garnie_fractionation.md](garnie_fractionation.md))

### Shared notation

| Symbol | Code / meaning |
|--------|----------------|
| `t` | Time in **minutes** from 16 h post-invasion |
| `[Hb]_HTV` | `conc_hb_htv` — Model 5–10 inner-vesicle cargo already in the DV (amount as M at `V_ref`; not pHrodo lumen) |
| `[Hb]_DV` | `conc_hb_dv` — protease-accessible lumen Hb as haem-equivalents (M); Model 4a: native tetramer; 4b/5–10: lumen Hb |
| `[Hb]_globin` | `conc_hb_globin` — Model 4a nicked globin (haem still protein-bound) |
| `[Hb]_tet` | `[Hb]_DV` / 4 — tetramer basis for MM |
| `[Hb]_RBC` | `conc_hb_rbc` — remaining host Hb (M, RBC basis) |
| `[Fe3]_aq`, `[Fe3]_lip` | Model 7–10 aqueous and bulk-lipid Fe(III) (lumen-basis M) |
| `[Fe3]_xtal` | Model 8–10 interfacial / crystal-competent Fe(III); assay Hm = aq + lip + xtal |
| `v_up` | Uptake rate into DV (M haem-eq / min on `V_DV`) |
| `v_dig` | Digestion / haem-release rate (M haem-eq / min) |
| `v_ox` | Fe(II)→Fe(III) oxidation rate |
| `v_hz` | Haemozoin formation rate |
| `f_exp(t)` | Empirical `a · b · exp(b · t)` uptake; first-order remaining host. Model 2a = one phase; Model 2b / 3–10 = two-phase at Fig. 5B 29 h. Fit to cumulative **fractionation** Fe, not pHrodo |
| `a_e, b_e, b_l` | Model 2b two-phase `f_exp` (still × remaining host); break at Fig. 5B 29 h; `a_l` from continuity of `f_exp` at the join |
| `k_release` | Model 5 constant inner-vesicle lysis (Klemba `t½` < 20 min bound). Models 6–9: `k_htv_release · s_PM(t)`. Model 10: that clock with a fitted plateau scale |
| `V_DV(t)` | Variable DV lumen (`variable_dv_volume`); current schedule: Garnie Dd2 Gompertz then linear collapse |
| `s_PM(t)` | Garnie Fig. 3 relative PM amount (Model 3–10 enzyme amount; Model 6–10 also lysis); plateau 40–44 h = 1 |

Host mass balance:

```text
# Remaining host RBC Hb
d[Hb]_RBC / dt = − v_up · V_DV(t) / V_RBC
```

Concentration ODEs include dilution `−C · (dV/dt) / V` so `C·V` is conserved when the lumen grows or collapses. This is geometry, not extra chemistry. Uptake is **not** `C·dV/dt` and **not** assay `dF/dt`.

### Shared physical / volume constants

| Constant | Value | Units | Description |
|----------|------:|-------|-------------|
| `V_RBC` | 90×10⁻¹⁵ | L | Volume of the host red blood cell |
| `V_DV,ref` | 1×10⁻¹⁵ | L | Reference DV volume (init API and PaxDB `n_E`) |
| `V_DV(t)` | variable | L | Aqueous lumen (`variable_dv_volume`); current schedule is Garnie Dd2 (`variable_dv_volume_L`) |
| `f_lip` | 0.016 | — | Fractional volume of lipid nanospheres relative to the DV |
| `N_A` | 6.022×10²³ | mol⁻¹ | Avogadro's number |
| `N_prot` | 1.9×10⁸ | — | Average number of proteins per *P. falciparum* cell |
| `[Hb]_RBC,0` | ≈ 0.0211 | M | Uninfected RBC haemoglobin concentration (haem-equivalents) |
| Total Fe budget | ≈ 106 | fg/cell | Total iron inventory per infected RBC (`[Hb]_RBC,0 · V_RBC × 55.85`) |

### Shared enzyme inputs (PaxDB ppm → derived DV `[E]`)

**Do not treat `[E]` as an independent constant.** PaxDB reports abundance in ppm; molar DV concentration is derived:

`[E] = ppm × 10⁻⁶ × N_prot / (N_A · V_DV(t))`

PaxDB sets enzyme **amount** `n_E`; molarity follows the lumen. Default ppm source: PaxDB **P. falciparum 3D7 — Whole organism, Dd2, SC (Tao, MCP, 2014)** (no DV-specific tissue). `fp_2` is falcipain-2a. Table `kcat` is s⁻¹; code uses `kcat [min⁻¹] = 60 × kcat [s⁻¹]`.

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

where `[E]_i,eff` is `n_E,i / V_DV(t)` on Models 1–2 and `s_PM(t) · n_E,i / V_DV(t)` on Model 3–10. Degradation’s sandbox still uses `f_exp`-scaled PMs. Models 1–3 sum all six peptide terms on DV Hb. Model 4 uses `k_enzymes_native` on the tetramer (PM I, PM II, FP-2 only); 4a applies the peptide table to nicked globin; 4b and 5–10 lump haem release with that native rate.

---

## How to run

```bash
pip install -e .
python examples/run.py   # plots + RMSE / χ²_red vs Garnie Dd2
```

```python
from haem_kinetics.models.model3 import Model3

model = Model3()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model3.png',
)
```

Legacy imports: `from haem_kinetics.models.legacy import LegacyModel3` (etc.).

---

## Cross-cutting issues

1. **`variable_dv_volume` is shared bookkeeping**, not a ladder win: Garnie et al. 2025 measure a dynamic Dd2 **aqueous lumen** (~3.7 fL peak) by pHrodo, which currently supplies `V_DV(t)`. Amount-linear rates make fg scores insensitive to that volume. Uptake is not `C·dV/dt`, not pHrodo intensity, and not assay `dF/dt` ([garnie_fractionation.md](garnie_fractionation.md)).
2. **Standing DV Hb stayed collapsed through Model 4** because lumen `Vmax` ≫ uptake. Model 5 puts internalized Hb in inaccessible **inner-vesicle** cargo already in the DV first ([model5.md](models/model5.md)); lumen native still collapses. That is not a re-targeting of `f_exp` from parasite to DV. A literature native-Hb `kcat_app` was attempted and still lacks moles of enzyme ([enzyme_kinetics.md](enzyme_kinetics.md)). Do not invent that number from Garnie fg or park haem. Model 6 clocks lysis with `s_PM(t)` ([model6.md](models/model6.md)). Model 7 addresses Hm via lipid partition ([model7.md](models/model7.md)), not by slowing `k_hz` with `φ`. Model 8 puts `k_hz` on interfacial Fe only ([model8.md](models/model8.md)).
3. **NF54 digits** in `experimental_data.py` may not match Garnie 2025; prefer Dd2. `s_PM` itself is from NF54 Fig. 3 blots (the published PM time course).
4. Peptide `kcat`/`Km` (Banerjee/Luker; Ramjee) applied to DV Hb is an approximation — see [`enzyme_kinetics.md`](enzyme_kinetics.md).
5. ppm are Tao 2014 Dd2 whole-organism (not DV-specific).
6. Model 2b / 3–10 `f_exp` is empirical (cumulative fractionation Fe, first-order remaining host, two-phase at Fig. 5B 29 h, continuous at the join), not cytostome kinetics and not pHrodo Fig. 2C ([model2.md](models/model2.md), [garnie_fractionation.md](garnie_fractionation.md)). Model 2a is the one-phase comparison. Do **not** retune `a, b` in Model 10.
7. **Late Hz on this ladder is mostly inventory, not crystal lag.** With Model 2b, ~5 fg remains in the host at 44 h (DV Fe RMSE 2.56). Ranked next ladder step if that last gulp is still short: Garnie Fig. 5B Dd2 phases as `v_up` (0.9 then 4.8 fg/h); see [model10.md](models/model10.md#next-mechanisms-ranked-not-this-page).

---

## References

- Assay vs interpretation: [garnie_fractionation.md](garnie_fractionation.md) (saponin pellet, not isolated DVs; pHrodo is a different experiment).
- Garnie LF, Egan TJ, Wicht KJ. Heme processing in the malaria parasite, *Plasmodium falciparum*: a time-dependent basal-level analysis. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z) · [Nature full text](https://www.nature.com/articles/s42003-025-08991-z) · [PDF](https://www.nature.com/articles/s42003-025-08991-z.pdf) · [Figshare raw data](https://doi.org/10.6084/m9.figshare.28801805)
- Combrinck JM, Fong KY, Gibhard L, Smith PJ, Wright DW, Egan TJ. Optimization of a multi-well colorimetric assay to determine haem species in *Plasmodium falciparum*. *Malar. J.* (2015) 14:253. [doi:10.1186/s12936-015-0729-9](https://doi.org/10.1186/s12936-015-0729-9)
- Egan TJ, Chen JY, de Villiers KA, et al. Haemozoin (β-haematin) biomineralization requires both a lipid medium and an accelerating structure to promote haem dimerization. *Malaria Journal* (2012) 11:337. [doi:10.1186/1475-2875-11-337](https://doi.org/10.1186/1475-2875-11-337)
- Enzyme kcat/Km: [`enzyme_kinetics.md`](enzyme_kinetics.md) (Banerjee 2002; Luker 1996; Ramjee 2006 — DOIs listed there).
