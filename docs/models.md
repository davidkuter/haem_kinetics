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
| [Model 9a/9b/9c](models/model9.md) | `model9a.py`, `model9b.py`, `model9c.py` | Model 8 + crystal-area growth: exponent 2/3 (sphere), 1/2 (rod), 1/3 (needle) |
| [Model 10](models/model10.md) | `model10.py` | Model 9a + NLB/Hz compartmentalization (lip/xtal/Hz not concentrated by lumen collapse) |
| [Model 11](models/model11.md) | `model11.py` | Model 9a + explicit **two-compartment volumes** (aqueous lumen shrinks, NLB constant); identical fg curves to 9a — volume-invariant |
| [Model 12a/12b/12c](models/model12.md) | `model12a.py`, `model12b.py`, `model12c.py` | **Myburgh 2023 Hm/Hz (Hz from aqueous hematin) on our HTV Hb**; 12a best combined fit. 12b/12c = Myburgh empirical uptake |
| [Model 13](models/model13.md) | `model13.py` | Model 12b uptake on an **upper-range MCHC Fe budget** (36 g/dL → ~112 fg); host survives → removes 12b's 44 h host-exhaustion cliff, one cited parameter |
| [Model 14a / 14b](models/model14.md) | `model14a.py`, `model14b.py` | Model 13 with **HTV release decoupled from `s_PM`** (constant lysis); removes the early Hb hump. 14a t½ 20 min, 14b t½ 30 min — bracket shows release must accelerate |
| [Model 15a / 15b](models/model15.md) | `model15a.py`, `model15b.py` | Lower Model 13's **early Hb amounts** (not flatten the shape); **not** Garnie Fig. 5B. 15a softer lysis clock (best Hb χ² 0.43); 15b Elliott 24–30 h cytostome gate (starves early Fe) |
| [Model 99](models/model99.md) | `model99.py` | **What-if** (not mechanistic): Model 9a with Garnie-tuned `k_release` and `K_xtal` (`f_exp` unchanged) |

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
| 9a | same as 8 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | same | `k_hz × [Fe3]_xtal × (n_Hz / n_Hz_start)^{2/3}` (sphere) |
| 9b | same as 8 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | same | `k_hz × [Fe3]_xtal × (n_Hz / n_Hz_start)^{1/2}` (rod) |
| 9c | same as 8 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | same | `k_hz × [Fe3]_xtal × (n_Hz / n_Hz_start)^{1/3}` (needle) |
| 10 | same as 9a | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | **lip/xtal/Hz as amounts** | same area law as 9a |
| 11 | same as 9a | 4b lumen chemistry; **[E] on collapse-capped V** | **two-compartment: aq lumen `V_DV(t)`, NLB `V_nlb` const** | same area law as 9a |
| 12a | same as 6 (`f_exp` × host) | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | **aq ⇄ lip** (Kp=398, Myburgh) | **`k_hz × [Fe3]_aq`** (Myburgh eq. 4.53) |
| 12b | **Myburgh empirical `A·B·exp(B·t)`** (host-conserving) | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 12c | Myburgh empirical; **constant `[Hb_RBC]`**; const **1 fL** | **6-parallel peptide MM, constant `[E]`** | constant `[E]` (no `s_PM`) | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 13 | Myburgh empirical (as 12b), host-conserving; **upper-range MCHC budget (~112 fg)** | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 14a | as 13; **`k_release` constant (t½ 20 min), decoupled from `s_PM`** | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 14b | as 14a; **constant `k_release`, t½ 30 min** (slower, provisional) | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 15a | as 13; **`k_release = k_Klemba · [½ + ½ s_PM]`** (constitutive + blot) | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 15b | as 13; Myburgh `v_up` × **Elliott 24–30 h cytostome smoothstep** | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | aq ⇄ lip | `k_hz × [Fe3]_aq` |
| 99 | same as 9a; **slower plateau `k_release` (fit)** | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | same; **`K_xtal` fit** | same area law as 9a |

**Incremental ladder:** one mechanistic change per step.

| Step | Problem in previous model | What this model changes |
|------|---------------------------|-------------------------|
| 1 | Need a minimal closed Fe path | Linear uptake + PMs + FP2/3 + ox + Hz |
| 2a / 2b | Linear uptake leaves most Fe in host | Empirical `f_exp(t)` × remaining host; enzyme **amount** unchanged. **2a** = one phase (comparison). **2b** = two phases, break at Garnie Fig. 5B 29 h, `f_exp` continuous at the join (ladder default). Models 3–9c inherit Model 2b |
| 3 | Full PaxDB amount from `t` = 0 collapses DV Hb | Garnie Fig. 3 `s_PM(t)` on that amount; `f_exp` unchanged |
| 4a / 4b | Peptide `kcat` on native tetramer; PM II dominates `Vmax` | Goldberg ordered pathway: native-competent set vs peptide MM on globin |
| 5 | Lumen `Vmax` still ≫ uptake; standing Hb ~0 | Inaccessible **inner-vesicle** cargo in the DV before 4b lumen proteases; `k_release` from Klemba lysis bound |
| 6 | Constant `k_release` makes assay Hb track `v_up(t)` | `k_release(t) = k_Klemba · s_PM(t)` (Fig. 3 clock; plateau keeps Klemba `t½`) |
| 7 | Single Fe(III) pool drains Hm at lipid-assay `k_hz` | Aqueous ⇄ lipid Fe(III); `v_hz` on lipid pool at full `k_hz` (no `φ`) |
| 8 | Lipid pool still drains Hm at full `k_hz` | Interfacial / crystal-competent Fe(III); `v_hz` on xtal at full `k_hz` |
| 9a/9b/9c | First-order `k_hz` does not grow Egan’s accelerating structure | `v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^α` — α = 2/3 (sphere), 1/2 (rod), 1/3 (needle) |
| 10 | Model 9 Hm/Hb dip at 39h; experiment flat/rising | All Fe species as `AMOUNT_SPECIES`; **no change in fg outputs** (see model10.md) |
| 11 | Is the late dip a lumen-collapse concentration artifact? | Explicit **two-compartment volumes** (aq lumen `V_DV(t)` shrinks, NLB constant), ratio-scaled fluxes; **identical fg curves to 9a** → volume-invariant, the late Hm dip is the `f_exp` taper not geometry (see model11.md) |
| 12a | Model 9 needs an area-growth patch to hold Hm | Crystallise from **aqueous** hematin (Myburgh eq. 4.53), lipid as buffer; drops the `(n/n₀)^α` patch. Best combined fit (see model12.md) |
| 12b | Isolate whether the late Hm/Hb dip is host-depletion | Myburgh's empirical uptake, **host-conserving**. Its magnitude over-delivers vs a conserved 106 fg budget → host exhausts ~43.5 h → honest late cliff |
| 12c | Faithful Myburgh replication | Myburgh's **constant `[Hb_RBC]`** BC + const vol + const `[E]` 6-parallel digestion. Monotonic through 44 h (no cliff), but host+DV Fe not conserved (→194 fg) and Hz overshoots ~5 fg |
| 13 | 12b's 44 h cliff is host exhaustion, not a rate/pool problem | The 106 fg budget is a **population-mean** cell (MCHC 34 g/dL); Garnie's 44 h DV total (~105.5 fg) already equals it yet Hb still turns over → these cells carried **upper-range** Hb. Use MCHC 36 g/dL (~112 fg): host survives, cliff gone, one cited parameter (see model13.md) |
| 14a / 14b | Model 13's early assay-Hb hump (low, noisy early `s_PM` slows release → cargo piles up) | Decouple inner-vesicle lysis from plasmepsin amount: **constant `k_release`** (vesicle rupture isn't protease-gated). 14a t½ 20 min (Klemba) fits late/low early; 14b t½ 30 min fits early/overshoots late → release must accelerate (see model14.md) |
| 15a / 15b | Model 13's *shape* is fine; early Hb *amounts* sit too high. Must not use Garnie Fig. 5B (scoring inventory) | **15a:** constitutive + `s_PM` lysis (`f = ½`); early cargo down, late unchanged, Hb χ² 0.43 — dip damps. **15b:** Elliott 2008 cytostome window gates Myburgh `v_up` (0 before 24 h); starves early Hm/Hz (see model15.md) |
| 99 | What-if: can Hb/Hm match if two knobs are freed? | **Not chemistry** — slower plateau `k_release` and larger `K_xtal` vs Garnie Dd2; inherits Model 9a area; `f_exp` unchanged |

Model 2a’s `f_exp` is a provisional first-order remaining-host clock with `a`, `b` fit to cumulative **fractionation** Fe (`Hb+Hm+Hz`; see [model2.md](models/model2.md) and [garnie_fractionation.md](garnie_fractionation.md)) — not lumen volume, not pHrodo Fig. 2C, not cytostome kinetics, and not “delivery to the parasite” with a later DV step. On Dd2, ODE DV Fe **R² = 0.612**, RMSE 15.89 fg (~36 fg still in host at 44 h). Model 2b is the same constitutive class with a faster late phase after Garnie Fig. 5B’s 29 h break, `f_exp` continuous at the join (DV Fe RMSE 2.56, R² = 0.990, host ~5 fg at 44 h); **Models 3–9c inherit 2b**. The **form** is general; `a` and `b` can be refit to another strain. Variable `V_DV(t)` (`variable_dv_volume`) is **shared bookkeeping** from Model 1 onward (not a numbered ladder step). The present numerical schedule is Garnie Dd2 **aqueous lumen**. Model 3’s one change is the blot-derived amount clock (see [model3.md](models/model3.md)). Model 4’s one change is native vs peptide substrate (two encodings: [model4.md](models/model4.md)). Model 5’s one change is inaccessible **inner-vesicle** cargo already in the DV on the **4b** pathway (see [model5.md](models/model5.md)); 4a’s globin pool is not carried because it never accumulated. A native-Hb `kcat_app` was attempted and still lacks moles of enzyme ([enzyme_kinetics.md](enzyme_kinetics.md)); that gap is not filled with a Garnie-fitted scalar. Model 6’s one change is `k_release(t) ∝ s_PM(t)` so lysis follows the same Fig. 3 clock (see [model6.md](models/model6.md)). Model 7’s one change is Egan lipid partition of Fe(III) (see [model7.md](models/model7.md)), not `φ` on `k_hz`. Model 8’s one change is interfacial (crystal-competent) Fe(III) so bulk NLB haem is not the `k_hz` substrate (see [model8.md](models/model8.md)). Model 9’s one change is crystal-area growth of that `v_hz` (see [model9.md](models/model9.md)): amount, not lumen concentration; `n_Hz_start` is the 20 fg seed; three exponents test crystal habit (9a sphere, 9b rod, 9c needle). Model 10's one change is amount encoding for all Fe species ([model10.md](models/model10.md)): all Fe and Hz become `AMOUNT_SPECIES` but this doesn't change fg outputs because the volume conversion compensates. Model 99 is a **what-if** on Model 9a topology ([model99.md](models/model99.md)), not a cited mechanism: slower plateau `k_release` and slightly larger `K_xtal`, chosen on the **2a** inventory and not refit. It does **not** retune `f_exp`.

`v_dig` includes only proteases that liberate haem from Hb / haem-bearing globin (PMs + falcipains). Downstream peptidases are omitted. Falcipains are present from Model 1 onward (Degradation remains PMs-only as a sandbox).

---

## Fit vs Garnie Dd2 (tracking)

Scores are **diagnostics**, not an objective to minimize with extra terms. A ladder step counts as an improvement only if the change is mechanistic **and** the relevant series move toward the assay.

**Protocol:** `t = [0, 1700]` min from 16 h, `t_eval` step 20 min, init `[0.018, 0, 0, 0.36]`. Model interpolated onto Dd2 ages 20–44 h (`n` = 9). Hm is scored as `conc_fe3pp` (Models 1–6), aq + lip (Model 7), or aq + lip + xtal (Models 8–99). Model 4a scores assay Hb as native + globin; Models 5–99 as inner-vesicle + lumen Hb. Init Hb seeds the inner-vesicle pool in Models 5–99. What the tubes measure: [garnie_fractionation.md](garnie_fractionation.md).

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
| 9a | 0.33 | 0.74 | 1.14 | 53.3 | 2.30 | 0.08 | 2.56 | 0.09 |
| 9b | 0.33 | 0.74 | 1.18 | 59.7 | 2.36 | 0.09 | 2.56 | 0.09 |
| 9c | 0.33 | 0.74 | 1.61 | 69.3 | 2.67 | 0.10 | 2.56 | 0.09 |
| 10 | 0.33 | 0.74 | 1.14 | 53.3 | 2.30 | 0.08 | 2.56 | 0.09 |
| 11 | 0.33 | 0.74 | 1.14 | 53.3 | 2.30 | 0.08 | 2.56 | 0.09 |
| 12a | 0.33 | 0.74 | **0.93** | **4.29** | **2.18** | 0.10 | 2.56 | 0.09 |
| 12b | 0.73 | 4.95 | 0.72 | 3.41 | 5.76 | 0.87 | 6.16 | 1.13 |
| 12c | 0.50 | 2.01 | 0.49 | 6.19 | 6.41 | 1.15 | 6.29 | 1.13 |
| 13 | 0.67 | 4.83 | 0.58 | 3.34 | 5.77 | 0.87 | 6.29 | 1.13 |
| 14a | 0.50 | 2.01 | 0.49 | 6.19 | 6.41 | 1.15 | 6.29 | 1.13 |
| 14b | 0.88 | 3.82 | 0.50 | 4.76 | 5.72 | 0.97 | 6.29 | 1.13 |
| 15a | **0.32** | **0.43** | 0.51 | 5.31 | 6.19 | 1.05 | 6.29 | 1.13 |
| 15b | 0.92 | 9.80 | 1.36 | 210 | 12.51 | 2.63 | 13.22 | 2.91 |
| 99 | 2.33 | 52.7 | 1.24 | 1.29 | 2.82 | 0.11 | 2.56 | 0.09 |

**Reading the ladder:** Model 2a improved Hz and internalized Fe relative to Model 1, but still leaves ~36 fg in the host at 44 h (inventory **R² = 0.612**). Model 2b (ladder default; Models 3–9c inherit it) uses a faster late `f_exp` after Fig. 5B’s 29 h break, still × remaining host, continuous at the join: DV Fe RMSE 15.89 → 2.56 (R² = 0.990), host ~5 fg at 44 h, Hz RMSE 11.63 → 5.32. The last gulp is still a little short. Hb stayed collapsed through Model 4 (lumen `Vmax` ≫ uptake). Shared `variable_dv_volume` did not move M1/M2 fg scores (amount-linear rates). Model 3 matches Model 2b: lag `s_PM(20 h) ≈ 0.41` of plateau is already enough enzyme to collapse DV Hb. Models 4a and 4b also match Model 3. Model 5 is the first step that moves Hb (χ²_red 37 → 4.4): assay Hb tracks inner-vesicle cargo while lumen Hb stays ~0. DV Fe is unchanged (same Model 2b internalized inventory). Model 6 clocks lysis with `s_PM(t)` so standing cargo is `v_up / (k · s_PM)` rather than `v_up / k` ([model6.md](models/model6.md)): Hb χ²_red 4.4 → 0.74. Model 7 splits Fe(III) aq ⇄ lip and crystallizes from the lipid pool at full literature `k_hz` (no `φ`). Hm barely moves (χ²_red 217 → 207). Model 8 puts `k_hz` on interfacial Fe only: Hm signed error flips from −2.88 to +2.60 fg (high, not drained). `K_xtal = 3δ/R = 0.08` from cited NLB radius and film thickness ([model8.md](models/model8.md)) — do not retune it. Model 9 scales that `v_hz` by growing crystal area ([model9.md](models/model9.md)): Hb and DV Fe stay at Model 8; Hm RMSE 3.31 → 1.14 (signed +2.60 → 0.00), Hz RMSE 4.10 → 2.30. Models 9b (1/2) and 9c (1/3) test elongated crystal habit. Model 12a ([model12.md](models/model12.md)) instead crystallises from the **aqueous** pool (Myburgh 2023 eq. 4.53) with the lipid as a buffer, dropping Model 9's area-growth patch: Hm RMSE 1.14 → 0.93 (χ²_red 53.3 → 4.3), Hz RMSE 2.30 → 2.18 (signed −0.04), Hb unchanged — the best combined fit on the ladder, with no retuned constants. Models 12b/12c both swap our `f_exp × host` uptake for Myburgh's empirical exponential but differ in the host boundary condition. **12b keeps our conserving ladder:** the exponential over-delivers vs a conserved 106 fg budget, so the finite host is exhausted near 43.5 h and assay Hb/Hm fall sharply — an honest late cliff, not a bug. **12c uses Myburgh's own constant-`[Hb_RBC]` reservoir** (plus constant 1 fL + constant `[E]` 6-parallel digestion): with the infinite RBC reservoir Myburgh assumes, Hb/Hm/Hz rise smoothly through 44 h with **no cliff** (Hm signed ≈ 0), at the cost that host+DV Fe is not conserved (model total → ~194 fg) and Hz overshoots ~5 fg. 12a's residual dip, by contrast, is not the Hm/Hz mechanism but our `f_exp × host` uptake tapering as the host depletes (leaving ~5 fg at 44 h). The 12a-taper and 12b-cliff both point to the same next step: a host-conserving late-phase delivery schedule (Garnie Fig. 5B phased `v_up`). Model 13 ([model13.md](models/model13.md)) isolates the *cause* of 12b's cliff instead: decomposing 12b shows lumen Hb ≈ 0 (the assay Hb is the HTV cargo at `v_up/(k·s_PM)`) and the 44 h drop is **host exhaustion** — Myburgh's exponential drains the ~85 fg host to zero, `v_up` stops dead, and the standing cargo decays in its ~20 min half-life. The 106 fg budget is a *population-mean* cell (MCHC 34 g/dL × MCV 90 fL), but Garnie's 44 h DV total (~105.5 fg) already equals that whole budget while ~2 fg is still Hb-form and turning over — impossible unless the cells carried upper-range Hb. Model 13 therefore uses the clinical upper MCHC (36 g/dL → ~112 fg): the host survives (~1 fg at 44 h), `v_up` never stops dead, and the cliff vanishes with one cited parameter and no reshaping — keeping 12b's smooth Hm/Hz/Hb (Hz still ~5 fg high, a separate Myburgh-exponential over-delivery). The remaining early Hb hump is the low-early-`s_PM` release clock, decoupled from the terminal behaviour. Models 14a/14b ([model14.md](models/model14.md)) address that hump: since inner-vesicle lysis is a membrane-rupture event (not gated by plasmepsin amount, which only digests Hb *after* it reaches the lumen), release is decoupled from `s_PM` and made constant. 14a (Klemba t½ 20 min) removes the hump and halves Hb χ²_red (4.83 → 2.01) but runs low at 20–29 h (a fast pool can't hold ~2 fg on ~2 fg/h uptake); 14b (t½ 30 min) fits the early points but overshoots to ~4 fg at 44 h. Neither constant rate fits both ends, so inner-vesicle lysis must **accelerate** with development — the direction Model 6's `s_PM` encoded, only with a gentler early magnitude than the low, noisy blot. An ad hoc `s_PM` floor fits best but is rejected as a plot-driven cap. Models 15a/15b ([model15.md](models/model15.md)) treat Model 13's early Hb as an *amount* offset on a good shape, and **do not** use Garnie Fig. 5B (the scoring inventory). **15a** splits lysis into constitutive Klemba rupture plus an `s_PM` term (`f = ½`): early cargo falls (20 h 1.96 → 1.17 fg), late is Model 13, Hb χ²_red 4.83 → 0.43 — the 29 h dip damps because it was the same blot lever. **15b** gates Myburgh `v_up` with Elliott 2008's 24–30 h cytostome window (0 before 24 h): that starves early Hm/Hz (Hm χ²_red 210) — Elliott said the cytostome *increases* then, not that continuous feed is off. Model 99 is a **what-if** on the same ODEs ([model99.md](models/model99.md)): plateau `k_release` and `K_xtal` were chosen on **2a** and are not refit; on 2b plus the Model 6 clock they overshoot Hb (χ²_red 0.74 → 53). Ranked next *ladder* step if 2b’s last ~5 fg is still short: Garnie Fig. 5B Dd2 phases as `v_up` (0.9 then 4.8 fg/h), not a new `f_exp` `b`.

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
- Lumen DV species are integrated in **M at variable `V_DV(t)`** (`variable_dv_volume`) and converted to fg/cell as `C·V_DV(t)`. `constants.vol_dv = 1 fL` is the **reference** volume for the init API and PaxDB amount `n_E = [E]_1fL · 1 fL`. Init `[0.018, 0, 0, 0.36]` stays 1 fL-reference M so the fg seed is unchanged; `run()` rescales lumen species to true M at `V(t=0)`. Model 5–99 inner-vesicle cargo (`conc_hb_htv`) is an **amount** (`AMOUNT_SPECIES`): inside the DV but not pHrodo aqueous lumen, encoded as M at `V_ref`, not diluted by `dV_DV/dt`, converted as `C·V_ref`. The current lumen schedule is Garnie Dd2 (`variable_dv_volume_L`).
- `[O2−]` = 0 (SOD) → Fe(III) reduction is off.
- Parameters: [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py)
- Experiment tables: [`experimental_data.py`](../haem_kinetics/components/experimental_data.py) ([garnie_fractionation.md](garnie_fractionation.md))

### Shared notation

| Symbol | Code / meaning |
|--------|----------------|
| `t` | Time in **minutes** from 16 h post-invasion |
| `[Hb]_HTV` | `conc_hb_htv` — Model 5–99 inner-vesicle cargo already in the DV (amount as M at `V_ref`; not pHrodo lumen) |
| `[Hb]_DV` | `conc_hb_dv` — protease-accessible lumen Hb as haem-equivalents (M); Model 4a: native tetramer; 4b/5–99: lumen Hb |
| `[Hb]_globin` | `conc_hb_globin` — Model 4a nicked globin (haem still protein-bound) |
| `[Hb]_tet` | `[Hb]_DV` / 4 — tetramer basis for MM |
| `[Hb]_RBC` | `conc_hb_rbc` — remaining host Hb (M, RBC basis) |
| `[Fe3]_aq`, `[Fe3]_lip` | Model 7–12 aqueous and bulk-lipid Fe(III) (lumen-basis M). Model 12: Hz forms from `[Fe3]_aq`, lipid is a buffer |
| `[Fe3]_xtal` | Model 8–99 interfacial / crystal-competent Fe(III); assay Hm = aq + lip + xtal |
| `v_up` | Uptake rate into DV (M haem-eq / min on `V_DV`) |
| `v_dig` | Digestion / haem-release rate (M haem-eq / min) |
| `v_ox` | Fe(II)→Fe(III) oxidation rate |
| `v_hz` | Haemozoin formation rate |
| `f_exp(t)` | Empirical `a · b · exp(b · t)` uptake; first-order remaining host. Model 2a = one phase; Model 2b / 3–9c = two-phase at Fig. 5B 29 h. Fit to cumulative **fractionation** Fe, not pHrodo |
| `a_e, b_e, b_l` | Model 2b two-phase `f_exp` (still × remaining host); break at Fig. 5B 29 h; `a_l` from continuity of `f_exp` at the join |
| `k_release` | Model 5 constant inner-vesicle lysis (Klemba `t½` < 20 min bound). Models 6–9c: `k_htv_release · s_PM(t)`. Model 99: that clock with a fitted plateau scale |
| `V_DV(t)` | Variable DV lumen (`variable_dv_volume`); current schedule: Garnie Dd2 Gompertz then linear collapse |
| `s_PM(t)` | Garnie Fig. 3 relative PM amount (Model 3–9c enzyme amount; Model 6–9c also lysis); plateau 40–44 h = 1 |

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
6. Model 2b / 3–9c `f_exp` is empirical (cumulative fractionation Fe, first-order remaining host, two-phase at Fig. 5B 29 h, continuous at the join), not cytostome kinetics and not pHrodo Fig. 2C ([model2.md](models/model2.md), [garnie_fractionation.md](garnie_fractionation.md)). Model 2a is the one-phase comparison. Do **not** retune `a, b` in Model 99.
7. **Late Hz on this ladder is mostly inventory, not crystal lag.** With Model 2b, ~5 fg remains in the host at 44 h (DV Fe RMSE 2.56). Ranked next ladder step if that last gulp is still short: Garnie Fig. 5B Dd2 phases as `v_up` (0.9 then 4.8 fg/h); see [model99.md](models/model99.md#next-mechanisms-ranked-not-this-page).

---

## References

- Assay vs interpretation: [garnie_fractionation.md](garnie_fractionation.md) (saponin pellet, not isolated DVs; pHrodo is a different experiment).
- Garnie LF, Egan TJ, Wicht KJ. Heme processing in the malaria parasite, *Plasmodium falciparum*: a time-dependent basal-level analysis. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z) · [Nature full text](https://www.nature.com/articles/s42003-025-08991-z) · [PDF](https://www.nature.com/articles/s42003-025-08991-z.pdf) · [Figshare raw data](https://doi.org/10.6084/m9.figshare.28801805)
- Combrinck JM, Fong KY, Gibhard L, Smith PJ, Wright DW, Egan TJ. Optimization of a multi-well colorimetric assay to determine haem species in *Plasmodium falciparum*. *Malar. J.* (2015) 14:253. [doi:10.1186/s12936-015-0729-9](https://doi.org/10.1186/s12936-015-0729-9)
- Egan TJ, Chen JY, de Villiers KA, et al. Haemozoin (β-haematin) biomineralization requires both a lipid medium and an accelerating structure to promote haem dimerization. *Malaria Journal* (2012) 11:337. [doi:10.1186/1475-2875-11-337](https://doi.org/10.1186/1475-2875-11-337)
- Enzyme kcat/Km: [`enzyme_kinetics.md`](enzyme_kinetics.md) (Banerjee 2002; Luker 1996; Ramjee 2006 — DOIs listed there).
