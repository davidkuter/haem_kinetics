# Haem kinetics models

Simulates haemoglobin (Hb) uptake into the *Plasmodium falciparum* digestive vacuole (DV), enzymatic release of haem, Fe(II)→Fe(III) oxidation, and detoxification to haemozoin (Hz). Experimental targets are Garnie et al. (*Commun. Biol.* 2025) heme fractionation time courses (fg Fe/cell), primarily **Dd2**.

**Time base:** simulation `t` is in **minutes** from a trophozoite offset of **16 h** post-invasion (`parasite age = 16 + t/60` hours).

## Modeling principles

Every addition to the ladder must be **mechanistically accountable** (chemistry / cell biology / cited measurement). We do **not** add numerical patches, efficiency fudges, or ad hoc rate caps to improve plots or stabilize the solver. Solver tolerances exist only to integrate the **stated** ODEs accurately. If the solution goes negative or total Fe drifts, that is a signal to fix the model or the integration settings — not to clip or reshape the RHS.

---

## Model pages (active ladder)

| Model | File | One-line summary |
|-------|------|------------------|
| [Degradation](models/degradation.md) | `degradation.py` | Uptake / Fe(II) release sandbox |
| [Model 1](models/model1.md) | `model1.py` | Linear uptake; PMs + FP2/3; full Fe speciation; no lipid |
| [Model 2](models/model2.md) | `model2.py` | Model 1 + empirical `f_exp` uptake (PaxDB amount unchanged) |
| [Model 3](models/model3.md) | `model3.py` | Model 2 + Garnie Fig. 3 `s_PM(t)` amount schedule |
| [Model 4a / 4b](models/model4.md) | `model4a.py`, `model4b.py` | Native-Hb-competent enzymes (ordered pathway vs single pool) |
| [Model 5](models/model5.md) | `model5.py` | Model 4b + inaccessible HTV / inner-vesicle Hb cargo |
| [Model 6](models/model6.md) | `model6.py` | Model 5 + aqueous ⇄ lipid Fe(III); Hz from lipid at `k_hz` |

**Archived prior ladder** (φ → f_exp+[E] → lipid → logistic → xtal): [`docs/models/legacy/`](models/legacy/README.md) and `haem_kinetics/models/legacy/`.

---

## Quick comparison

| Model | Hb uptake | Proteases | Enzyme schedule | Lipid / Fe(III) | Hz rate |
|-------|-----------|-----------|-----------------|-----------------|---------|
| Degradation | Exponential | PMs | `f_exp` | — | — |
| 1 | Linear | **PMs + FP2/3** | `n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 2 | empirical `f_exp(t)` | PMs + FP2/3 | `n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 3 | empirical `f_exp(t)` | PMs + FP2/3 | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 4a | empirical `f_exp(t)` | **PM I/II/FP-2 on native**; all six peptide MM on globin | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 4b | empirical `f_exp(t)` | **PM I/II/FP-2 native rate only** (lumped) | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 5 | `f_exp` into **HTV cargo**; first-order release to lumen | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | None | `k_hz × [Fe3]` |
| 6 | same as 5 | 4b lumen chemistry | `s_PM(t) · n_E / V_DV(t)` | **aq ⇄ lip** | `k_hz × [Fe3]_lip` |

**Incremental ladder:** one mechanistic change per step.

| Step | Problem in previous model | What this model changes |
|------|---------------------------|-------------------------|
| 1 | Need a minimal closed Fe path | Linear uptake + PMs + FP2/3 + ox + Hz |
| 2 | Linear uptake leaves most Fe in host | Empirical `f_exp` uptake; enzyme **amount** unchanged |
| 3 | Full PaxDB amount from `t` = 0 collapses DV Hb | Garnie Fig. 3 `s_PM(t)` on that amount; `f_exp` unchanged |
| 4a / 4b | Peptide `kcat` on native tetramer; PM II dominates `Vmax` | Goldberg ordered pathway: native-competent set vs peptide MM on globin |
| 5 | Lumen `Vmax` still ≫ uptake; standing Hb ~0 | Inaccessible HTV cargo before 4b lumen proteases; `k_release` from Klemba vesicle→lumen bound |
| 6 | Single Fe(III) pool drains Hm at lipid-assay `k_hz` | Aqueous ⇄ lipid Fe(III); `v_hz` on lipid pool at full `k_hz` (no `φ`) |

Model 2’s `f_exp` is a provisional Fe-delivery schedule fit to cumulative DV Fe (see [model2.md](models/model2.md)) — not lumen volume, not cytostome/HTV kinetics, and not an independent prediction of total Fe. Variable `V_DV(t)` (`variable_dv_volume`) is **shared bookkeeping** from Model 1 onward (not a numbered ladder step). The present numerical schedule is Garnie Dd2 lumen. Model 3’s one change is the blot-derived amount clock (see [model3.md](models/model3.md)). Model 4’s one change is native vs peptide substrate (two encodings: [model4.md](models/model4.md)). Model 5’s one change is inaccessible HTV cargo on the **4b** pathway (see [model5.md](models/model5.md)); 4a’s globin pool is not carried because it never accumulated. A native-Hb `kcat_app` was attempted and still lacks moles of enzyme ([enzyme_kinetics.md](enzyme_kinetics.md)); that gap is not filled with a Garnie-fitted scalar. Model 6’s one change is Egan lipid partition of Fe(III) (see [model6.md](models/model6.md)), not `φ` on `k_hz`.

`v_dig` includes only proteases that liberate haem from Hb / haem-bearing globin (PMs + falcipains). Downstream peptidases are omitted. Falcipains are present from Model 1 onward (Degradation remains PMs-only as a sandbox).

---

## Fit vs Garnie Dd2 (tracking)

Scores are **diagnostics**, not an objective to minimize with extra terms. A ladder step counts as an improvement only if the change is mechanistic **and** the relevant series move toward the assay.

**Protocol:** `t = [0, 1700]` min from 16 h, `t_eval` step 20 min, init `[0.018, 0, 0, 0.36]`. Model interpolated onto Dd2 ages 20–44 h (`n` = 9). Hm is scored as `conc_fe3pp` (Models 1–5) or aq + lip (Model 6). Model 4a scores assay Hb as native + globin; Models 5–6 as HTV + lumen Hb. Init Hb seeds HTV in Models 5–6.

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
| 2 | 1.91 | 37.13 | 3.38 | 231 | 11.63 | 0.57 | 15.89 | 1.17 |
| 3 | 1.91 | 37.13 | 3.38 | 231 | 11.63 | 0.57 | 15.89 | 1.17 |
| 4a | 1.91 | 37.13 | 3.38 | 231 | 11.63 | 0.57 | 15.89 | 1.17 |
| 4b | 1.91 | 37.13 | 3.38 | 231 | 11.63 | 0.57 | 15.89 | 1.17 |
| 5 | 1.03 | 10.43 | 3.38 | 232 | 12.27 | 0.62 | 15.89 | 1.17 |
| 6 | 1.03 | 10.43 | 3.34 | 225 | 12.30 | 0.62 | 15.89 | 1.17 |

**Reading the ladder:** Model 2 improved Hz and internalized Fe (`DV_Fe` = model DV pools vs assay Hb+Hm+Hz) relative to Model 1. Hb stayed collapsed through Model 4 (lumen `Vmax` ≫ uptake). Hm is still far outside SEM. Shared `variable_dv_volume` did not move M1/M2 fg scores (amount-linear rates). Model 3 matches Model 2: lag `s_PM(20 h) ≈ 0.41` of plateau is already enough enzyme to collapse DV Hb. Models 4a and 4b also match Model 3. Model 5 is the first step that moves Hb (χ²_red 37 → 10): assay Hb tracks HTV cargo (~0.5–1.2 fg, rising with `f_exp`) while lumen Hb stays ~0. Hz is slightly delayed vs 4b (vesicle→lumen lag). DV Fe is unchanged (same internalized inventory). Model 6 splits Fe(III) aq ⇄ lip and crystallizes from the lipid pool at full literature `k_hz` (no `φ`). Hm barely moves (χ²_red 232 → 225): lipid Fe³⁺ still drains into Hz. That is the accountable result of this chemistry, not a reason to add `xtal` or slow `k_hz` in this step. Hb and DV Fe match Model 5; Hz is slightly more delayed.

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
- Lumen DV species are integrated in **M at variable `V_DV(t)`** (`variable_dv_volume`) and converted to fg/cell as `C·V_DV(t)`. `constants.vol_dv = 1 fL` is the **reference** volume for the init API and PaxDB amount `n_E = [E]_1fL · 1 fL`. Init `[0.018, 0, 0, 0.36]` stays 1 fL-reference M so the fg seed is unchanged; `run()` rescales lumen species to true M at `V(t=0)`. Model 5–6 HTV cargo is an **amount** (`AMOUNT_SPECIES`): encoded as M at `V_ref`, not diluted by `dV_DV/dt`, converted as `C·V_ref`. The current lumen schedule is Garnie Dd2 (`variable_dv_volume_L`).
- `[O2−]` = 0 (SOD) → Fe(III) reduction is off.
- Parameters: [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py)
- Experiment tables: [`experimental_data.py`](../haem_kinetics/components/experimental_data.py)

### Shared notation

| Symbol | Code / meaning |
|--------|----------------|
| `t` | Time in **minutes** from 16 h post-invasion |
| `[Hb]_HTV` | `conc_hb_htv` — Model 5–6 inaccessible cargo (amount as M at `V_ref`) |
| `[Hb]_DV` | `conc_hb_dv` — DV Hb as haem-equivalents (M); Model 4a: native tetramer; 4b/5/6: lumen Hb |
| `[Hb]_globin` | `conc_hb_globin` — Model 4a nicked globin (haem still protein-bound) |
| `[Hb]_tet` | `[Hb]_DV` / 4 — tetramer basis for MM |
| `[Hb]_RBC` | `conc_hb_rbc` — remaining host Hb (M, RBC basis) |
| `[Fe3]_aq`, `[Fe3]_lip` | Model 6 aqueous and lipid-associated Fe(III) (lumen-basis M); assay Hm = aq + lip |
| `v_up` | Uptake rate into DV (M haem-eq / min on `V_DV`) |
| `v_dig` | Digestion / haem-release rate (M haem-eq / min) |
| `v_ox` | Fe(II)→Fe(III) oxidation rate |
| `v_hz` | Haemozoin formation rate |
| `f_exp(t)` | Empirical `a · b · exp(b · t)` uptake (Model 2–6); fit to cumulative DV Fe |
| `k_release` | Model 5–6 first-order HTV → lumen (Klemba vesicle→lumen `t½` = 20 min bound; provisional) |
| `V_DV(t)` | Variable DV lumen (`variable_dv_volume`); current schedule: Garnie Dd2 Gompertz then linear collapse |
| `s_PM(t)` | Garnie Fig. 3 relative PM amount (Model 3–6); plateau 40–44 h = 1 |

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

where `[E]_i,eff` is `n_E,i / V_DV(t)` on Models 1–2 and `s_PM(t) · n_E,i / V_DV(t)` on Model 3–6. Degradation’s sandbox still uses `f_exp`-scaled PMs. Models 1–3 sum all six peptide terms on DV Hb. Model 4 uses `k_enzymes_native` on the tetramer (PM I, PM II, FP-2 only); 4a applies the peptide table to nicked globin; 4b, 5, and 6 lump haem release with that native rate.

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

1. **`variable_dv_volume` is shared bookkeeping**, not a ladder win: Garnie et al. 2025 measure a dynamic Dd2 lumen (~3.7 fL peak), which currently supplies `V_DV(t)`. Amount-linear rates make fg scores insensitive to that volume. Uptake is not `C·dV/dt` and not assay `dF/dt`.
2. **Standing DV Hb stayed collapsed through Model 4** because lumen `Vmax` ≫ uptake. Model 5 puts internalized Hb in inaccessible HTV cargo first ([model5.md](models/model5.md)); lumen native still collapses. A literature native-Hb `kcat_app` was attempted and still lacks moles of enzyme ([enzyme_kinetics.md](enzyme_kinetics.md)). Do not invent that number from Garnie fg or park haem. Model 6 addresses Hm via lipid partition ([model6.md](models/model6.md)), not by slowing `k_hz` with `φ`.
3. **NF54 digits** in `experimental_data.py` may not match Garnie 2025; prefer Dd2. `s_PM` itself is from NF54 Fig. 3 blots (the published PM time course).
4. Peptide `kcat`/`Km` (Banerjee/Luker; Ramjee) applied to DV Hb is an approximation — see [`enzyme_kinetics.md`](enzyme_kinetics.md).
5. ppm are Tao 2014 Dd2 whole-organism (not DV-specific).
6. Model 2–6 `f_exp` is empirical (cumulative DV Fe fit), not cytostome/HTV kinetics; see [model2.md](models/model2.md).

---

## References

- Garnie LF, Egan TJ, Wicht KJ. Heme processing in the malaria parasite, *Plasmodium falciparum*: a time-dependent basal-level analysis. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z) · [Nature full text](https://www.nature.com/articles/s42003-025-08991-z) · [PDF](https://www.nature.com/articles/s42003-025-08991-z.pdf) · [Figshare raw data](https://doi.org/10.6084/m9.figshare.28801805)
- Egan TJ, Chen JY, de Villiers KA, et al. Haemozoin (β-haematin) biomineralization requires both a lipid medium and an accelerating structure to promote haem dimerization. *Malaria Journal* (2012) 11:337. [doi:10.1186/1475-2875-11-337](https://doi.org/10.1186/1475-2875-11-337)
- Enzyme kcat/Km: [`enzyme_kinetics.md`](enzyme_kinetics.md) (Banerjee 2002; Luker 1996; Ramjee 2006 — DOIs listed there).
