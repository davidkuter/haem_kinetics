# Model 10 (what-if)

**Code:** [`haem_kinetics/models/model10.py`](../../haem_kinetics/models/model10.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 9](model9.md)

**This is not a mechanistic ladder step.** Same ODEs as Model 9 (including crystal-area growth). Two knobs are freed against Garnie Dd2 to ask: *if* inner-vesicle lysis were slower and *if* the interfacial fraction were slightly larger, could Hb and Hm approach the assay?

They are **not** Klemba’s `t½ < 20 min` and **not** `K_xtal = 3δ/R`. **`f_exp` is unchanged** (Model 2b two-phase): do not retune `a, b` here. The `k_release` / `K_xtal` numbers below were chosen on the **Model 2a** inventory and are not refit. Model 8 keeps the citations; Model 9 keeps the area law. This page is a diagnostic.

---

## What changed vs Model 9

| Item | Model 9 | Model 10 |
|------|---------|----------|
| ODEs | aq ⇄ lip ⇄ xtal → Hz with area factor; inner-vesicle cargo | **Unchanged** |
| `f_exp` | Model 2b two-phase | **Unchanged** |
| `s_PM` on lysis | `k_release(t) = k_Klemba · s_PM(t)` (Model 6) | **Unchanged form** |
| plateau `k_release` | ln(2)/20 min⁻¹ (Klemba bound) | **0.477 × that** (plateau `t½ ≈ 41.9 min`) |
| `K_xtal` | `3δ/R = 0.08` | **0.10** (grid vs Hm χ²_red) |
| Area factor | `(n_Hz / n_Hz_start)^{2/3}` | **Inherited** (not refit) |

So `k_10(t) = 0.4774 · k_Klemba · s_PM(t)`. The scale was not refit after Model 6 or Model 9.

**Not claimed:** a new trafficking paper, a new NLB size, or that Garnie prefers these numbers as biology. `k_hz` is still the Egan lipid-assay value.

---

## How the knobs were chosen

Standing inner-vesicle cargo is QSS: `n_HTV ≈ v_up / k_release(t)`. Scaling the plateau `k` by `s` is `k_10 = s · k_Klemba · s_PM(t)`. Least-squares `s` vs Garnie Dd2 Hb (ages 20–44 h) is **0.4774** → plateau `t½ ≈ 41.9 min`. That is slower than Klemba’s **upper** bound on half-time (`t½ < 20 min`). Klemba’s figure is PM II trafficking, not inner-vesicle Hb lysis; this scale is still a diagnostic, not a replacement citation. It was chosen on **Model 2a** and is not refit.

With that `k_release` held, `K_xtal` was gridded `{0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.24}`:

| `K_xtal` | Hm RMSE | Hm χ²_red | Hm signed |
|--------:|--------:|----------:|----------:|
| 0.08 (M8 geometry) | 0.69 | 15.8 | +0.19 |
| **0.10** | 0.85 | **2.4** | −0.43 |
| 0.12 | 1.18 | 3.5 | −0.85 |

**0.10** is the lowest Hm χ²_red on that grid (RMSE prefers 0.08; early points have small SEM). `k_xtal_ex` is unchanged (does not set standing Hm).

Starting the ODE at 0 h post-invasion would **not** replace this `k_release` tweak: standing Hb is `v_up(now)/k`, not integrated ring-stage uptake ([timing note](../models.md)).

---

## Why the what-if knobs miss on 2b (do not retune)

`k_release` and `K_xtal` were least-squares / grid chosen on the interfacial model **with Model 2a uptake**. Models 3–10 now inherit 2b, Models 6–10 inherit `k_release ∝ s_PM`, and Model 10 also inherits Model 9 area growth. Late `v_up` is larger, so standing cargo overshoots Dd2 Hb. **Do not refit these two knobs** to recover the old χ² — that would be a second Garnie fit, not chemistry.

Hz is no longer an inventory problem: DV Fe RMSE is 2.56 (host ~5 fg at 44 h). The leftover ~5 fg is the last gulp 2b still damps.

---

## What this can and cannot do

- **Could (on 2a, without area):** show the interfacial topology had enough freedom for Hb χ²_red ~1 and Hm χ²_red ~2 if `k_release` and `K_xtal` were free.
- **Cannot (on 2b, with those same knobs):** recover that Hb/Hm match. Late delivery is larger; the 2a-tuned scale overshoots.
- **Cannot** flatten Hb vs age by a constant scale. Model 6 already clocks lysis with `s_PM`; this page only scales that clock.
- **Cannot** restore Klemba or NLB geometry. A later mechanistic step would be a cited lysis time or a cited `(δ, R)`, not these two numbers.

---

## Next mechanisms (ranked; not this page)

One mechanistic change per numbered model. Do not invent a new `a, b`.

1. **Recommended if Model 2b’s last gulp is still short:** Garnie Fig. 5B Dd2 delivery phases as **`v_up`**, not a new `f_exp` `b` — **0.9 fg/h** (20–29 h) then **4.8 fg/h** (29–44 h), clipped by remaining host. Model 3 deferred those numbers as `v_dig` because they are the scoring inventory. Using them as uptake is the observation `f_exp` claimed to encode but does not (wrong *shape* for the last ~5 fg). **Not** pHrodo Fig. 2C (standing acidic-lumen probe, not `dFe/dt`).
2. **Later, if (1) is too circular:** cytostome / DV-surface-limited fusion (Garnie’s `r²` vs `r³` argument). Rate set by contact area, not `f(t) × n_host`. Needs a cited area or vesicle-fusion law.
3. **Hb without slowing `k_release`:** a cited inner-vesicle lysis time (Klemba is PM II trafficking), or a measured native-Hb `kcat` with moles of enzyme ([enzyme_kinetics.md](../enzyme_kinetics.md)). No vesicle-number cap. Do not refit the 2a what-if scale on 2b. Model 6 already tests `k_release ∝ s_PM`.
4. **Hm without fitting `K_xtal`:** keep `3δ/R = 0.08`. Model 9 already tests `v_hz ∝ A_Hz`. Extra delivered Fe may still overfill assay Hm; if the *shape* is still wrong, a cited NLB size/count vs time is the next chemistry, not a fitted exponent. Do not raise `k_hz` or restore `φ`.
5. **Do not bother for late Hz while ~5 fg stays in the host:** lipid split, xtal fraction, and `k_hz` cannot supply that iron.

---

## Fit vs Garnie Dd2

Protocol: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

| Series | RMSE (fg/cell) | MAE | mean signed error | χ²_red | n |
|--------|---------------:|----:|-----:|-------:|--:|
| Hb | 2.33 | 2.26 | +2.26 | 52.7 | 9 |
| Hm | 1.24 | 0.75 | −0.63 | 1.29 | 9 |
| Hz | 2.82 | 2.58 | −1.84 | 0.11 | 9 |
| DV Fe | 2.56 | 1.98 | −0.21 | 0.09 | 9 |

**Vs Model 9:** the 2a-tuned knobs overshoot Hb (χ²_red 0.74 → 52.7, signed +2.26). Hm χ²_red 53 → 1.29 because this page inherits the area law and a larger `K_xtal`; that is not a reason to treat 0.10 as chemistry. Hz RMSE 2.30 → 2.82. DV Fe unchanged. Do not refit `k_release` / `K_xtal` here.

---

## How to run

```python
from haem_kinetics.models.model10 import Model10

model = Model10()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model10.png',
)
```

---

## References

- Klemba bound (what Model 10 **violates** at plateau): [model5.md](model5.md)
- Model 6 `s_PM` lysis clock (what Model 10 **scales**): [model6.md](model6.md)
- Geometric `K_xtal` (what Model 10 **replaces**): [model8.md](model8.md)
- Crystal-area growth (what Model 10 **inherits**): [model9.md](model9.md)
- Model 2b `f_exp` (what Model 10 **does not retune**): [model2.md](model2.md)
- Assay vs pHrodo: [garnie_fractionation.md](../garnie_fractionation.md)
- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z)
