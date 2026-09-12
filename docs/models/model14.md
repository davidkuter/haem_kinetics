# Model 14 (14a / 14b)

**Code:** [`model14a.py`](../../haem_kinetics/models/model14a.py) · [`model14b.py`](../../haem_kinetics/models/model14b.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 13](model13.md) · **Next:** [Model 15](model15.md)

**Addressing Model 13's early assay-Hb hump by decoupling inner-vesicle release from plasmepsin amount.** Both variants make the same mechanistic change from Model 13 — release is no longer `∝ s_PM(t)` — and differ only in the lysis timescale. Together they bracket the constant-rate hypothesis and show that inner-vesicle lysis must *accelerate* through development.

---

## Why Model 13 had an early hump

Assay Hb is the HTV cargo (lumen Hb ≈ 0) at the quasi-steady value `n_HTV ≈ v_up / (k_release · s_PM)`. Model 6 set `k_release(t) = k_Klemba · s_PM(t)` — flagged even there as "a modelling convenience, not a claim." Garnie's early plasmepsin blot is both **low** (~0.41 of plateau) and **noisy** (a dropped-outlier dip at 24 h), so the release rate collapses exactly where the data is flat and the cargo piles into a spurious 24–26 h hump (~2.9 fg vs experiment ~1.9). Denoising the blot (monotone `s_PM`) barely helps — it is the whole early level, not just the 24 h point.

## The one mechanistic change (13 → 14)

Inner-vesicle membrane lysis (delivering Hb to the DV lumen) is a fusion/rupture event; the plasmepsins digest Hb only *after* it reaches the lumen, so they do not gate vesicle lysis. Model 14 therefore **decouples release from `s_PM`**. Digestion still scales with enzyme amount (`s_PM`); only release is decoupled.

| | Release `k_release` | Rationale |
|--|---------------------|-----------|
| Model 13 | `k_Klemba · s_PM(t)` | Model 6 convenience coupling |
| **14a** | constant `k_Klemba` (t½ ≈ 20 min) | Klemba cytostomal-delivery bound, taken as the lysis rate |
| **14b** | constant, t½ = 30 min (~1.5× Klemba) | Klemba's <20 min bounds PM *trafficking*, not necessarily vesicle-membrane lysis — provisional slower timescale |

---

## Fit vs Garnie Dd2

Protocol and definitions: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

**Model 14a** — constant Klemba rate (mean signed: Hb −0.25, Hm −0.02, Hz +5.79, DV Fe +5.52):

| Series | RMSE | MAE | signed | χ²_red | n |
|--------|-----:|----:|-------:|-------:|--:|
| Hb | 0.50 | 0.42 | −0.25 | 2.01 | 9 |
| Hm | 0.49 | 0.40 | −0.02 | 6.19 | 9 |
| Hz | 6.41 | 5.79 | 5.79 | 1.15 | 9 |
| DV Fe | 6.29 | 5.52 | 5.52 | 1.13 | 9 |

**Model 14b** — slower constant t½ = 30 min (mean signed: Hb +0.53, Hm −0.06, Hz +5.05, DV Fe +5.52):

| Series | RMSE | MAE | signed | χ²_red | n |
|--------|-----:|----:|-------:|-------:|--:|
| Hb | 0.88 | 0.64 | 0.53 | 3.82 | 9 |
| Hm | 0.50 | 0.40 | −0.06 | 4.76 | 9 |
| Hz | 5.72 | 5.05 | 5.05 | 0.97 | 9 |
| DV Fe | 6.29 | 5.52 | 5.52 | 1.13 | 9 |

(Assay Hb at the Dd2 points, fg/cell — 20 → 44 h. Experiment: 1.22, 1.44, 2.11, 1.65, 1.68, 2.16, 2.17, 2.46, 1.98.)

| | 20h | 26h | 32h | 38h | 44h |
|--|--:|--:|--:|--:|--:|
| Model 13 (`∝ s_PM`) | 1.96 | **2.94** | 2.13 | 2.21 | 2.73 |
| 14a (t½ 20) | 0.83 | 1.12 | 1.51 | 2.03 | 2.74 |
| 14b (t½ 30) | 1.23 | 1.66 | 2.24 | 3.01 | **4.06** |

---

## Known behaviour / findings

1. **Decoupling removes the hump (14a).** With release at the constant Klemba rate, assay Hb rises smoothly and monotonically — no 24–26 h pile-up. Hb χ²_red falls 4.83 → 2.01 and Hb/Hm RMSE both improve. Digestion still scales with `s_PM`, so Hm/Hz are essentially unchanged.

2. **A fast pool runs low early (14a).** At t½ ≤ 20 min the protected pool cannot retain the observed ~2 fg on the low early uptake (~2 fg/h), so 14a sits below the 20–29 h points. This is an honest residual about the lysis timescale, not a defect to patch.

3. **A slower constant rate fixes early but overshoots late (14b).** t½ = 30 min matches the early points (20 h ≈ experiment) but, since `n_HTV ≈ v_up / k_release` and late uptake is large, the standing pool climbs to ~4 fg at 44 h (experiment ~2). So Hb χ²_red is *worse* (3.82) despite the better early fit.

4. **No single constant rate fits both ends → release accelerates.** The 14a↔14b bracket (20 min fits late/misses early; 30 min fits early/misses late) shows inner-vesicle lysis must speed up through development. That is exactly the *direction* Model 6's `s_PM` coupling encoded; only its early magnitude — from the low, noisy blot — was too extreme.

**Interpretation.** The hump was a release-clock artifact, and decoupling from the plasmepsin blot is the correct mechanistic move. But the early-vs-late bracket says release is genuinely developmentally accelerated, not constant. The next mechanistic step is a *cited* maturation schedule for vesicle lysis / cytostome fusion (DV-surface or membrane-turnover kinetics), rather than either a constant rate or the noisy PM blot as a proxy. What is **not** the answer: an ad hoc floor on `s_PM` (fits best, χ²_Hb ~1.5, but it is a cap chosen to move the plot).

---

## How to run

```python
from haem_kinetics.models.model14a import Model14a  # or Model14b

model = Model14a()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model14a.png',
)
```

---

## References

- Klemba M, Beatty W, Gluzman I, Goldberg DE. Trafficking of plasmepsin II to the food vacuole. *J. Cell Biol.* (2004) 164:47–56. [doi:10.1083/jcb.200307147](https://doi.org/10.1083/jcb.200307147) — cytostomal-delivery t½ < 20 min (PM trafficking bound).
- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z) — Fig. 3 plasmepsin blot (`s_PM`), Dd2 fractionation.
