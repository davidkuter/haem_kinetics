# Model 15 (15a / 15b)

**Code:** [`model15a.py`](../../haem_kinetics/models/model15a.py) · [`model15b.py`](../../haem_kinetics/models/model15b.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 14](model14.md) · **Next:** [Model 99](model99.md) (what-if)

**Lower Model 13's early assay-Hb *amounts* without treating the trajectory as a hump to flatten, and without using Garnie Fig. 5B (the scoring inventory).** Two independent hypotheses from Model 13:

| | Lever | Source |
|--|--------|--------|
| **15a** | Softer lysis clock: constitutive rupture + `s_PM`-modulated term | Klemba bound vs blot-as-rate |
| **15b** | Gate Myburgh `v_up` by cytostome maturity | Elliott *PNAS* 2008 (24–30 h), **not** Garnie Fe |

---

## Why not Fig. 5B

Garnie Fig. 5B (0.9 then 4.8 fg/h; 20 / 29 h knots) **is** the Dd2 Fe inventory this ladder is scored against. Using those numbers as `v_up` or as a maturation gate is fitting the model to the data it is judged on. Model 15 does not do that.

---

## 15a — two-component release

Model 6 set `k_release = k_Klemba · s_PM(t)`. Early `s_PM ≈ 0.41` then makes `t½ ≈ 49 min`, slower than Klemba's `t½ < 20 min` bound, and the standing cargo `n_HTV ≈ v_up / (k · s_PM)` sits too high. The blot clock still gives the useful wiggle; only the early *scale* is wrong.

15a splits lysis:

```text
k_release(t) = k_Klemba · [f + (1 − f) · s_PM(t)]    with f = 1/2
```

Constitutive inner-membrane rupture (Klemba trafficking) does not wait on the NF54 blot; a protease-modulated term still follows `s_PM`, so late (`s_PM → 1`) is unchanged. Early effective `t½ ≈ 28 min` instead of 49 min. **Provisional:** `f = 1/2` is an equal-weight split, not fitted to Garnie Hb. The strict Klemba reading (`t½` never > 20 min) forces `f = 1` and collapses to Model 14a.

---

## 15b — Elliott cytostome window

Elliott et al. (*PNAS* 2008) 105:2463: rings take up Hb mainly as a one-shot **Big Gulp**; the cytostome's small-vesicle path only “increases its contribution to total hemoglobin uptake” once the parasite is a trophozoite (**24–30 h**). Myburgh's exponential is a continuous cytostome-like feed from invasion — the wrong early process.

15b multiplies Model 13's Myburgh rate by a Hermite smoothstep that is **0 before 24 h** and **1 after 30 h**. No free amplitude. Release stays `∝ s_PM`. This is Elliott's window, not Garnie 29 h.

---

## Fit vs Garnie Dd2

Protocol: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

**Model 15a** (mean signed: Hb −0.01, Hm −0.01, Hz +5.54, DV Fe +5.52):

| Series | RMSE | MAE | signed | χ²_red | n |
|--------|-----:|----:|-------:|-------:|--:|
| Hb | 0.32 | 0.21 | −0.01 | 0.43 | 9 |
| Hm | 0.51 | 0.41 | −0.01 | 5.31 | 9 |
| Hz | 6.19 | 5.54 | 5.54 | 1.05 | 9 |
| DV Fe | 6.29 | 5.52 | 5.52 | 1.13 | 9 |

**Model 15b** (mean signed: Hb −0.35, Hm −0.78, Hz −11.65, DV Fe −12.78):

| Series | RMSE | MAE | signed | χ²_red | n |
|--------|-----:|----:|-------:|-------:|--:|
| Hb | 0.92 | 0.66 | −0.35 | 9.80 | 9 |
| Hm | 1.36 | 0.98 | −0.78 | 210 | 9 |
| Hz | 12.51 | 11.65 | −11.65 | 2.63 | 9 |
| DV Fe | 13.22 | 12.78 | −12.78 | 2.91 | 9 |

Assay Hb (fg/cell) at the Dd2 points:

| | 20h | 26h | 29h | 35h | 44h |
|--|--:|--:|--:|--:|--:|
| Experiment | 1.22 | 2.11 | 1.65 | 2.16 | 1.98 |
| Model 13 | 1.96 | 2.94 | 2.49 | 2.12 | 2.73 |
| **15a** | **1.17** | **1.61** | **1.69** | **1.92** | 2.74 |
| **15b** | 0.03 | 0.30 | 1.82 | 2.12 | 2.73 |

---

## Known behaviour / findings

1. **15a does the job asked of it.** Early amounts come down (20 h 1.96 → 1.17; 26 h 2.94 → 1.61) and sit on the data; late Hb/Hm/Hz are almost Model 13 (same uptake, `s_PM → 1`). Hb χ²_red 4.83 → **0.43**. Hm stays a rising curve (signed ≈ 0). Conserved, no cliff.

2. **Softening the clock also damps the 29 h dip.** Experiment dips 2.11 → 1.65 at 29 h; Model 13 had that wiggle (too high). 15a's Hb is a smooth rise. The early overshoot and the blot-driven dip were the same lever (`s_PM` in the denominator). You cannot lower early amounts a lot via this clock and keep the dip at full strength.

3. **15b is the honest Elliott gate, and it fails.** Zero continuous uptake until 24 h starves HTV, Hm and Hz through the early trophozoite (Hm χ²_red 210; host still 24 fg at 44 h). After 30 h it is Model 13 again. Elliott said the cytostome *increases* its contribution in that window, not that continuous feed is off until 24 h (rings already have the Big Gulp). A hard 0→1 switch on Myburgh's law is too crude — and we will not put a fitted floor under it to rescue the plot.

4. **Hz ~5 fg high in 15a is still Myburgh's exponential**, not the release clock (same DV Fe as Model 13).

**Interpretation.** Early assay Hb was an *amount* problem on Model 13's shape, not a reason to flatten release (14) or to copy Garnie Fig. 5B. A two-component lysis clock (15a) lowers that amount and is the better of these two tests. The remaining 29 h dip in the data, if we want it back, needs a cited process that is *not* the same `s_PM` we just softened — not a Garnie-inventory ramp.

---

## How to run

```python
from haem_kinetics.models.model15a import Model15a  # or Model15b

model = Model15a()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model15a.png',
)
```

---

## References

- Klemba M, Beatty W, Gluzman I, Goldberg DE. *J. Cell Biol.* (2004) 164:47–56. [doi:10.1083/jcb.200307147](https://doi.org/10.1083/jcb.200307147) — cytostomal-delivery `t½` < 20 min.
- Elliott DA, et al. Four distinct pathways of hemoglobin uptake in *P. falciparum*. *PNAS* (2008) 105:2463–2468. [doi:10.1073/pnas.0711067105](https://doi.org/10.1073/pnas.0711067105) — Big Gulp vs cytostome SHV; 24–30 h trophozoite window.
- Gligorijevic B, et al. *Biochemistry* (2006) — Hb digestion / Hz onset ~18–20 h (not used as a knot here).
