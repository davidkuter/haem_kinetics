# Model 6

**Code:** [`haem_kinetics/models/model6.py`](../../haem_kinetics/models/model6.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 5](model5.md) · **Next:** [Model 7](model7.md)

Model 5 chemistry (inner-vesicle cargo, 4b lumen proteases, `f_exp`, `k_hz · [Fe3]`) plus **lysis that follows the Garnie Fig. 3 plasmepsin amount clock**. Inner-vesicle breakdown is first-order in cargo as in Model 5, but `k_release` is no longer constant.

---

## What changed vs Model 5 (and why)

**Problem in Model 5:** assay Hb is QSS inner-vesicle cargo, `n_HTV ≈ v_up / k_release`, with a **constant** `k_release` from Klemba’s `t½ < 20 min` bound. Continuous 2b uptake is slow early and faster after 29 h, so the Hb trace copies that shape (low through the 20s, then a late rise). Dd2 Hb is flatter (~1.2–2.5 fg).

Klemba’s figure is PM II **arriving** in the DV. Garnie Fig. 3 already supplies a cited relative amount clock `s_PM(t)` (used for enzyme amount since Model 3). Tying lysis to that clock is the same trafficking biology, not a fit to Dd2 Hb.

**Change (one mechanism):**

```text
k_release(t) = k_htv_release · s_PM(t)
```

`k_htv_release` stays `ln(2)/20` min⁻¹. `s_PM` is [`garnie_pm_amount_scale`](../../haem_kinetics/models/helpers.py) (NF54 Fig. 3; plateau 40–44 h = 1). Early `s_PM(20 h) ≈ 0.41` → `t½ ≈ 49 min`. Late `s_PM → 1` → `t½ = 20 min` (Klemba bound).

| Item | Model 5 | Model 6 |
|------|---------|---------|
| Uptake | Model 2b `f_exp` | Unchanged |
| Lumen proteases | 4b, `s_PM(t) · n_E` | Unchanged |
| `k_release` | `ln(2)/20` (constant) | **`k_htv_release · s_PM(t)`** |
| Fe³⁺ / Hz | `k_hz · [Fe3]` | Unchanged |
| Assay Hb | inner vesicle + lumen | Unchanged |

**Not this step:**

- A scale fitted so Dd2 Hb sits at 1.9 fg (that is [Model 10](model10.md)).
- Flattening 2b to decorate this curve.
- A second unfused extra-DV pool.
- Lipid / interfacial Fe(III) ([Model 7](model7.md), [Model 8](model8.md)).

---

## Process schematic

```mermaid
flowchart LR
  Host["conc_hb_rbc"] -->|"f_exp"| InnerVesicle["conc_hb_htv"]
  InnerVesicle -->|"k_release s_PM"| Lumen["conc_hb_dv"]
  Lumen -->|"PM I/II/FP-2"| Fe2["conc_fe2pp"]
  Fe2 -->|"k_ox x O2"| Fe3["conc_fe3pp"]
  Fe3 -->|"k_hz"| Hz["conc_hz"]
```

---

## Governing equations

Same states and ODEs as [Model 5](model5.md), with

```text
k_release(t) = k_htv_release · s_PM(t)

d C_htv / dt   = v_up,mol / V_ref − k_release(t) · C_htv
d [Hb]_lumen/dt = k_release(t) · n_HTV / V_DV(t) − v_dig + dil([Hb]_lumen)
```

QSS: `n_HTV ≈ v_up / (k_htv_release · s_PM)`. Early cargo is larger than Model 5 at the same `v_up`; late `k` is the Klemba bound.

---

## Parameters (Model 6–specific)

| Constant | Value | Units | Source |
|----------|------:|-------|--------|
| `k_htv_release` | ln(2)/20 ≈ 0.0347 | min⁻¹ | Same Klemba plateau as Model 5; **not** refit |
| `s_PM(t)` | Garnie Fig. 3 | — | Same blot clock as Model 3 enzyme amount ([model3.md](model3.md)) |

**Early `t½` cost:** at 20 h, `t½ ≈ 20 / 0.41 ≈ 49 min`, slower than Klemba’s *upper* bound on half-time. That is an accountable consequence of using the blot as a relative clock with the bound at plateau, not a second fitted delay.

---

## Assumptions

- Inner-membrane lysis scales with the amount of DV-resident PM (Klemba: PM II trafficking into the FV).
- NF54 Fig. 3 is the relative clock on Dd2 models (same approximation as Model 3).
- Lumen `Vmax` remains ≫ uptake, so assay Hb is still the inner-vesicle pool.

---

## Known behaviour / issues

- This lifts early Hb (~1.60 fg at 20 h vs Model 5 ~0.67; assay ~1.2) and χ²_red 4.40 → 0.74. The late 2b peak is still there (~2.7 fg near 40 h). It is not fully flat: `s_PM` rises ~2.4× while late `v_up` rises more than that.
- Hm remains drained by `k_hz · [Fe3]` — [Model 7](model7.md) is the next accountable step (lipid partition, not `φ`).
- Success for this step is the cited clock on lysis, not a still-off Hb χ² used to justify a fitted scale.

---

## Fit vs Garnie Dd2

Protocol and definitions: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

| Series | RMSE (fg/cell) | MAE | mean signed error | χ²_red | n |
|--------|---------------:|----:|-----:|-------:|--:|
| Hb | 0.33 | 0.25 | +0.14 | 0.74 | 9 |
| Hm | 3.20 | 2.95 | −2.95 | 217 | 9 |
| Hz | 3.59 | 2.60 | +2.60 | 0.22 | 9 |
| DV Fe | 2.56 | 1.98 | −0.21 | 0.09 | 9 |

**Vs Model 5:** Hb moves (RMSE 0.65 → 0.33, χ²_red 4.40 → 0.74; signed −0.45 → +0.14). Early cargo is higher because `s_PM` is below 1. DV Fe is unchanged (same 2b inventory). Hm is unchanged (still drained). Hz is slightly less high (RMSE 3.97 → 3.59) because cargo spends longer in the vesicle while `s_PM` is low.

---

## How to run

```python
from haem_kinetics.models.model6 import Model6

model = Model6()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model6.png',
)
```

---

## References

- Klemba M, Beatty W, Gluzman I, Goldberg DE. Trafficking of plasmepsin II to the food vacuole of the malaria parasite *Plasmodium falciparum*. *J. Cell Biol.* (2004) 164:47–56. [doi:10.1083/jcb.200307147](https://doi.org/10.1083/jcb.200307147)
- Garnie Fig. 3 `s_PM(t)`: [model3.md](model3.md)
- Model 5 inner-vesicle topology: [model5.md](model5.md)
- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z)
