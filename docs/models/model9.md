# Model 9a / 9b / 9c

**Code:** [`model9a.py`](../../haem_kinetics/models/model9a.py), [`model9b.py`](../../haem_kinetics/models/model9b.py), [`model9c.py`](../../haem_kinetics/models/model9c.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 8](model8.md) · **Next:** [Model 10](model10.md)

Model 8 chemistry (HTV cargo, `k_release(t) ∝ s_PM`, 4b lumen proteases, aq ⇄ lip ⇄ xtal) plus **crystal-area growth**. Haemozoin still forms from the interfacial pool at literature `k_hz`, now scaled by growing Hz surface. Three variants explore crystal-habit geometry:

| Variant | Exponent | Geometry | Rate |
|---------|----------|----------|------|
| **9a** | 2/3 | sphere | `v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{2/3}` |
| **9b** | 1/2 | rod/needle | `v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{1/2}` |
| **9c** | 1/3 | extreme elongation | `v_hz = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^{1/3}` |

`n_Hz = [Hz] · V_DV(t)` (amount). Lumen concentration is not used: dilution must not fake a change in crystal size. `n_Hz_start` is the protocol seed, `0.36` M at `V_ref` ≈ **20 fg**. At `t = 0` the area factor is 1 (same rate as Model 8). The exponent is geometry, not a Garnie fit.

---

## What changed vs Model 8 (and why)

**Problem in Model 8:** `v_hz = k_hz · [Fe3]_xtal` is first-order in the interface only. Assay Hm sits high (signed **+2.60 fg**) and Hz low (signed **−2.95 fg**): the same iron is parked in non-crystalline Fe(III) instead of crystal. Egan's lipid-mediated β-haematin needs a lipid medium **and** an accelerating structure. A first-order `k_hz` does not grow that structure.

**Change (one mechanism):** crystallization rate scales with crystal surface area.

| Item | Model 8 | Model 9 |
|------|---------|---------|
| Fe(III) | aq ⇄ lip ⇄ xtal → Hz | Unchanged topology |
| `v_hz` | `k_hz · [Fe3]_xtal` | **`k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^α`** |
| `k_hz`, `K_xtal` | literature / `3δ/R` | **Unchanged** |
| `f_exp` | Model 2b two-phase | Unchanged |
| Assay Hm | aq + lip + xtal | Unchanged |

**Not this step:**

- Fitting `k_hz`, `K_xtal`, or the exponent to Garnie Hm.
- Changing `f_exp`.
- Restoring `φ`.
- The what-if plateau `k_release` scale and `K_xtal = 0.10` ([Model 99](model99.md)).

---

## Crystal habit and exponent

β-haematin crystals are **elongated faceted prisms**, not spheres. The exponent `α` depends on crystal habit:

| Geometry | Area scaling | Exponent |
|----------|--------------|----------|
| Sphere | A ∝ V^{2/3} | **2/3** |
| Rod (fixed diameter) | A ∝ length ∝ V^{1/2} | **1/2** |
| Needle (length only) | A ∝ V^{1/3} | **1/3** |

Kapishnikov TEM and Egan in vitro studies show β-haematin is elongated. The 2/3 exponent (sphere) is the geometric default; 1/2 and 1/3 test whether elongated-crystal geometry better matches the late Hm/Hz drain.

A lower exponent means **slower acceleration** late → Hm drains more slowly → higher Hm at 44 h.

---

## Process schematic

```mermaid
flowchart LR
  Host["conc_hb_rbc"] -->|"f_exp"| HTV["conc_hb_htv"]
  HTV -->|"k_release"| Lumen["conc_hb_dv"]
  Lumen -->|"PM I/II/FP-2"| Fe2["conc_fe2pp"]
  Fe2 -->|"k_ox x O2"| Fe3aq["conc_fe3pp_aq"]
  Fe3aq <--> Fe3lip["conc_fe3pp_lip"]
  Fe3lip <--> Fe3xtal["conc_fe3pp_xtal"]
  Fe3xtal -->|"k_hz x A_Hz"| Hz["conc_hz"]
```

Assay **Hb**: `conc_hb_htv + conc_hb_dv`. Assay **Hm**: `conc_fe3pp_aq + conc_fe3pp_lip + conc_fe3pp_xtal`.

---

## Governing equations

HTV, lumen digestion, oxidation, aq ⇄ lip ⇄ xtal: [model8.md](model8.md). Model 9 replaces only `v_hz`:

```text
n_Hz       = [Hz] · V_DV(t)
n_Hz_start = 0.36 · V_ref
v_hz       = k_hz · [Fe3]_xtal · (n_Hz / n_Hz_start)^α

d [Fe3]_xtal / dt = v_lip↔xtal − v_hz + dil([Fe3]_xtal)
d [Hz] / dt       = v_hz + dil([Hz])
```

`Hz_start` is the existing init Hz (~20 fg), not a fitted nucleation seed. Early half-time cost is none: the factor is 1 at start.

---

## Parameters (Model 9–specific)

| Constant | Value | Units | Source |
|----------|------:|-------|--------|
| `n_Hz_start` | `0.36 · V_ref` | mol | Protocol init Hz (~20 fg); not a fit |
| α (9a) | 2/3 | — | Sphere geometry |
| α (9b) | 1/2 | — | Rod/needle geometry |
| α (9c) | 1/3 | — | Extreme elongation |
| `k_hz` | 0.12 | min⁻¹ | Egan lipid-mediated β-haematin (unchanged) |
| `K_xtal` | 0.08 | — | `3δ/R` from Model 8 (unchanged) |

---

## Assumptions

- Growing Hz is treated as a geometric body for the area–amount relation. Real β-haematin is faceted; exponents are idealized limits.
- Amount (`C · V_DV(t)`), not lumen concentration, is the size proxy.
- The 20 fg seed is already present at the 16 h offset (same init as every other model).

---

## Known behaviour / issues

- Hb and DV Fe match Model 8 (uptake and lysis unchanged).
- Hm and Hz move because `v_hz` grows with crystal area.
- Lower exponent (9b, 9c) → slower late acceleration → higher Hm at 44 h.
- At ~90 fg Hz: 9a factor ≈ 2.7; 9b ≈ 2.1; 9c ≈ 1.7.

---

## Fit vs Garnie Dd2

Protocol and definitions: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

| Model | RMSE Hb | χ² Hb | RMSE Hm | χ² Hm | RMSE Hz | χ² Hz | RMSE DV Fe |
|-------|--------:|------:|--------:|------:|--------:|------:|-----------:|
| 9a | 0.33 | 0.74 | 1.14 | 53.3 | 2.30 | 0.08 | 2.56 |
| 9b | 0.33 | 0.74 | 1.18 | 59.7 | 2.36 | 0.09 | 2.56 |
| 9c | 0.33 | 0.74 | 1.61 | 69.3 | 2.67 | 0.10 | 2.56 |

**Vs Model 8:** Hb and DV Fe match. 9a Hm RMSE 3.31 → 1.14, signed +2.60 → 0.00, χ²_red 102 → 53. Hz RMSE 4.10 → 2.30. 9b and 9c test whether elongated geometry improves the late-drain mismatch.

---

## How to run

```python
from haem_kinetics.models.model9a import Model9a
from haem_kinetics.models.model9b import Model9b
from haem_kinetics.models.model9c import Model9c

model = Model9a()  # or Model9b(), Model9c()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model9a.png',
)
```

---

## References

- Egan accelerating structure (what Model 8 deferred): [model8.md](model8.md)
- Egan TJ, Chen JY, de Villiers KA, et al. Haemozoin (β-haematin) biomineralization requires both a lipid medium and an accelerating structure to promote haem dimerization. *Malaria Journal* (2012) 11:337. [doi:10.1186/1475-2875-11-337](https://doi.org/10.1186/1475-2875-11-337)
- Kapishnikov S, Weiner A, Shimoni E, et al. Oriented nucleation of hemozoin at the digestive vacuole membrane in *Plasmodium falciparum*. *Proc. Natl. Acad. Sci. USA* (2012) 109:11188–11193. [doi:10.1073/pnas.1118120109](https://doi.org/10.1073/pnas.1118120109)
- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z)
