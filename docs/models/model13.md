# Model 13

**Code:** [`model13.py`](../../haem_kinetics/models/model13.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 12](model12.md) · **Next:** [Model 14](model14.md)

**Model 12b (Myburgh empirical uptake, host-conserving) on an upper-range red-cell Fe budget.** The single change from Model 12b is the *size of the host Fe pool*: it uses the clinical upper bound of MCHC instead of the population mean. That one cited-parameter change removes 12b's spurious 44 h cliff — with no new rate term, no delivery reshaping, and none of the ugly 29 h step of the earlier Fig. 5B draft.

---

## Why 12b had a cliff, and what actually fixes it

12b keeps the conserving ladder on Myburgh's monotonic exponential uptake. Decomposing its trajectory shows two things:

- Lumen Hb (`conc_hb_dv`) is ≈ 0 at every timepoint — the assay Hb is **entirely the HTV inner-vesicle cargo**, sitting at the quasi-steady value `n_HTV ≈ v_up / (k_release · s_PM)`.
- The 44 h drop is **host exhaustion**: Myburgh's exponential drains the ~85 fg host to exactly zero at ~43.5 h, `v_up` stops dead, and the standing cargo (2.6 fg) decays within its ~20 min release half-life. Reducing the *early* standing pool (e.g. 4× faster phase-1 release) leaves the 44 h collapse byte-for-byte identical; delivering ~12 % less Fe (so the host survives) removes it entirely. The cliff is therefore a host-budget problem, not a rate-law or standing-pool problem.

**How the 106 fg budget was built — from population means.** `Constants.compute_conc_hb_rcb()` uses MCHC = 34 g/dL and MCV = 90 fL, i.e. MCH ≈ 30.6 pg Hb/cell → 106 fg Fe. Both inputs are averages.

**Why the upper end is data-motivated.** Garnie's Dd2 fractionation at 44 h is Hz 97.7 + Hm 5.8 + Hb 2.0 ≈ **105.5 fg** — essentially the *whole* average budget is already inside the DV — yet ~2 fg is still Hb-form and turning over. An average-106-fg cell cannot hold 105.5 fg in the DV **and** keep delivery running; it is fully drained, which is exactly what forces the host to zero and crashes the pool. For the balance to close with ongoing turnover, these cells must have carried **upper-range** haemoglobin.

## The one change (12b → 13)

| Quantity | Model 12b (mean cell) | Model 13 (upper range) | Clinical range |
|--|--:|--:|--:|
| MCHC | 34 g/dL | **36 g/dL** | 32–36 |
| MCV (`vol_rbc`) | 90 fL | 90 fL | 80–100 |
| MCH → Fe budget | 30.6 pg → 106 fg | 32.4 pg → **112 fg** | 27–33 pg |

Implemented by overriding `_full_hb_rbc_m()` (and the matching budget diagnostic) to the upper-MCHC concentration. Everything else — Myburgh's exponential uptake, aqueous-hematin crystallisation, aqueous ⇄ lipid partition, HTV cargo, host conservation — is exactly Model 12b. **Provisional:** MCHC is a per-donor/per-cell distribution; 36 g/dL is a cited upper bound, not a fitted number.

---

## Fit vs Garnie Dd2

Protocol and definitions: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

Mean signed error: Hb +0.52, Hm +0.01, Hz +4.99, DV Fe +5.52.

| Series | RMSE (fg/cell) | MAE | mean signed error | χ²_red | n |
|--------|---------------:|----:|------------------:|-------:|--:|
| Hb | 0.67 | 0.54 | 0.52 | 4.83 | 9 |
| Hm | 0.58 | 0.45 | 0.01 | 3.34 | 9 |
| Hz | 5.77 | 4.99 | 4.99 | 0.87 | 9 |
| DV Fe | 6.29 | 5.52 | 5.52 | 1.13 | 9 |

---

## Known behaviour / findings

1. **The cliff is gone, with no added complexity.** Host is drawn down to ~1 fg at 44 h but never exhausts, so `v_up` never stops dead and the standing HTV cargo holds at ~2.7 fg through 44 h (12b collapsed it to 0.8). Hb, Hm and Hz are smooth and monotonic — the good 12b shape, terminal crash removed. This is a single cited-parameter change (host budget), not a new mechanism.

2. **The residual overshoot is Myburgh's uptake, not the budget.** Myburgh's exponential delivers ~88 fg over 16–44 h regardless of how much host is available, so on the larger budget it fills the DV to ~109 fg vs Garnie's ~105.5, and Hz runs ~5 fg high (signed +4.99; still within the 12–20 fg Hz SEM). That over-delivery is a property of the empirical exponential (same as 12b/12c), separate from the budget question.

3. **The early Hb hump is a separate, smaller issue.** Assay Hb peaks ~2.9 fg near 26 h vs experiment ~1.9 because `n_HTV ≈ v_up/(k_release·s_PM)` and early `s_PM` is low (Garnie's blots: few plasmepsins early), so release lags and cargo piles up. It is decoupled from the terminal behaviour and drives most of the Hb χ².

**Interpretation.** The 44 h cliff was an artifact of scoring an average-106-fg cell against data whose DV inventory already equals that whole average budget. Correcting the host budget to the upper reference bound — the value the data implies — removes the cliff cleanly. What remains (Hz ~5 fg high, the early Hb hump) is not about the budget: the first is Myburgh's exponential over-delivering, the second is the low-early-`s_PM` release clock. The early hump is taken up by [Model 14](model14.md) (decoupling inner-vesicle lysis from the plasmepsin blot).

---

## How to run

```python
from haem_kinetics.models.model13 import Model13

model = Model13()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model13.png',
)
```

---

## References

- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z) — Dd2 fractionation (44 h DV total ≈ 105.5 fg).
- Myburgh KDV. PhD thesis, University of Cape Town (2023), Chapter 4 — Table 4.5 empirical uptake. Extract: [`docs/KDV_2023_Myburgh_model-extract.pdf`](../KDV_2023_Myburgh_model-extract.pdf).
- MedlinePlus. Mean corpuscular hemoglobin concentration (MCHC) reference range. [medlineplus.gov/ency/article/003648.htm](https://medlineplus.gov/ency/article/003648.htm)
