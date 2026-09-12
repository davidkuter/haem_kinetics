# Model 11

**Code:** [`model11.py`](../../haem_kinetics/models/model11.py)  
**Up:** [Model index](../models.md) · **Prev:** [Model 10](model10.md) · **Next:** [Model 12](model12.md)

Model 9a chemistry (crystal-area growth, sphere 2/3) with an explicit **two-compartment volume model**: the aqueous DV lumen shrinks on the Garnie schedule while the NLB/pigment compartment stays constant. Where Model 10 encoded the lipid/crystal/Hz species as bare amounts, Model 11 gives them their own physical volume and converts inter-compartment fluxes by the volume ratio.

**Result:** identical fg time courses to Model 9a (and Model 10). Splitting the volumes does **not** change first-order mass profiles — the same volume-invariance conclusion as Model 10, reached a second way.

---

## The change (9a → 11)

| Compartment | Volume | Species | Dilution |
|-------------|--------|---------|----------|
| aqueous lumen | `V_aq = V_DV(t)` (Garnie; collapses 36–46 h) | `conc_hb_dv`, `conc_fe2pp` | yes (`−C·V̇/V`) |
| NLB / pigment | `V_nlb = 1 fL` (constant) | `conc_fe3pp_aq`*, `conc_fe3pp_lip`, `conc_fe3pp_xtal`, `conc_hz` | no |

*Fe(III)_aq is placed on the NLB basis because it equilibrates rapidly with the lipid phase at the NLB–water interface. Oxidation flux Fe(II)→Fe(III)_aq crosses the lumen→NLB boundary and is scaled by `V_aq/V_nlb`. Enzyme concentration uses a collapse-capped volume (`max(V_DV, 1 fL)`) since the proteases are membrane-associated and don't concentrate as the lumen shrinks.

Mechanistic basis: NLB lipid droplets don't shrink with the aqueous lumen (Jackson 2004); Hz crystals are solid and excluded from the pHrodo volume (Garnie 2025); the aq⇄lip exchange is an interfacial process.

---

## Fit vs Garnie Dd2

Protocol: [models.md](../models.md#fit-vs-garnie-dd2-tracking). Recompute with `python examples/run.py`.

| Series | RMSE (fg/cell) | MAE | mean signed | χ²_red | n |
|--------|---------------:|----:|------------:|-------:|--:|
| Hb | 0.33 | 0.25 | 0.14 | 0.74 | 9 |
| Hm | 1.14 | 0.88 | 0.00 | 53.28 | 9 |
| Hz | 2.30 | 1.92 | −0.35 | 0.08 | 9 |
| DV Fe | 2.56 | 1.98 | −0.21 | 0.09 | 9 |

Identical to Model 9a and Model 10.

---

## Known behaviour / findings

1. **Volume compartmentalization is mass-invariant.** Giving the NLB/pigment species their own constant volume, with ratio-scaled inter-compartment fluxes, reproduces Model 9a's fg curves exactly. Model 10 showed this by encoding those species as amounts; Model 11 shows the same by giving them a physical constant volume. For first-order kinetics the lumen-collapse bookkeeping cancels either way (Myburgh's own Model 4 reached the same conclusion).

2. **The late Hb *and* Hm decrease (38–44 h) is inherited from Model 9a's uptake, not the volume model.** Both pools peak near 38 h then fall (Hb 2.74 → 1.5 fg; Hm 5.1 → 3.0 fg) because `f_exp × host` uptake **tapers** as the host depletes (33 fg at 38 h → 5 fg at 44 h): the Fe feed drops while Hz crystallisation keeps draining free haem. 
   - The **Hb** decrease is plausibly *real* — Garnie's Dd2 Hb dips 2.46 → 1.98 fg over 41 → 44 h as the parasite finishes internalised Hb.
   - The **Hm** decrease is *wrong* — Garnie's Hm keeps rising (5.17 → 5.83 fg). Free haem should accumulate, not drain, as Hz outpaces the tapering feed.

3. **The fix lives in the uptake law, not the geometry.** The late Hm drain is the same `f_exp × host` taper that Models 12a and 13 addressed on the uptake side: [Model 13](model13.md)'s non-tapering delivery (Myburgh uptake on an upper-range host budget) sustains the feed and gives a *rising* late Hm (5.66 fg at 44 h). So Model 11's late-Hm defect is not a compartmentalization problem — the two-compartment split is correct but inert for mass, and the real lever is keeping the Fe feed alive late.

---

## How to run

```python
from haem_kinetics.models.model11 import Model11

model = Model11()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model11.png',
)
```

---

## References

- Jackson KE, et al. Food vacuole-associated lipid bodies / neutral lipid in *Plasmodium falciparum*. (NLB persistence during lumen collapse.)
- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z) — pHrodo lumen volume excludes solid pigment.
