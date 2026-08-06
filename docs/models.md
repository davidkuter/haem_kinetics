# Haem kinetics models

This package simulates haemoglobin (Hb) uptake into the *Plasmodium falciparum* digestive vacuole (DV), enzymatic release of haem, Fe(II)→Fe(III) oxidation, and detoxification to haemozoin (Hz). Experimental targets are Combrink-style heme fractionation time courses (fg Fe/cell for residual Hb, free haem, and Hz), primarily **Dd2**.

Simulation time `t` is in **minutes** from a trophozoite offset of **16 h** post-invasion (`parasite age = 16 + t/60` hours), unless a model states otherwise.

---

## Quick comparison

| Model | Hb uptake | Proteases | Lipid / Fe(III) | Hz rate | States | Notes |
|-------|-----------|-----------|-----------------|---------|--------|-------|
| **Degradation** | Exponential (sandbox) | PMs (+ fudge) | None | None (Hz=0) | Hb, Fe2 | Uptake/degradation only |
| **1** | Linear | PMs; fudge **÷** [E] (default 1) | None | `k_hz × [Fe3]` | Hb, Fe2, Fe3, Hz | Baseline |
| **2** | Linear | PMs; fudge **÷** [E] (default 1) | `φ` multiplies Fe3 rates | `k_hz × φ × [Fe3]` | same as 1 | First lipid attempt |
| **3** | Exponential `f_exp(t)` | PMs grow with `f_exp`; fudge **×** [E] (default 1) | same φ penalty | same as 2 | same as 1 | Best legacy visual fit |
| **4** | Same as 3 | Same as 3 | Same as 3 | Same as 3 | same as 1 | Currently identical to Model 3 |
| **5** | Same as 3 | Same as 3 | **Aqueous ⇄ lipid Fe3 pools** | `k_hz × [Fe3]_lipid` (**no φ**) | + Fe3_aq, Fe3_lip | Fixes lipid–Hz logic |
| **6** | Combrink-like sigmoidal; **depleting** host Fe | Same PMs; [E] tracks lumen growth | Same pools as 5 | Same as 5 | **fg/cell** + host Fe | Volume + mass balance |
| **7** | Same as 6 | PMs **+ falcipain-2/3** | Same as 5 | Same as 5 | same as 6 | Falcipains added; fudge unused |

Legacy Models 1–4 keep concentrations in **M** (DV volume basis) and convert to fg/cell for plotting. Models 6–7 integrate **fg Fe/cell** directly.

---

## Shared biology (all full models)

```text
Hb (RBC) --uptake--> Hb (DV) --proteases--> Fe(II)PPIX --O2--> Fe(III)PPIX --> Hz
```

- Superoxide is taken as **0** (SOD), so Fe(III)→Fe(II) reduction is inactive.
- Fe(II) oxidation is fast on the trophozoite timescale (`k_fe2pp_ox × [O2]`).
- Parameters live in [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py).
- Dd2/NF54 tables: [`experimental_data.py`](../haem_kinetics/components/experimental_data.py).

---

## Existing models

### Degradation (`degradation.py`)

Sandbox for Hb **uptake** and enzymatic release of Fe(II). No Fe(III)/Hz ODEs; after integration `conc_hz` is set to 0. Useful for testing transport laws. Docstrings note that early full models could not explain basal free haem (Combrink), which motivated Model 2’s lipid term.

### Model 1 (`model1.py`)

- Linear uptake: `k_hb_trans × [Hb]_RBC`.
- All four DV plasmepsins (PM1, PM2, HAP, PM4) with Michaelis–Menten kinetics.
- `fudge` **lowers** effective enzyme concentration.
- No lipid partitioning; Hz is first-order in total Fe(III).

### Model 2 (`model2.py`)

Same transport and proteases as Model 1, plus a lipid sequestration factor:

```text
φ = (1 − f_lip) / (1 + f_lip + f_lip × K_partition) ≈ 0.133
```

Fe(III)-involving rates (reduction, Hz) are multiplied by `φ`. Intent: leave a basal free-haem pool. Docstring says “HAP-limited”; the code currently sums **all** plasmepsins (HAP-only loop is commented out).

**Limitation:** `k_hz` was taken from lipid-mediated β-haematin assays, then slowed by `φ`—chemically inconsistent with lipids *catalysing* Hz formation.

### Model 3 (`model3.py`)

- Uptake and enzyme levels scale with fractional exponential growth  
  `f_exp(t) = a·b·e^(b t)` (`a=0.1578`, `b=0.001102` for t0 = 16 h).
- `fudge` **raises** effective [E].
- Same `φ` penalty on Hz as Model 2.
- Host `[Hb]_RBC` is reduced once for initial DV contents, then **not** depleted during integration.

This is the main legacy comparison model against Dd2 ([`examples/model3-1.png`](../examples/model3-1.png)): late Hz and free haem tend to under-predict; DV Hb can run slightly high.

### Model 4 (`model4.py`)

Intended as the next WIP after Model 3. In the current tree it is **functionally the same** as Model 3 (exponential uptake, lipid `φ` on Hz). Prefer Model 5+ for new work.

---

## New models (recommendation ladder)

### Model 5 (`model5.py`) — corrected lipid / Hz mechanism

Builds on Model 3’s uptake and protease schedule, but replaces the single Fe(III) pool + `φ` rate penalty with:

| Species | Role |
|---------|------|
| `conc_hb_dv` | Hb in DV (haem-equivalents, M) |
| `conc_fe2pp` | Fe(II)PPIX |
| `conc_fe3pp_aq` | Aqueous Fe(III)PPIX |
| `conc_fe3pp_lip` | Lipid-associated non-Hz Fe(III)PPIX |
| `conc_hz` | Haemozoin |

- Exchange: aqueous ⇄ lipid with effective equilibrium ratio from `K_partition` and `vol_fract_lip`.
- **Hz forms only from the lipid pool:** `d[Hz]/dt = k_hz × [Fe3]_lip` (literature `k_hz`, **not** multiplied by `φ`).
- Plotted “free haem” = `conc_fe3pp_aq + conc_fe3pp_lip` (matches fractionation: non-Hb, non-Hz Fe).

Still uses fixed DV volume and non-dynamic host Hb (same limitations as Model 3 for total Fe budget).

### Model 6 (`model6.py`) — volume + mass balance

Extends Model 5’s speciation with:

1. **States in fg Fe/cell** (extensive), so changing volume does not invent/destroy Fe.
2. **Time-dependent DV lumen volume** `V_DV(t)` (Gompertz rise toward ~3.7 fL, then late collapse)—shape inspired by Combrink et al. 2025; parameters are approximate and documented in code.
3. **Host Hb Fe depleted** as uptake proceeds (hard mass balance toward ~106 fg/cell total Fe).
4. **Sigmoidal uptake rate** of remaining host Fe (Combrink-like Dd2 cargo delivery), instead of Model 3’s pure exponential in concentration space.

Michaelis–Menten steps convert fg→M using `V_DV(t)` when computing rates. Enzyme abundance still follows a lumen-linked growth factor (no falcipains yet; `fudge` retained but should become unnecessary in Model 7).

### Model 7 (`model7.py`) — falcipains

Same transport, volume, and Fe(III)/Hz scheme as Model 6, but Hb degradation includes **falcipain-2 and falcipain-3** in addition to the four plasmepsins. Provisional kcat/Km and abundances are in `Constants` and should be refined against literature / proteomics. Default `fudge` for this model is **1.0** (identity).

---

## How to run

```python
from haem_kinetics.models.model5 import Model5
# Model6, Model7 similarly

model = Model5()
model.run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.0, 0.36],  # Model5: 5 states (M)
    t_eval=range(0, 1700, 20),
    plot='model5.png',
)
```

Model 6/7 initial conditions are **fg/cell** (and include remaining host Hb Fe as the last state)—see each module’s docstring and [`examples/run.py`](../examples/run.py).

---

## Design notes / known issues in legacy models

1. **`φ` on `k_hz` (Models 2–4):** lipids mediate β-haematin; slowing an already lipid-derived `k_hz` by `φ` is the wrong sign for chemistry (fixed in Model 5+).
2. **`vol_dv = 1 fL` vs comment “4 fL”:** inconsistent; real lumen volume is dynamic (~3.7 fL peak in Dd2). Addressed in Model 6+.
3. **Missing falcipains (Models 1–6):** major Hb degraders; formerly motivated `fudge` ≠ 1 (now default **1.0**).
4. **NF54 table in `experimental_data.py`:** legacy digits (Hb starting ~26 fg) do not match Combrink 2025 NF54 (~0.3–1.3 fg). Prefer Dd2 for calibration until NF54 is updated.
5. **No formal optimizer:** comparison is visual overlay in `KineticsModel._plot`.

## Expected behaviour of Models 5–7 (first pass)

These models implement the *mechanistic* fixes; they are **not** yet re-tuned to overlay Combrink curves.

- With literature `k_hz = 0.12 min⁻¹` on the lipid pool and fast aqueous⇄lipid exchange, **free haem often undershoots** (~6 fg target) because almost all Fe(III) crystallises. Next tuning knobs: slower exchange, a non-crystallisable lipid-associated sub-pool, or L/H-dependent `k_hz`.
- Model 6/7 conserve total Fe (~106 fg). Sigmoidal uptake + strong protease activity can drive **DV Hb near zero** late; ease uptake/`kcat` or raise basal Hb if matching the ~2 fg Dd2 residual.
- Falcipain kcat/Km/ppm in Model 7 are **provisional placeholders**.

---

## References (entry points)

- Combrink et al., *Communications Biology* (2025): DV volume, uptake, basal Hb/Hm/Hz time courses.
- Egan et al., *Malaria Journal* 11:337 (2012): lipid-mediated β-haematin kinetics (`k_hz`).
- Plasmepsin / falcipain Hb digestion reviews (Goldberg lab and others): redundant aspartic + cysteine protease pathways.
