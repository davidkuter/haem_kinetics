# Legacy Model 3

**Code:** [`haem_kinetics/models/legacy/model3.py`](../../../haem_kinetics/models/legacy/model3.py)  
**Up:** [Legacy index](README.md) · [Active models](../../models.md) · **Prev:** [Legacy Model 2](model2.md) · **Next:** [Legacy Model 4](model4.md)

Replaces linear uptake with **exponential fractional growth** of remaining host Hb, and scales haem-releasing protease levels (PMs + FP2/3) with the same `f_exp(t)`. Keeps Model 2’s `φ` multiplier on Fe(III) rates.

---

## What changed vs Model 2 (and why)

**Problems in Models 1–2:**

1. **Uptake timing.** Linear `k_hb_trans · [Hb]_RBC` does not reproduce the accelerating Hb delivery through mid–late trophozoite that Combrink-style curves imply. Host→DV cargo should track parasite growth / DV expansion, not a constant first-order coefficient.
2. **Protease timing.** With PaxDB [E] present at full strength from `t` = 0, Models 1–2 digest DV Hb almost as soon as it arrives (Hb→0). Real DV proteases accumulate / mature over the trophozoite window, so early digestion capacity should be low.

**Changes (coupled on purpose — one growth clock):**

| Item | Model 2 | Model 3 |
|------|---------|---------|
| Uptake | `k_hb_trans · [Hb]_RBC` | `f_exp(t) · [Hb]_RBC · V_RBC / V_DV` |
| Effective [E] | Constant PaxDB [E] | `f_exp(t) × [E]_pax` |
| Fe³⁺ / Hz | `φ` on rates | Unchanged |

**Why share `f_exp` for uptake and enzymes:** a single exponential fraction of remaining host Fe is a compact way to gate both cargo delivery and proteolytic capacity to the same developmental schedule. It is a modelling convenience, not a claim that protease expression literally equals the uptake prefactor.

---

## Process schematic

```mermaid
flowchart LR
  Host["conc_hb_rbc"] -->|"f_exp(t) x host"| HbDV["conc_hb_dv"]
  HbDV -->|"PMs+FP2/3 x f_exp(t)"| Fe2["conc_fe2pp"]
  Fe2 -->|"k_ox x O2"| Fe3["conc_fe3pp"]
  Fe3 -->|"k_hz x phi"| Hz["conc_hz"]
```

with `a` = 0.1578, `b` = 0.001102 (`t0` = 16 h).

---

## State variables

| Symbol | Meaning |
|--------|---------|
| `conc_hb_dv` | DV haemoglobin (haem-eq, M) |
| `conc_fe2pp` | Fe(II)PPIX |
| `conc_fe3pp` | Fe(III)PPIX (single pool) |
| `conc_hz` | Haemozoin |
| `conc_hb_rbc` | Remaining host Hb (appended by `run()`) |

Init: `[Hb_DV, Fe2, Fe3, Hz]`.

---

## Governing equations

```text
# Fractional exponential growth (uptake / enzyme clock)
f_exp(t) = a · b · exp(b · t)     a = 0.1578,  b = 0.001102

# Host → DV Hb uptake (M/min on V_DV)
v_up = f_exp(t) · [Hb]_RBC · V_RBC / V_DV

# Effective protease concentration (PMs + falcipains)
[E]_i,eff = f_exp(t) · [E]_i
i ∈ {plm_1, plm_2, hap, plm_4, fp_2, fp_3}

# Haem release from Hb (MM sum; 4 haem-eq per tetramer)
v_dig = 4 · Σ_i  (60 · kcat_i) · [E]_i,eff · [Hb]_tet
                     / (Km_i + [Hb]_tet)

# Fe(II) → Fe(III) oxidation
v_ox  = k_fe2_ox · [Fe(II)] · [O2]

# Aqueous fraction of Fe(III) (lipid sequestration factor)
φ = (1 − f_lip) / (1 + f_lip + f_lip · K_partition)

# Fe(III) → Fe(II) reduction (off: [O2−] = 0; φ-scaled)
v_red = k_fe3_red · φ · [Fe(III)] · [O2−]

# Haemozoin formation (φ-scaled)
v_hz  = k_hz · φ · [Fe(III)]
```

ODEs:

```text
# DV haemoglobin
d[Hb]_DV / dt   = v_up − v_dig

# Free Fe(II)PPIX
d[Fe(II)] / dt  = v_dig + v_red − v_ox

# Free Fe(III)PPIX
d[Fe(III)] / dt = v_ox − v_red − v_hz

# Haemozoin
d[Hz] / dt      = v_hz

# Remaining host RBC Hb
d[Hb]_RBC / dt  = − v_up · V_DV / V_RBC
```

---

## Constants used

| Constant | Value | Units | Description |
|----------|------:|-------|-------------|
| `a` | 0.1578 | — | Prefactor in fractional exponential growth `f_exp(t)` |
| `b` | 0.001102 | min⁻¹ | Rate constant in fractional exponential growth `f_exp(t)` |
| `V_RBC` | 90×10⁻¹⁵ | L | Volume of the host red blood cell |
| `V_DV` | 1×10⁻¹⁵ | L | Fixed digestive-vacuole volume used for M ↔ fg conversion |
| `N_A` | 6.022×10²³ | mol⁻¹ | Avogadro's number |
| `N_prot` | 1.9×10⁸ | — | Average number of proteins per *P. falciparum* cell |
| `k_fe2_ox` | 193800 | min⁻¹ | Rate constant for Fe(II)PPIX oxidation by O₂ |
| `[O2]` | 1×10⁻³ | M | Dissolved oxygen concentration in the DV |
| `k_fe3_red` | 180×10⁻⁹ | — | Rate constant for Fe(III)PPIX reduction by O₂⁻ (inactive when `[O2−]` = 0) |
| `[O2−]` | 0 | M | Superoxide concentration (taken as zero due to SOD) |
| `k_hz` | 0.12 | min⁻¹ | First-order rate constant for haemozoin formation from Fe(III) |
| `f_lip` | 0.016 | — | Fractional volume of lipid nanospheres relative to the DV |
| `K_partition` | 398 | — | Equilibrium partition coefficient of Fe(III)PPIX into lipid |
| `φ` | ≈ 0.133 | — | Aqueous fraction of Fe(III) at partition equilibrium; multiplies Fe(III) rates |

Enzyme inputs for `v_dig`. `[E]` is **derived** (`ppm × 10⁻⁶ × N_prot / (N_A · V_DV)`), not an independent constant; code converts `kcat` to min⁻¹ as `60 × kcat[s⁻¹]`. Full citations: [`docs/enzyme_kinetics.md`](../enzyme_kinetics.md).

| Constant | Value | Units | Description |
|----------|------:|-------|-------------|
| `ppm_plm_1` | 752 | — | PaxDB Tao 2014 Dd2 abundance of plasmepsin-1 (input; `[E]` derived) |
| `kcat_plm_1` | 2.3 | s⁻¹ | Luker/Banerjee α33–34 peptide (native PM I) |
| `Km_plm_1` | 0.49×10⁻⁶ | M | Luker/Banerjee α33–34 peptide (native PM I) |
| `ppm_plm_2` | 1204 | — | PaxDB Tao 2014 Dd2 abundance of plasmepsin-2 (input; `[E]` derived) |
| `kcat_plm_2` | 11 | s⁻¹ | Luker/Banerjee α33–34 peptide (native PM II) |
| `Km_plm_2` | 2.6×10⁻⁶ | M | Luker/Banerjee α33–34 peptide (native PM II) |
| `ppm_hap` | 1373 | — | PaxDB Tao 2014 Dd2 abundance of HAP (input; `[E]` derived) |
| `kcat_hap` | 0.05 | s⁻¹ | Banerjee 2002 Table 1 (native HAP, α33–34) |
| `Km_hap` | 0.29×10⁻⁶ | M | Banerjee 2002 Table 1 (native HAP, α33–34) |
| `ppm_plm_4` | 3139 | — | PaxDB Tao 2014 Dd2 abundance of plasmepsin-4 (input; `[E]` derived) |
| `kcat_plm_4` | 1.05 | s⁻¹ | Banerjee 2002 Table 1 (recombinant PM IV, α33–34) |
| `Km_plm_4` | 0.33×10⁻⁶ | M | Banerjee 2002 Table 1 (recombinant PM IV, α33–34) |
| `ppm_fp_2` | 20.2 | — | PaxDB Tao 2014 Dd2 abundance of falcipain-2a (input; `[E]` derived) |
| `kcat_fp_2` | 0.79 | s⁻¹ | Ramjee 2006 best FP-2 FRET peptide |
| `Km_fp_2` | 0.9×10⁻⁶ | M | Ramjee 2006 best FP-2 FRET peptide |
| `ppm_fp_3` | 23.5 | — | PaxDB Tao 2014 Dd2 abundance of falcipain-3 (input; `[E]` derived) |
| `kcat_fp_3` | 0.204 | s⁻¹ | Ramjee 2006 FP-3 Leu-Arg FRET peptide |
| `Km_fp_3` | 4.0×10⁻⁶ | M | Ramjee 2006 FP-3 Leu-Arg FRET peptide |

## Assumptions

- One exponential law drives both cargo delivery and protease capacity (convenient coupling, not a claim of shared molecular control).
- Lipid effect remains the Model 2 `φ` rate penalty — still chemically inconsistent with lipid-mediated Hz.
- Fixed DV volume.

---

## Behaviour notes

- With host depletion, total Fe stays ≈ 106 fg; late host Fe remains if ∫ `f_exp` does not exhaust the pool by the end of the window.
- Relative to Models 1–2, Hz typically rises more realistically because uptake and digestion ramp together — but free-haem chemistry is still the Model 2 rate adjustment.

---

## Example

```python
from haem_kinetics.models.legacy.model3 import Model3

Model3().run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model3.png',
)
```
