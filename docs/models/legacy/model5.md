# Legacy Model 5

**Code:** [`haem_kinetics/models/legacy/model5.py`](../../../haem_kinetics/models/legacy/model5.py)  
**Up:** [Legacy index](README.md) · [Active models](../../models.md) · **Prev:** [Legacy Model 4](model4.md) · **Next:** [Legacy Model 6](model6.md)

**One change vs Model 4:** protease abundance follows a **logistic** maturation curve vs parasite age, independent of uptake `f_exp`.

Enzymes (PMs + FP2/3), peptide-table `kcat`, and Fe(III) path are unchanged from Model 4.

---

## What changed vs Model 4 (and why)

**Problem in Models 3–4:** `[E]_eff = f_exp(t) × [E]_pax` forces digestion capacity to track the uptake prefactor. Protease maturation is better described as a delayed rise toward a plateau during the trophozoite window — not the same clock as cargo delivery.

**Concrete change:**

| Item | Model 4 | Model 5 |
|------|---------|---------|
| Proteases | PMs + FP2/3 | **Same** |
| `[E](t)` | `f_exp(t) × [E]_pax` | `s_log(t) × [E]_pax` |
| Uptake | `f_exp` | **Same** (now **decoupled** from [E]) |
| Fe(III) / Hz | aq ⇄ lip → Hz | **Same** |

```text
s_log(t) = 1 / (1 + exp(−σ · (τ_age − t_mid)))
τ_age = 16 + t/60                 σ = 0.35 h⁻¹,  t_mid = 26 h
```

---

## Process schematic

```mermaid
flowchart LR
  Host["conc_hb_rbc"] -->|"f_exp(t)"| HbDV["conc_hb_dv"]
  HbDV -->|"PMs+FP x s_log"| Fe2["conc_fe2pp"]
  Fe2 -->|"k_ox x O2"| Fe3aq["conc_fe3pp_aq"]
  Fe3aq <--> Fe3lip["conc_fe3pp_lip"]
  Fe3lip -->|"k_hz"| Hz["conc_hz"]
```

---

## State variables

| Symbol | Meaning |
|--------|---------|
| `conc_hb_dv` | DV Hb (haem-eq) |
| `conc_fe2pp` | Fe(II)PPIX |
| `conc_fe3pp_aq` | Aqueous Fe(III)PPIX |
| `conc_fe3pp_lip` | Lipid-associated non-Hz Fe(III) |
| `conc_hz` | Haemozoin |
| `conc_hb_rbc` | Remaining host Hb |

Init: `[Hb_DV, Fe2, Fe3_aq, Fe3_lip, Hz]`.

---

## Governing equations

Uptake / digestion / oxidation (logistic enzyme clock):

```text
# Fractional exponential growth (uptake only)
f_exp(t) = a · b · exp(b · t)     a = 0.1578,  b = 0.001102

# Logistic enzyme maturation vs parasite age
s_log(t) = 1 / (1 + exp(−σ · (τ_age − t_mid)))
τ_age = 16 + t/60                 σ = 0.35 h⁻¹,  t_mid = 26 h

# Host → DV Hb uptake (M/min on V_DV)
v_up = f_exp(t) · [Hb]_RBC · V_RBC / V_DV

# Effective protease concentration (PMs + falcipains; logistic)
[E]_i,eff = s_log(t) · [E]_i
i ∈ {plm_1, plm_2, hap, plm_4, fp_2, fp_3}

# Haem release from Hb (MM sum; 4 haem-eq per tetramer)
v_dig = 4 · Σ_i  (60 · kcat_i) · [E]_i,eff · [Hb]_tet
                     / (Km_i + [Hb]_tet)

# Fe(II) → Fe(III) oxidation
v_ox  = k_fe2_ox · [Fe(II)] · [O2]

# Fe(III) → Fe(II) reduction (off: [O2−] = 0; aqueous pool)
v_red = k_fe3_red · [Fe(III)]_aq · [O2−]
```

Lipid partition target and exchange:

```text
# Aqueous fraction at lipid/aqueous partition equilibrium
φ = (1 − f_lip) / (1 + f_lip + f_lip · K_partition)

# Effective lip/aq equilibrium ratio
K_eff = (1 − φ) / φ = [Fe(III)]_lip / [Fe(III)]_aq |eq

# Aqueous ⇄ lipid Fe(III) exchange
v_ex = k_lipid_ex · ( [Fe(III)]_aq − [Fe(III)]_lip / K_eff )

# Haemozoin formation (from lipid pool; full k_hz)
v_hz = k_hz · [Fe(III)]_lip
```

ODEs:

```text
# DV haemoglobin
d[Hb]_DV / dt       = v_up − v_dig

# Free Fe(II)PPIX
d[Fe(II)] / dt      = v_dig + v_red − v_ox

# Aqueous Fe(III)PPIX
d[Fe(III)]_aq / dt  = v_ox − v_red − v_ex

# Lipid-associated Fe(III)PPIX
d[Fe(III)]_lip / dt = v_ex − v_hz

# Haemozoin
d[Hz] / dt          = v_hz

# Remaining host RBC Hb
d[Hb]_RBC / dt      = − v_up · V_DV / V_RBC
```

---

## Constants used

| Constant | Value | Units | Description |
|----------|------:|-------|-------------|
| `a` | 0.1578 | — | Prefactor in fractional exponential growth `f_exp(t)` (uptake only) |
| `b` | 0.001102 | min⁻¹ | Rate constant in fractional exponential growth `f_exp(t)` (uptake only) |
| `σ` | 0.35 | h⁻¹ | Steepness of the logistic enzyme-maturation curve |
| `t_mid` | 26 | h | Midpoint age (h post-invasion) of the logistic enzyme-maturation curve |
| `t0` | 16 | h | Simulation time origin in hours post-invasion |
| `V_RBC` | 90×10⁻¹⁵ | L | Volume of the host red blood cell |
| `V_DV` | 1×10⁻¹⁵ | L | Fixed digestive-vacuole volume used for M ↔ fg conversion |
| `N_A` | 6.022×10²³ | mol⁻¹ | Avogadro's number |
| `N_prot` | 1.9×10⁸ | — | Average number of proteins per *P. falciparum* cell |
| `k_fe2_ox` | 193800 | min⁻¹ | Rate constant for Fe(II)PPIX oxidation by O₂ |
| `[O2]` | 1×10⁻³ | M | Dissolved oxygen concentration in the DV |
| `k_fe3_red` | 180×10⁻⁹ | — | Rate constant for Fe(III)PPIX reduction by O₂⁻ (inactive when `[O2−]` = 0) |
| `[O2−]` | 0 | M | Superoxide concentration (taken as zero due to SOD) |
| `f_lip` | 0.016 | — | Fractional volume of lipid nanospheres relative to the DV |
| `K_partition` | 398 | — | Equilibrium partition coefficient of Fe(III)PPIX into lipid |
| `K_eff` | ≈ 6.50 | — | Equilibrium ratio `[Fe(III)]_lip / [Fe(III)]_aq` |
| `k_lipid_ex` | 50 | min⁻¹ | Rate constant for aqueous ⇄ lipid Fe(III) exchange |
| `k_hz` | 0.12 | min⁻¹ | First-order rate constant for haemozoin formation from lipid Fe(III) |

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

## Behaviour vs Model 4

- Early digestion weaker, later stronger → HbDV/Hz *shape* can change even if end Hz is similar.
- Free-haem chemistry remains aq ⇄ lip → Hz at full `k_hz` on the lipid pool.

---

## Example

```python
from haem_kinetics.models.legacy.model5 import Model5

Model5().run(
    t=[0, 1700],
    init=[0.018, 0.0, 0.0, 0.0, 0.36],
    t_eval=range(0, 1700, 20),
    plot='examples/model5.png',
)
```
