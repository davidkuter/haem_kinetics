# Enzyme kinetic parameters (kcat, Km)

Shared assignment for `k_enzymes` in [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py). All values are **peptide-substrate** Michaelis–Menten parameters used as the best available proxy for haem release from Hb in the DV. They are **not** native-Hb turnover numbers.

## Modeling use

`v_dig` sums MM terms over PMs + falcipains (Models 1–6; Degradation is PMs-only):

```text
v_dig = 4 · Σ_i  (60 · kcat_i) · [E]_i,eff · [Hb]_tet / (Km_i + [Hb]_tet)
```

with `kcat` stored in s⁻¹ and converted to min⁻¹ in code. Peptide `kcat`/`Km` applied to tetrameric Hb haem-equivalents is an explicit approximation; if digestion is too fast relative to Combrink Dd2, the accountable next step is a dedicated Hb-efficiency factor grounded in protein-vs-peptide data — not an uncited rescale of `kcat`.

## Assigned values (code defaults)

| Key | Enzyme | kcat (s⁻¹) | Km (M) | Substrate / conditions | Primary source |
|-----|--------|-----------:|-------:|------------------------|----------------|
| `plm_1` | Plasmepsin I | **2.3** | **0.49×10⁻⁶** | Fluorogenic Hb α33–34 peptide (“33–34”); native enzyme (Luker) | Luker et al. 1996; Banerjee et al. 2002 Table 1 |
| `plm_2` | Plasmepsin II | **11** | **2.6×10⁻⁶** | Same peptide family; native enzyme (Luker) | Luker et al. 1996; Banerjee et al. 2002 Table 1 |
| `hap` | Histoaspartic protease (PM III) | **0.05** | **0.29×10⁻⁶** | Same 33–34 peptide; **native** FV-purified HAP; pH ≈ 5.7, 37 °C | Banerjee et al. 2002 Table 1 |
| `plm_4` | Plasmepsin IV | **1.05** | **0.33×10⁻⁶** | Same 33–34 peptide; recombinant PM IV; pH 5.4, 37 °C | Banerjee et al. 2002 Table 1 |
| `fp_2` | Falcipain-2 | **0.79** | **0.9×10⁻⁶** | Abz-Leu-Arg▼Phe-Pro-Tyr(NO₂)-Glu-NH₂ (best FP-2 FRET in study) | Ramjee et al. 2006 Table 1 |
| `fp_3` | Falcipain-3 | **0.204** | **4.0×10⁻⁶** | Abz-Leu-Leu-Arg▼Ala-Tyr(NO₂)-Glu-NH₂ (Leu-Arg motif; same study) | Ramjee et al. 2006 |

### Why Banerjee / Luker for the four plasmepsins

Banerjee et al. (*PNAS* 2002) report a single Table 1 for **PM I, PM II, HAP, and PM IV** on the **same** fluorogenic α33–34 substrate. PM I/II entries there are taken from Luker et al. (*Mol. Biochem. Parasitol.* 1996); HAP and PM IV were measured in that work. Using that table keeps all four FV aspartic proteases on a common substrate basis.

**Caveats from Banerjee (same paper):**

- HAP is ~20× less efficient than PM I/II on the α33–34 **peptide**, but faster than PM II/IV on **acid-denatured globin**; HAP does not efficiently cleave **native** Hb alone and synergizes with PM II.
- PM IV prefers globin over native Hb and is slower than PM II on protein substrates.
- Peptide MM parameters therefore overstate (or mis-order) relative contributions when the model substrate is treated as native Hb.

### Alternative recombinant PM I / PM II dataset (not used as default)

Siripurkpong et al. (*J. Biol. Chem.* 2002) Table I, recombinant enzymes, DABCYL-GABA-Glu-Arg-Met-Phe*Leu-Ser-Phe-Pro-GABA-EDANS, pH 5.0, 37 °C:

| Enzyme | Km (μM) | kcat (s⁻¹) |
|--------|--------:|----------:|
| PM I-WT | 0.917 ± 0.144 | 0.473 ± 0.172 |
| PM II-WT | 0.489 ± 0.080 | 2.313 ± 0.092 |

These differ substantially from Luker / Banerjee native values (especially PM II). Prefer Banerjee’s consistent four-enzyme table unless the model is deliberately switched to recombinant-only parameters.

Lolupiman et al. (*PLoS ONE* 2014) Table 1 compares further PM-I preparations (e.g. their own: Km ≈ 0.092 μM, kcat ≈ 0.344 s⁻¹ on a related FRET peptide at pH 4.5).

### HAP — values that must **not** be used

| Claimed source | Problem |
|----------------|---------|
| Bjelic & Åqvist, *Biochemistry* 2004 (`bi048252q`) | **Computational** mechanism study, not an experimental kcat/Km assignment |
| Xiao et al., *Protein Expr. Purif.* 2006 | Recombinant HAP: Km = 3.4 μM, kcat = 1.6×10⁻³ s⁻¹ on a related IQF peptide — valid data, but **not** native Banerjee HAP and far slower; do not mix with Banerjee PM I/II/IV without documenting the switch |
| Prior code comment citing Luker for HAP `kcat = 0.1` | Mis-citation; Banerjee Table 1 gives **0.05 s⁻¹** (Km 0.29 μM matches the old `2.98e-7` within rounding) |

### Falcipains

Ramjee et al. (*Biochem. J.* 2006) map FRET substrates for FP-2 and FP-3. Defaults use the highest-efficiency Leu-Arg motif peptide reported for FP-2 and a Leu-Arg containing peptide for FP-3 from the same tables (not Z-Leu-Arg-AMC-only Km without kcat).

Sijwali et al. (*Biochem. J.* 2001) report FP-3 Km on AMC peptides (e.g. Z-Leu-Arg-AMC ≈ 72 μM) — useful cross-check, but not the default pair.

**Limitation:** FRET/AMC peptide kinetics ≠ native Hb haemoglobinase rates (see also Subramanian et al., *PLoS ONE* 2009, Hb cleavage-site mapping).

## References

1. Luker KE, Francis SE, Gluzman IY, Goldberg DE. Kinetic analysis of plasmepsins I and II… *Mol. Biochem. Parasitol.* (1996) 79:71–78. [doi:10.1016/0166-6851(96)02651-5](https://doi.org/10.1016/0166-6851(96)02651-5)
2. Banerjee R, Liu J, Beatty W, Pelosof L, Klemba M, Goldberg DE. Four plasmepsins are active in the *P. falciparum* food vacuole… *PNAS* (2002) 99:990–995. [doi:10.1073/pnas.022630099](https://doi.org/10.1073/pnas.022630099) — **Table 1** (PM I/II from Luker; HAP & PM IV measured)
3. Siripurkpong P, Yuvaniyama J, Wilairat P, Goldberg DE. Active site contribution to specificity of plasmepsins I and II. *J. Biol. Chem.* (2002) 277:41009–41013. [doi:10.1074/jbc.M204852200](https://doi.org/10.1074/jbc.M204852200)
4. Ramjee MK, Flinn NS, Pemberton TP, Quibell M, Wang Y, Watts JP. Substrate mapping… falcipain-2, falcipain-3… *Biochem. J.* (2006) 399:47–57. [doi:10.1042/BJ20060422](https://doi.org/10.1042/BJ20060422)
5. Xiao H, Sinkovits AF, Bryksa BC, Ogawa M, Yada RY. Recombinant… histo-aspartic protease… *Protein Expr. Purif.* (2006) 49:88–94. [doi:10.1016/j.pep.2006.02.022](https://doi.org/10.1016/j.pep.2006.02.022) — alternative recombinant HAP only
6. Sijwali PS, Shenai BR, Gut J, Singh A, Rosenthal PJ. Expression and characterization of… falcipain-3. *Biochem. J.* (2001) 360:481–489. [doi:10.1042/bj3600481](https://doi.org/10.1042/bj3600481)
