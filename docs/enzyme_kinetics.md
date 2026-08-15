# Enzyme kinetic parameters (kcat, Km)

Shared assignment for `k_enzymes` in [`haem_kinetics/components/constants.py`](../haem_kinetics/components/constants.py). All values are **peptide-substrate** Michaelis–Menten parameters used as the best available proxy for haem release from Hb in the DV. They are **not** native-Hb turnover numbers.

## Modeling use

`v_dig` sums MM terms over PMs + falcipains on Models 1–3 (Degradation is PMs-only). Models 4a/4b restrict the **native tetramer** to PM I, PM II, and FP-2 and keep this peptide table for **globin / fragments** only (see [Native tetramer vs peptide](#native-tetramer-vs-peptide-models-4a4b)).

```text
v_dig = 4 · Σ_i  (60 · kcat_i) · [E]_i,eff · [Hb]_tet / (Km_i + [Hb]_tet)
```

with `kcat` stored in s⁻¹ and converted to min⁻¹ in code. Peptide `kcat`/`Km` applied to tetrameric Hb haem-equivalents is an explicit approximation used on Models 1–3. Model 4 encodes the substrate distinction instead of an uncited rescale of `kcat`.

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

## Native tetramer vs peptide (Models 4a/4b)

Assay DV Hb is native tetramer. The table above is **peptide** (Banerjee/Luker α33–34; Ramjee FRET). Goldberg’s ordered pathway: PM I nicks native Hb at α33–34; PM II prefers denatured globin; HAP and PM IV prefer globin over native Hb; FP-2 *does* hydrolyze native Hb at vacuolar pH (Shenai 2000, updating Gluzman’s cysteine-protease gel for this gene product).

**Native-competent set:** PM I, PM II, FP-2. Code: `NATIVE_HB_ENZYMES` in [`native_hb.py`](../haem_kinetics/models/native_hb.py); parameters in `k_enzymes_native`.

**Globin / fragments:** all six proteases with `k_enzymes` (this peptide table).

No cited in vitro native-Hb digest (Gluzman 1994 Fig. 3 methods; Banerjee Hb gels) states enzyme amount and [Hb] in a form that converts to `kcat` in s⁻¹. The native table is therefore **not** fit to Garnie standing Hb or Dd2 digestion rates:

| Key | Native kcat (s⁻¹) | Source of the number |
|-----|------------------:|----------------------|
| `plm_1` | 2.3 | Peptide `kcat`, labelled provisional (initiator) |
| `plm_2` | 2.3 **not** 11 | Gluzman: equal globin-degrading units, PM I ≫ PM II on native; peptide 11 s⁻¹ would invert that |
| `fp_2` | 0.79 | Peptide `kcat` kept (Shenai: FP-2 cuts native Hb; peptide `Vmax` already ~4 fg h⁻¹) |

HAP, PM IV, and FP-3 have no native-table entry. Full encoding: [model4.md](models/model4.md).

## Attempted native-Hb `kcat_app` reconstruction

Model 4 left PM I on peptide `kcat` 2.3 s⁻¹ because no paper reports native-tetramer turnover in s⁻¹. An attempted ladder step was: convert a **cited native-Hb assay** to `kcat_app` and put that number on the 4a pathway — not a Garnie-fitted `f_Hb`. The conversion requires enzyme amount (mol, nM, or μg + MW), [Hb], time, and a cleaved fraction. That reconstruction is **not** the current Model 5 (HTV cargo on 4b: [model5.md](models/model5.md)).

**Diagnostic (not a parameter).** Saturated PM I peptide `Vmax` at PaxDB amount is ~440 fg Fe h⁻¹; at `s_PM(20 h) ≈ 0.41` that is still ~180 fg h⁻¹. Uptake / Garnie Dd2 digestion are ~1–5 fg h⁻¹. Matching `v_nick ≈ v_up` for ~1 fg standing Hb would need native PM I `kcat` roughly 10²–10³× below 2.3 s⁻¹. That bound is recorded here so a future measured `kcat` can be compared to it; it was **not** used to choose a number.

### What each paper supplies

**Gluzman et al. *J. Clin. Invest.* (1994)** — PDF methods / Fig. 3.

| Input | Value | Usable for `kcat`? |
|-------|-------|-------------------|
| Enzyme in Fig. 3 | **75 mU** (equal globin-units of AH I and AH II) | Activity units, not moles |
| Substrate | unlabeled human Hb **1.6 mg/ml**, pH 5.0, **2 h** | [Hb] and time yes |
| Cleaved fraction | silver-stained SDS-PAGE; “considerably more active” for AH I vs II | no densitometry / no *f* |
| Unit definition | **1 U = 1 μg/h [¹⁴C]globin → TCA-soluble fragments**, pH 5.0, 37 °C | denatured globin, not native nick |
| Protein mass / U/mg | homogeneous by SDS-PAGE; **no U/mg, no μg in 75 mU** | **[E] missing** |

If a purification table had given specific activity *S* (U/mg), enzyme mass in the Fig. 3 lane would be `0.075 mg / S` and moles would follow from PM I MW (~37 kDa). That table is not in Gluzman. Treating 75 mU as if it were native-Hb nick rate would also be wrong: the unit is TCA-soluble **globin** fragments.

LC-MS incubations (20 mU + 3.2 μg Hb, 15 h) are overnight endpoint maps, not initial rates.

**Francis et al. *EMBO J.* (1994)** — full PDF.

Cloning, localization, and SC-50083 inhibition. Aspartic hemoglobinase I “was purified as previously described (Goldberg et al., 1991).” Characterization of the purified proteases “will be described elsewhere” (the Gluzman paper). No μg enzyme × μg Hb time course.

**Goldberg et al. *J. Exp. Med.* (1991)** — methods OCR (DocsLib / JEM PDF text).

Purification to a single major 40 kDa band; **282-fold** over starting material (Table 1). “Human hemoglobin was incubated with the purified enzyme for 30 min” for N-terminal sequencing of fragments. Hemoglobinase assay “details … previously (16)” = Goldberg *PNAS* 1990 (vacuole assay), not a native-Hb `kcat`. Table 1’s protein (mg) and specific activity (U/mg) are **not stated in the running text**; inventing them from a silver-stained gel is not a conversion.

**Shenai et al. *J. Biol. Chem.* (2000).**

Shows FP-2 hydrolyzes native Hb at vacuolar pH 5.5 / 1 mM GSH (updates Gluzman’s cysteine-protease gel for this gene product). Peptide `kcat`/`Km` are the kinetic table. Native-Hb evidence is gel hydrolysis, not `[E]`, `[Hb]`, time, and *f* in a form that converts to s⁻¹. Later papers (e.g. Hanspal 2002: 100 nM rFP-2 + 3 μg Hb in 25 μl, 60 min; Marques et al.: 20 nM + 100 μg/ml Hb, 90 min) are **not** Shenai’s methods and were not mixed in.

**Vander Jagt et al. *BBA* (1992)** 1122:256.

Three vacuole aspartic activities (M1, M2, M3). Native and denatured Hb are both substrates. The fractions are not identified as PM I, and no nmol min⁻¹ mg⁻¹ table was available that could be assigned to PM I.

### Outcome

**`[E]` in moles is still missing.** No `kcat_app` was invented. That attempted step is closed: Model 5 is HTV inaccessible cargo on the 4b pathway, not a reconstructed native turnover. The remaining gap for a measured native-Hb `kcat` is recorded here so a future number can be compared to the diagnostic bound above. Lipid / basal Hm is [Model 6](models/model6.md) (aq ⇄ lip at literature `k_hz`), not a substitute for that number.

## References

1. Luker KE, Francis SE, Gluzman IY, Goldberg DE. Kinetic analysis of plasmepsins I and II… *Mol. Biochem. Parasitol.* (1996) 79:71–78. [doi:10.1016/0166-6851(96)02651-5](https://doi.org/10.1016/0166-6851(96)02651-5)
2. Banerjee R, Liu J, Beatty W, Pelosof L, Klemba M, Goldberg DE. Four plasmepsins are active in the *P. falciparum* food vacuole… *PNAS* (2002) 99:990–995. [doi:10.1073/pnas.022630099](https://doi.org/10.1073/pnas.022630099) — **Table 1** (PM I/II from Luker; HAP & PM IV measured)
3. Siripurkpong P, Yuvaniyama J, Wilairat P, Goldberg DE. Active site contribution to specificity of plasmepsins I and II. *J. Biol. Chem.* (2002) 277:41009–41013. [doi:10.1074/jbc.M204852200](https://doi.org/10.1074/jbc.M204852200)
4. Ramjee MK, Flinn NS, Pemberton TP, Quibell M, Wang Y, Watts JP. Substrate mapping… falcipain-2, falcipain-3… *Biochem. J.* (2006) 399:47–57. [doi:10.1042/BJ20060422](https://doi.org/10.1042/BJ20060422)
5. Xiao H, Sinkovits AF, Bryksa BC, Ogawa M, Yada RY. Recombinant… histo-aspartic protease… *Protein Expr. Purif.* (2006) 49:88–94. [doi:10.1016/j.pep.2006.02.022](https://doi.org/10.1016/j.pep.2006.02.022) — alternative recombinant HAP only
6. Sijwali PS, Shenai BR, Gut J, Singh A, Rosenthal PJ. Expression and characterization of… falcipain-3. *Biochem. J.* (2001) 360:481–489. [doi:10.1042/bj3600481](https://doi.org/10.1042/bj3600481)
7. Gluzman IY, Francis SE, Oksman A, Smith CE, Dole K, Goldberg DE. Order and specificity of the *Plasmodium falciparum* hemoglobin degradation pathway. *J. Clin. Invest.* (1994) 93:1602–1608. [doi:10.1172/jci117140](https://doi.org/10.1172/jci117140)
8. Goldberg DE, Slater AF, Beavis R, Chait B, Cerami A, Henderson GB. Hemoglobin degradation in the human malaria pathogen *Plasmodium falciparum*: a catabolic pathway initiated by a specific aspartic protease. *J. Exp. Med.* (1991) 173:961–969. [doi:10.1084/jem.173.4.961](https://doi.org/10.1084/jem.173.4.961)
9. Shenai BR, Sijwali PS, Singh A, Rosenthal PJ. Characterization of native and recombinant falcipain-2… *J. Biol. Chem.* (2000) 275:29000–29010. [doi:10.1074/jbc.M004459200](https://doi.org/10.1074/jbc.M004459200)
10. Liu J, Gluzman IY, Drew ME, Goldberg DE. The role of *Plasmodium falciparum* food vacuole plasmepsins. *J. Biol. Chem.* (2005) 280:25416–25424. [doi:10.1074/jbc.M412086200](https://doi.org/10.1074/jbc.M412086200)
11. Nasamu AS, Polino AJ, Istvan ES, Goldberg DE. Malaria parasite plasmepsins: more than just plain old degradative pepsins. *J. Biol. Chem.* (2020) 295:8425–8441. [doi:10.1074/jbc.REV120.009309](https://doi.org/10.1074/jbc.REV120.009309)
12. Francis SE, Gluzman IY, Oksman A, Knickerbocker A, Mueller R, Bryant ML, Sherman DR, Russell DG, Goldberg DE. Molecular characterization and inhibition of a *Plasmodium falciparum* aspartic hemoglobinase. *EMBO J.* (1994) 13:306–317. [doi:10.1002/j.1460-2075.1994.tb06263.x](https://doi.org/10.1002/j.1460-2075.1994.tb06263.x)
13. Vander Jagt DL, Hunsaker LA, Campos NM, Scaletti JV. Localization and characterization of hemoglobin-degrading aspartic proteinases from the malarial parasite *Plasmodium falciparum*. *Biochim. Biophys. Acta* (1992) 1122:256–264. [doi:10.1016/0167-4838(92)90401-X](https://doi.org/10.1016/0167-4838(92)90401-X)
