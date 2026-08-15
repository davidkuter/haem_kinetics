# Garnie heme fractionation vs pHrodo

**Scoring target:** Dd2 Hb / Hm / Hz time courses in [`experimental_data.py`](../haem_kinetics/components/experimental_data.py).  
**Paper:** Garnie, Egan & Wicht, *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z)  
**Assay:** Combrinck *et al.*, *Malar. J.* (2015) 14:253. [doi:10.1186/s12936-015-0729-9](https://doi.org/10.1186/s12936-015-0729-9)

Garnie reports **two different experiments**. They are not interchangeable. The ODE ladder is scored against **fractionation**, not pHrodo intensity.

---

## Fractionation (Fig. 4 — Hb, Hm, Hz)

Not isolated digestive vacuoles. Protocol (Garnie methods; Combrinck 2015):

1. Harvest infected RBCs; add saponin (~0.15%); centrifuge; **PBS-wash the pellet** until host cytosol haemoglobin is gone.
2. Freeze the trophozoite pellet; thaw; add MilliQ water; **sonicate**.
3. Sequential extraction:
   - aqueous / HEPES supernatant → **Hb** (protein-bound haem);
   - SDS + pyridine on the remaining pellet → **free Hm**;
   - NaOH on what is left → **Hz**.
4. Quantify each fraction as Fe(III)PPIX–pyridine; divide by flow-cytometry cell count → fg Fe/cell.

Chemically, the Hb fraction is **all aqueous-extractable protein-bound haem in the saponin pellet**: protease-accessible DV lumen, inner (PVM-derived) vesicles already inside the DV, any unfused cytostomal / HTV cargo still in the parasite, and any other soluble Hb that survived the washes. The tubes do not spatially resolve those pools.

**How Garnie uses the numbers:** they treat Hb, Hm, and Hz as **DV-localized**. The paper states that the heme-containing species localize to the DV, divides Hb and Hm by pHrodo lumen volume to get relative lumen concentrations (Fig. 5A), and excludes Hz from that volume because crystals are not aqueous space. That is an interpretation from known DV biology, not from having isolated DVs.

Host RBC haemoglobin is **not** in the table: saponin plus PBS washes remove it. Remaining Hb is parasite-associated. Total Fe in the three fractions is the **DV-associated inventory** as Garnie uses it (Hz-dominated on Dd2).

---

## pHrodo confocal (Fig. 2 — lumen volume and “uptake”)

A separate live-cell experiment: RBC ghosts are loaded with pHrodo-dextran; parasites re-invade; fluorescence is recorded only at low pH. That reports **acidic aqueous DV lumen**.

It is **not** `dFe/dt`:

- Neutral cytostomal vesicles (Klonis 2007: ER-like pH) do not fluoresce.
- Hz crystals are excluded from the measured lumen (Garnie: pHrodo volume is aqueous space, not pigment).
- Fig. 2C is standing probe intensity in that lumen, not a delivery rate.

For Dd2, pHrodo intensity plateaus then falls while fractionation total Fe keeps rising. That is expected if delivery is converted to Hz, which pHrodo cannot see. Do **not** replace empirical `f_exp` with Fig. 2C. A later numbered model could test a cited fluid-phase delivery law; that is not a silent swap.

`V_DV(t)` in the code is the Garnie **lumen** schedule (Gompertz then collapse), used as shared molar bookkeeping — not as the uptake law.

---

## What this means for `f_exp` and Model 5

**`f_exp` is already host → DV-associated inventory.** Model 2a fitted `a, b` to cumulative fractionation Fe (`Hb + Hm + Hz`) as first-order remaining host (ODE **R² = 0.612**). Model 2b is the same law with a faster late phase after Fig. 5B’s 29 h break (still × leftover host). Neither is a cytostome assay, pHrodo, or “delivery to the parasite” with a later DV step. See [model2.md](models/model2.md).

**Inner-vesicle cargo (`conc_hb_htv`) is not that later DV step.** Model 5 exists because lumen `Vmax` ≫ uptake, so protease-accessible standing Hb collapses (~0 vs assay ~1–2 fg). The pool is **inner vesicles already inside the DV** after outer-membrane fusion (Yayon 1984; Klemba *JCB* 2004): still protein-bound Fe, not mixed with soluble proteases until the inner membrane lyses. Sonication puts that Hb in the assay Hb fraction, so plots score `conc_hb_htv + conc_hb_dv`. The cargo is **not** in pHrodo `V_DV(t)` (same reason Hz is excluded: not aqueous lumen). See [model5.md](models/model5.md).

Do not drop the pool (returns to Model 4: assay Hb ~0). Do not stop counting it as assay Hb (Garnie’s tubes would still see it).

---

## References

- Garnie LF, Egan TJ, Wicht KJ. *Commun. Biol.* (2025) 8:1564. [doi:10.1038/s42003-025-08991-z](https://doi.org/10.1038/s42003-025-08991-z)
- Combrinck JM, Fong KY, Gibhard L, Smith PJ, Wright DW, Egan TJ. Optimization of a multi-well colorimetric assay to determine haem species in *Plasmodium falciparum*. *Malar. J.* (2015) 14:253. [doi:10.1186/s12936-015-0729-9](https://doi.org/10.1186/s12936-015-0729-9)
- Klonis N, et al. *Biochem. J.* (2007) 407:343–354. [doi:10.1042/BJ20070934](https://doi.org/10.1042/BJ20070934)
- Klemba M, Beatty W, Gluzman I, Goldberg DE. *J. Cell Biol.* (2004) 164:47–56. [doi:10.1083/jcb.200307147](https://doi.org/10.1083/jcb.200307147)
