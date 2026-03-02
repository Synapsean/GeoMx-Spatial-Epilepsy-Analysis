# Targeted Intersection Analyses — Log
**GeoMx Spatial Transcriptomics | Chronic KA vs PBS | Mouse Hippocampus**

Script: `scripts/targeted_analyses.R`
Results: `results/targeted/`

Cross-reference: [SELECTED_COMPARISONS.md](SELECTED_COMPARISONS.md) | [INTERPRETATION.md](INTERPRETATION.md)

---

## Quick Reference Table

| ID | Question | Comparison 1 | Comparison 2 | Genes Found | Output |
|----|----------|-------------|-------------|------------|--------|
| [A01](#a01) | KA-specific genes in CA3 Ipsi Neurons | `Side_KA_CA3_Neuron` | `Trt_CA3_Ipsi_Neuron` | **1** | A01_CA3_Ipsi_Neuron_KA_specific.csv |
| [A02](#a02) | KA-specific genes in CA3 Ipsi Astrocytes | `Side_KA_CA3_Astro` | `Trt_CA3_Ipsi_Astro` | **66** | A02_CA3_Ipsi_Astro_KA_specific.csv |
| [A03a](#a03) | KA-amplified CA3/CA1 regional signal — Neurons | `Reg_KA_Neuron` vs `Reg_PBS_Neuron` | — | **148** | A03a_KA_amplified_regional_Neuron.csv |
| [A03b](#a03) | KA-amplified CA3/CA1 regional signal — Microglia | `Reg_KA_Micro` vs `Reg_PBS_Micro` | — | **39** | A03b_KA_amplified_regional_Micro.csv |
| [A03c](#a03) | KA-amplified CA3/CA1 regional signal — Astrocytes | `Reg_KA_Astro` vs `Reg_PBS_Astro` | — | **57** | A03c_KA_amplified_regional_Astro.csv |
| [A04](#a04) | Contralateral compensatory genes in CA3 Neurons | `Side_KA_CA3_Neuron` (Contra>Ipsi) | `Side_PBS_CA3_Neuron` (not sig) | **20** | A04_Contralateral_protective_Neuron_CA3.csv |
| [A05](#a05) | Interaction-term candidates cross-validated by regional contrast | `Int_Ipsi_Neuron` | `Reg_KA_Neuron` vs `Reg_PBS_Neuron` | **54** | A05_Interaction_Neuron_cross_validated.csv |
| [A06](#a06) | Pan-cell-type KA consensus at CA3 Ipsi | `Trt_CA3_Ipsi_Neuron/Micro/Astro` | — | **0** | A06_Pan_celltype_KA_consensus_CA3_Ipsi.csv |
| [A07](#a07) | PBS-baseline-only regional genes (lost CA3 identity) | `Reg_PBS_Neuron` (FDR-sig) | `Reg_KA_Neuron` (not FDR-sig) | **13** | A07_Lost_CA3_identity_PBS_only_Neuron.csv |

---

## Thresholds Used

| Parameter | Value | Notes |
|-----------|-------|-------|
| Nominal p-value (`P_NOM`) | 0.05 | Used for Trt_ comparisons (0 FDR-sig genes) |
| Wide p-value (`P_WIDE`) | 0.10 | Available but not used in final output |
| Min \|logFC\| | 0.5 | For intersection analyses |
| FDR for Reg_ comparisons | 0.05 | Global `FDR_THRESHOLD` from 01_setup.R |
| LFC for Reg_ comparisons | 1.0 | Global `LFC_THRESHOLD` from 01_setup.R |
| Concordance | Required | Same direction in both comparisons |

> **Why nominal p<0.05 for Treatment comparisons?**
> All `Trt_` (KA vs PBS) comparisons return 0 FDR-significant genes. This is expected:
> treatment effect is captured by laterality (Ipsi>Contra) not by direct KA/PBS comparison 
> once the random mouse effect and regional context are modeled. Nominal p<0.05 is used as 
> a hypothesis-generating threshold.

---

## A01 — CA3 Ipsi Neuron KA-specific {#a01}

**Question:** Which genes are specifically up/downregulated in KA CA3 ipsilateral neurons?  
**Strategy:** Dual-confirmation — must be significant in BOTH:
- `Side_KA_CA3_Neuron` (Ipsi vs Contra in KA CA3 neurons): confirms laterality
- `Trt_CA3_Ipsi_Neuron` (KA vs PBS in CA3 Ipsi neurons): confirms treatment specificity

**Thresholds:** p < 0.05 in both, |logFC| ≥ 0.5 in both, concordant direction

**Result: 1 gene**

| Gene | logFC (Ipsi/Contra KA) | logFC (KA/PBS CA3-Ipsi) | Direction |
|------|----------------------|------------------------|-----------|
| Zar1l | -0.59 | -0.83 | Down in Ipsi-KA |

**Interpretation:** Neurons show almost no cell-type-specific KA response in CA3 Ipsi that isn't shared with the contralateral side. The dominant gene expression changes in neurons are regional (CA3 vs CA1), not laterality-driven. `Zar1l` (Zygote Arrest-1 Like) is a mRNA stability regulator with limited known neuronal role.

**Biological note:** Contrast with A02 (astrocytes, 66 genes) — glial cells show far more lateralized KA responses than neurons in this dataset.

---

## A02 — CA3 Ipsi Astrocyte KA-specific {#a02}

**Question:** Which genes are specifically dysregulated in KA CA3 ipsilateral astrocytes?  
**Strategy:** Same dual-confirmation as A01, for astrocytes.

**Thresholds:** p < 0.05 in both, |logFC| ≥ 0.5 in both, concordant direction

**Result: 66 genes** (1 DOWN-only gene in A01 vs 66 here — astrocytes dramatically more reactive)

**Top genes by combined p-value rank score:**

| Gene | logFC (Side KA) | logFC (Trt CA3-Ipsi) | Mean |FC| | Known Role |
|------|----------------|---------------------|----------|------------|
| Them6 | -1.65 | -0.94 | 1.29 | Acyl-thioester hydrolase, lipid metabolism |
| Rgs14 | -1.81 | -1.21 | 1.51 | Synaptic plasticity suppressor (CA2-specific normally) |
| Becn1 | -1.26 | -1.30 | 1.28 | Autophagy initiation (BECLIN-1) |
| Tmem35a | -1.05 | -1.18 | 1.12 | TM protein, neuronal activity |
| Asph | -1.38 | -0.99 | 1.19 | Aspartyl hydroxylase |
| Hapln1 | +1.21 | +0.81 | 1.01 | ECM link protein, neuroinflammation |
| Gzma | +1.12 | +1.17 | 1.15 | Granzyme A — astrocytic immune effector |
| Mro | +1.29 | +0.77 | 1.03 | Maestro, oxidative stress |

**Direction split:** Most (≈50/66) are DOWN in KA Ipsi astrocytes  
**Notable:** `Gzma` (Granzyme A) in astrocytes suggests an innate immune/cytotoxic component. `Hapln1` upregulation may reflect reactive gliosis and ECM remodeling. `Becn1` (BECLIN-1) decrease suggests impaired autophagy.

**Follow-up priority:** Validate Gzma, Hapln1, Becn1 expression by ISH or IHC in this anatomical context.

---

## A03 — KA-Amplified Regional Differences {#a03}

**Question:** Which CA3/CA1 gene expression differences are LARGER after KA seizures than at baseline?  
**Strategy:** Compare regional logFC in KA vs PBS — same direction, |logFC_KA| > |logFC_PBS| + 0.5  
**Filter:** FDR < 0.05 in KA (signals must be robust in the epilepsy condition)

### A03a — Neurons (148 genes)

The largest class. These 148 genes have an enhanced CA3 vs CA1 transcriptional signature after KA.

**Top 10 by delta logFC (seizure amplification):**

| Gene | logFC KA | logFC PBS | Delta | Direction | Function |
|------|----------|-----------|-------|-----------|----------|
| Iyd | +2.30 | +0.84 | +1.46 | CA3>CA1 | Iodotyrosine deiodinase |
| Ephb1 | +2.15 | +0.98 | +1.17 | CA3>CA1 | Ephrin receptor, axon guidance |
| Fibcd1 | -3.67 | -2.53 | -1.14 | CA1>CA3 | Fibrinogen C domain, ECM |
| Islr2 | +3.27 | +2.14 | +1.13 | CA3>CA1 | Ig superfamily, neuronal |
| Cpne4 | +4.33 | +3.32 | +1.01 | CA3>CA1 | Copine, calcium binding |
| Wfs1 | -3.46 | -2.51 | -0.94 | CA1>CA3 | Wolframin, ER stress |
| Doc2b | -3.52 | -2.58 | -0.93 | CA1>CA3 | Ca2+-dependent exocytosis |
| Calb1 | -2.17 | -1.25 | -0.92 | CA1>CA3 | Calbindin-D28K |
| Ncald | +3.45 | +2.47 | +0.98 | CA3>CA1 | Neurocalcin delta |
| Grik4 | +3.11 | +2.24 | +0.87 | CA3>CA1 | Kainate receptor subunit 4 |

**Key finding:** `Calb1` (calbindin) is normally CA1-enriched and serves as a neuroprotective marker. Its CA1>CA3 gradient is AMPLIFIED in KA, suggesting loss of calbindin in KA-vulnerable CA3 neurons. `Grik4` (kainate receptor) being further amplified in CA3 after KA is biologically consistent with KA-induced excitotoxicity.

### A03b — Microglia (39 genes)

| Gene | logFC KA | logFC PBS | Delta | Function |
|------|----------|-----------|-------|----------|
| Pgm2l1 | +1.62 | +0.75 | +0.87 | Phosphoglucomutase-like |
| Bhlhe22 | +0.88 | +0.05 | +0.83 | Transcription factor, microglial identity |
| Hapln4 | +1.88 | +1.17 | +0.70 | ECM link protein |
| Pcp4 | +1.88 | +1.22 | +0.66 | PEP-19, calmodulin inhibitor |
| Cpne4 | +2.15 | +1.52 | +0.62 | Copine (shared with neurons) |
| Ncald | +2.21 | +1.61 | +0.60 | Neurocalcin delta (shared) |

**Note:** Several genes (Cpne4, Ncald, Hapln4, Pcp4) are shared with the neuron A03 list, suggesting these CA3 identity markers are adopted by or transferred to resident microglia in the epileptic brain.

### A03c — Astrocytes (57 genes)

| Gene | logFC KA | logFC PBS | Delta | Function |
|------|----------|-----------|-------|----------|
| Usp1 | -1.13 | -0.13 | -1.01 | Ubiquitin protease — DNA damage response |
| Sst | +1.30 | +0.34 | +0.96 | Somatostatin — interneuron marker in astrocytes? |
| Gng3 | +1.73 | +0.82 | +0.91 | G-protein gamma 3 |
| Acsl4 | +1.41 | +0.63 | +0.78 | Long-chain acyl-CoA ligase, ferroptosis |
| Il34 | +1.20 | +0.43 | +0.77 | CSF1R ligand, microglial activation |
| Calb1 | -1.06 | -0.35 | -0.71 | Calbindin (shared with Neuron A03 list) |

**Notable:** `Acsl4` upregulation in astrocytes is a canonical ferroptosis marker. `Il34` is a CSF1R ligand — astrocyte-to-microglia signaling pathway amplified in KA CA3.

---

## A04 — Contralateral Compensatory Genes {#a04}

**Question:** Which neuronal genes are elevated on the LESS-AFFECTED contralateral CA3 in KA, but not in PBS? These may represent neuroprotective/homeostatic responses.  
**Strategy:** `Side_KA_CA3_Neuron` with logFC < -0.5 (Contra > Ipsi) at p < 0.05, genes absent from PBS side comparisons

**Result: 20 genes** (CA3 KA-specific contralateral elevation)
- Strict (both CA3 and CA1): 0 genes
- CA3 only: 20 genes

**Top candidates by p-value:**

| Gene | LogFC (Contra vs Ipsi in KA) | p-value | Known Role |
|------|------------------------------|---------|------------|
| Nnat | -0.69 | 0.0067 | Neuronatin — neonatal brain, neuroprotection |
| Pabpn1l | -0.76 | 0.0087 | RNA binding, poly-A regulation |
| Syce2 | -0.76 | 0.0103 | Synaptonemal complex |
| Jun | -0.55 | 0.0151 | AP-1 transcription factor, immediate-early gene |
| Fam107a | -0.52 | 0.0110 | DRR1, stress-protective |
| Fmn1 | -0.50 | 0.0204 | Formin-1, cytoskeletal remodeling |

**Interpretation:** `Jun` (AP-1/c-Jun) elevation on the contralateral side is notable — it's an immediate-early gene associated with neuroprotective responses in non-injured neurons. `Nnat` (neuronatin) protects against ER stress. `Fam107a/DRR1` is a stress-response protein. Together these suggest the contralateral CA3 activates a protective transcriptional program in the absence of direct seizure damage.

---

## A05 — Interaction Term Cross-Validation {#a05}

**Question:** Which genes from the interaction term (CA3 × Ipsi × KA) are independently confirmed by comparing the CA3/CA1 regional fold-change between KA and PBS animals?  
**Strategy:** Top 300 from `Int_Ipsi_Neuron` where the difference (`logFC_KA - logFC_PBS`) is ≥0.5 and in the same direction as the interaction term

**Result: 54 genes confirmed by both approaches**

**Top confirmed interaction genes:**

| Gene | Int logFC | Delta Regional | logFC_KA | logFC_PBS | Interpretation |
|------|-----------|---------------|----------|-----------|----------------|
| Iyd | +1.72 | +1.46 | +2.30 | +0.84 | Strong CA3 amplification in KA |
| Zfhx4 | -1.37 | -0.87 | -1.59 | -0.71 | CA1>CA3 identity marker, amplified |
| Stac2 | +1.13 | +0.93 | +1.37 | +0.44 | Neuronal calcium sensor |
| Lats2 | +0.89 | +0.64 | +1.72 | +1.07 | Hippo pathway kinase |
| Itgb5 | -1.35 | -0.91 | -0.23 | +0.69 | Integrin beta-5, ECM |
| Fnbp1l | -1.17 | -0.68 | -1.23 | -0.55 | F-BAR domain, actin |
| Arhgef26 | +1.01 | +0.97 | +1.69 | +0.72 | Rho GEF, possibly neurotoxic |
| Cnr1 | -1.02 | -0.66 | +0.26 | +0.92 | CB1 receptor — LOSS of CA1 enrichment in KA |

**Biological note:** `Cnr1` (CB1 receptor) is normally CA1-enriched, but in KA this CA1 enrichment is LOST (logFC_PBS = +0.92 but logFC_KA = +0.26). Loss of CB1 in CA1 may reduce endocannabinoid-mediated inhibition and contribute to continued network hyperexcitability. This is a known finding in epilepsy research and serves as internal validation.

---

## A06 — Pan-Cell-Type KA Nominal Consensus {#a06}

**Question:** Are there any genes where KA treatment has the same directional effect in ALL THREE cell types at CA3 Ipsi?  
**Strategy:** p < 0.05 and same direction in `Trt_CA3_Ipsi_Neuron`, `Trt_CA3_Ipsi_Micro`, AND `Trt_CA3_Ipsi_Astro`

**Result: 0 genes**

**Interpretation:** No gene shows a statistically concordant KA vs PBS response across all three cell types in CA3 Ipsi. This confirms the lack of power in direct KA vs PBS comparisons (all 3 comparisons have 0 FDR-significant genes). Individual cell types show some nominal hits, but they do not overlap. The treatment effect, while biologically real, is heterogeneous across cell types and/or too small relative to mouse-to-mouse variability to reach even nominal significance simultaneously in all three.

**Alternative approach:** Consider deepening A01/A02 to also run `Trt_CA3_Ipsi_Micro` at a wide threshold (p < 0.10) and ask which genes are nominally sig in ANY two of three cell types.

---

## A07 — Lost CA3 Identity: PBS-Only Regional Genes {#a07}

**Question:** Which genes define normal CA3 vs CA1 identity in healthy (PBS) animals, but LOSE this regional distinction after KA seizures?  
**Strategy:** FDR < 0.05 AND |logFC| ≥ 1.0 in `Reg_PBS_Neuron`, but NOT FDR < 0.05 in `Reg_KA_Neuron`

**Context:**
- FDR-sig CA3/CA1 genes in PBS neurons: **211**
- FDR-sig CA3/CA1 genes in KA neurons: **390**
- PBS-only (lost after KA): **13**
- KA-gained (not in PBS): ~192 (see INTERPRETATION.md)

**The 13 lost-identity genes:**

| Gene | logFC PBS | FDR PBS | logFC KA | Note |
|------|-----------|---------|----------|------|
| Septin6 | -1.23 | 0.0010 | -0.86 | CA1>CA3 normally; gradient lost in KA |
| Cnih3 | -1.48 | 0.0010 | -0.97 | AMPA receptor auxiliary subunit — CA1 |
| Igsf21 | +1.00 | 0.0010 | +0.94 | Ig superfamily — CA3 |
| Atrnl1 | +1.11 | 0.0029 | +0.79 | Attractin-like — CA3 |
| Foxo1 | +1.65 | 0.0056 | +0.82 | FOXO1 transcription factor |
| Resp18 | +1.22 | 0.0078 | +0.55 | Regulated endocrine-specific protein |
| Ccdc136 | +1.16 | 0.0084 | +0.74 | Coiled-coil domain |
| Crip2 | -1.10 | 0.0093 | -0.85 | LIM domain protein — CA1 |
| Igfbp3 | +1.10 | 0.0140 | +0.66 | IGF binding protein-3 — pro-survival |
| Epha10 | -1.09 | 0.0155 | -0.84 | Ephrin receptor — CA1 |
| Gnas | +1.04 | 0.0173 | +0.97 | Gs alpha — cAMP signaling |
| Iqsec3 | +1.09 | 0.0216 | +0.88 | IQ motif, synaptic |
| Tmem25 | +1.12 | 0.0292 | +0.97 | Transmembrane protein |

**Interpretation:** Notably, most of these genes REMAIN in the same direction in KA — they simply lose FDR significance because the KA logFC is smaller. This pattern suggests partial erosion of regional identity in KA rather than complete reversal. `Cnih3` (AMPA receptor auxiliary subunit) normally enriched in CA1 is particularly interesting — loss of this CA1/CA3 difference may reflect redistribution of AMPA receptors in the epileptic hippocampus. `Foxo1` loss of regional distinction is notable given its roles in neuronal survival and apoptosis.

---

## Prior Targeted Analysis (Reference)

| ID | Script | Comparison 1 | Comparison 2 | Genes Found | Output |
|----|--------|-------------|-------------|------------|--------|
| intersect_micro | `scripts/intersect_CA3_ipsi_micro.R` | `Side_KA_CA3_Micro` | `Trt_CA3_Ipsi_Micro` | **4 (9 at p<0.10)** | `results/KA_CA3_Ipsi_Micro_intersection.csv` |

**Microglia-specific KA CA3 Ipsi genes:** Id4, Inpp4b, Nfe2, Abhd5 (p<0.05, |logFC|≥0.5)

---

## Outstanding Follow-Up Analyses

| Priority | Analysis | Rationale |
|----------|----------|-----------|
| HIGH | A06 variant: 2-of-3 cell types | 0 genes in all 3; 2-of-3 may reveal shared KA response |
| HIGH | A03 shared across cell types | Which KA-amplified regional genes appear in 2+ cell types? |
| MED | A04 for astrocytes | Do astrocytes also show contralateral compensation? |
| MED | A07 for microglia and astrocytes | Is loss of regional identity cell-type specific? |
| LOW | A05 for astrocytes and microglia | Interaction cross-validation in non-neuronal cells |

---

## File Index

```
results/targeted/
├── A01_CA3_Ipsi_Neuron_KA_specific.csv           (1 gene)
├── A02_CA3_Ipsi_Astro_KA_specific.csv            (66 genes)
├── A03a_KA_amplified_regional_Neuron.csv         (148 genes)
├── A03b_KA_amplified_regional_Micro.csv          (39 genes)
├── A03c_KA_amplified_regional_Astro.csv          (57 genes)
├── A04_Contralateral_protective_Neuron_CA3.csv   (20 genes)
├── A04_Contralateral_protective_Neuron_strict.csv (0 genes, both regions)
├── A05_Interaction_Neuron_cross_validated.csv    (54 genes)
├── A06_Pan_celltype_KA_consensus_CA3_Ipsi.csv    (0 genes)
└── A07_Lost_CA3_identity_PBS_only_Neuron.csv     (13 genes)

results/
└── KA_CA3_Ipsi_Micro_intersection.csv            (9 genes at p<0.10)
```

---

*Last updated: targeted_analyses.R run*
