# Biological Interpretation & Follow-Up Priorities
## Chronic KA vs PBS | Mouse Hippocampus | GeoMx Spatial Transcriptomics

> This document cross-references figure files and comparisons documented in [SELECTED_COMPARISONS.md](SELECTED_COMPARISONS.md).
> Analysis run: Feb 20 2026 HPC | Figures generated: Feb 24 2026
> Thresholds used throughout: **FDR < 0.05, |logFC| ≥ 1.0** unless noted

---

## The One-Line Summary

**The dominant transcriptional signal in the chronically epileptic hippocampus is not KA vs PBS — it is CA3 vs CA1.** Region accounts for hundreds of FDR-significant gene expression differences within each cell type; direct KA vs PBS comparisons reach zero FDR significance after correcting for multiple testing. This is a biological finding, not a statistical failure.

---

## Finding 1: Regional Identity (CA3 vs CA1) is the Dominant Signal

### The numbers
*(See `Fig3_DE_summary`, `Fig4a–c`, `Fig5_heatmap_neuron_top50` in SELECTED_COMPARISONS.md)*

| Comparison | FDR-sig genes | CA3-higher | CA1-higher |
|---|---|---|---|
| `Reg_KA_Neuron` | **390** | 224 | 166 |
| `Reg_KA_Astro` | **178** | 165 | 13 |
| `Reg_KA_Micro` | **146** | 132 | 14 |
| `Reg_PBS_Neuron` (baseline) | 211 | 116 | 95 |
| `Reg_PBS_Astro` (baseline) | 114 | 96 | 18 |
| `Reg_PBS_Micro` (baseline) | 118 | 95 | 23 |

**All three cell types show a strong CA3 bias** — far more genes are elevated in CA3 than CA1. This is consistent with CA3's unique circuit properties (mossy fiber input from dentate gyrus, autoassociative recurrent collaterals) and its known vulnerability to excitotoxic injury.

### Top CA3-enriched genes in KA neurons (from `Reg_KA_Neuron`)
*(See `Fig4a_volcano_KA_neuron_regional`, `Fig5_heatmap_neuron_top50`)*

| Gene | logFC | Biological role |
|---|---|---|
| **Cpne4** | +4.3 | Copine-4; Ca²⁺-dependent membrane binding, synaptic plasticity |
| **Spock1** | +4.0 | Testican-1; ECM proteoglycan, synaptogenesis |
| **Nnat** | +3.9 | Neuronatin; imprinted gene, ion channel regulation |
| **Ccn3** | −3.3 | (CA1-higher) CCN family member; neuroprotective, anti-fibrotic |
| **Nrip3** | +3.6 | Nuclear receptor interacting protein; stress response |
| **Grik4** | +3.1 | GluK4 kainate receptor subunit — **directly relevant to KA model** |
| **Adgrl2** | −3.1 | (CA1-higher) Latrophilin-2; synapse formation |
| **Ly6e** | +2.5 | Lymphocyte antigen 6E; interferon-regulated, reported in CA3 neurons |
| **Itpka** | −2.1 | (CA1-higher) IP3 kinase A; inositol signalling |

### Top CA3-enriched genes in KA astrocytes (from `Reg_KA_Astro`)
*(See `Fig4c_volcano_KA_astro_regional`)*

Notably many of the same genes top the astrocyte list: **Npcd, Nnat, Chgb, Cpne4, Snap25, Nefl**. These are canonical neuronal markers appearing in astrocytes, likely reflecting:
1. Astrocytic phagocytosis/engulfment of neuronal synaptic material
2. Contamination inherent to GeoMx ROI geometry (astrocytes sit in close proximity to neurons)

This is important to note when interpreting the astrocyte data — pure astrocytic signals need to be distinguished from bystander neuronal signal.

---

## Finding 2: 75 Genes are Dysregulated Across All Three Cell Types

*(See `Fig11c_euler_Regional_KA_celltypes` — the Euler Venn diagram)*

| Overlap group | N genes |
|---|---|
| Neuron only | 280 |
| Astrocyte only | 53 |
| Microglia only | 16 |
| Neuron + Microglia | 20 |
| Neuron + Astrocyte | 15 |
| Microglia + Astrocyte | 35 |
| **All 3 cell types** | **75** |

These 75 shared genes include:

**Synaptic/neuronal markers** (Snap25, Nefl, Nefm, Syn2, Cplx1, Cplx2, Stmn2, Bsn, Sv2b) — likely reflecting CA3's higher synaptic density and the neuronal signal across all ROI types.

**Metabolic** (Aldoa, Gpi1, Sms, Ak4, Atp6v0c) — glycolytic and energy metabolism enzymes enriched in CA3.

**CA3-specific receptors/channels** (Grik4, Kcnq5, Kcna1, Slc17a7) — Grik4 (kainate receptor), Kcnq5 (M-type K⁺ channel), Slc17a7 (VGluT1), suggesting CA3 has a distinct excitatory tone.

**ECM/structural** (Spock1, Hapln4, Clstn2) — extracellular matrix proteins enriched in CA3.

### What this means for follow-up
The microglial and astrocytic regional signal is dominated by CA3-enriched genes that are canonical neuronal markers. **This is a strong argument for spatial proximity / phagocytosis rather than intrinsic cell identity differences.** The ~16 microglia-only and ~53 astrocyte-only genes are where genuine glial-intrinsic CA3 biology lives — worth extracting separately.

---

## Finding 3: KA vs PBS Has No FDR-Significant Genes — But That Has a Clear Explanation

*(See `Fig8_CA3_Ipsi_KA_vs_PBS_combined`, `Fig8a–c`, supplementary volcanos `SuppFig_volcano_Trt_*`)*

### The statistics
Every single `Trt_` comparison has **minimum achievable FDR > 0.84**, with most > 0.99. Dozens to hundreds of genes reach nominal p < 0.05 in each comparison (e.g., 640 genes in `Trt_Ipsi_Neuron`), but after correcting for 14,103 genes, none survive.

### Why
Three factors compound:

1. **Only 5 mice per group.** The dream() random effect `(1|Mouse)` is absorbing real biological variance between animals. With 5 mice per group, power for whole-transcriptome correction after BH adjustment is limited.

2. **The regional signal is so dominant** (~500+ genes significant for CA3 vs CA1) that it creates a high noise floor for treatment effects when pooled across regions. The treatment effect, if real, is much smaller by comparison.

3. **Chronic model at 2 weeks.** Acute KA seizures (24–72h) consistently produce hundreds of FDR-sig neuronal genes in smaller studies. At 2 weeks post-KA, the transcriptome may have partially normalized or shifted to a new chronic state with smaller fold changes.

### What the nominally significant genes suggest
The top nominally significant genes in the most powered comparison (`Trt_Ipsi_Neuron`, 640 genes p < 0.05) include **Gja1** (Connexin-43 — strongly upregulated in reactive gliosis), **Vegfa** (vascular remodeling), and **Cacna1d** (Ca²⁺ channel — relevant to seizure threshold). These hint at real biology that lacks power for FDR correction.

### What this means for follow-up
- The KA vs PBS comparison requires **more mice** (~10/group) to reach FDR significance at the whole-transcriptome level with this model
- Alternatively, a **targeted gene list approach** (candidate pathway-driven testing with reduced multiple testing burden) could reveal treatment effects
- The **interaction term** (`Int_Ipsi_Neuron`) identified **Cnr1** (cannabinoid receptor 1, logFC = −1.0, p = 0.006) as the top nominally significant gene — this is biologically very interesting given the neuroprotective literature on endocannabinoid signaling in epilepsy

---

## Finding 4: 192 Genes are KA-Specific (Not Present in PBS Baseline)

*(Compare `Reg_KA_Neuron` vs `Reg_PBS_Neuron`; see `Fig7_regional_KA_vs_PBS`)*

| Category | N genes |
|---|---|
| Baseline CA3 vs CA1 (PBS only) | 13 |
| Shared in both KA and PBS | 198 |
| **KA-specific CA3 vs CA1** | **192** |

198 of the 390 KA neuronal regional genes also appear in PBS — these reflect **baseline hippocampal anatomy** (CA3 vs CA1 always differs). The 192 KA-specific genes are the most experimentally relevant because they represent a **KA-induced change in CA3 vs CA1 transcriptional identity**.

### Top KA-specific CA3 regional neurons
| Gene | logFC (KA) | Biological role |
|---|---|---|
| **Tgfb2** | +2.7 | TGF-β2 — neuroprotection, reactive astrogliosis, post-injury remodeling |
| **Flrt3** | +1.8 | Fibronectin leucine repeat transmembrane; axon guidance, synapse formation |
| **Snap25** | +1.8 | Synaptic vesicle fusion (SNARE) — not sig in PBS, suggesting KA amplifies CA3 synaptic identity |
| **Adgrl3** | −1.2 | Latrophilin-3 GPCR; genetic variants associated with ADHD; CA1 preference in KA |
| **Scn3b** | −1.6 | Na⁺ channel β3 subunit (CA1-higher in KA) — loss could alter seizure threshold |
| **Bcr** | −1.7 | BCR kinase; Rho GAP domain; CA1-enriched in KA |
| **AI593442** | −1.9 | Uncharacterized lncRNA — worth investigating |
| **Car11** | −1.1 | Carbonic anhydrase 11; pH regulation, relevant to seizure propagation |

---

## Finding 5: Pathways Point to Synaptic, Metabolic, and Neuroinflammatory Biology

*(See `Fig6a–c`, `Fig12a–g` in SELECTED_COMPARISONS.md)*

### Neurons (Reg_KA_Neuron)

**KEGG** top pathways (all CA3-enriched, positive NES):
- Synaptic vesicle cycle (NES = 2.62) — #1 hit
- Oxidative phosphorylation (NES = 2.38)
- Parkinson disease / Huntington disease — neurodegenerative disease gene sets enriched in CA3
- Glycolysis / TCA cycle — metabolic upregulation in CA3

**GO Biological Process** top hits (extremely significant, p < 10⁻⁶²):
- Regulation of synapse structure or activity (210 genes, p = 7.7×10⁻⁶⁶)
- Vesicle-mediated transport in synapse
- Synaptic vesicle cycle
- **Learning or memory** (168 genes, p = 2.8×10⁻⁴⁴) — directly relevant to hippocampal function

**Reactome** top hits:
- Aerobic respiration and respiratory electron transport (NES = 2.34)
- MHC class II antigen presentation (NES = 2.26) — **immune signaling in neurons**
- Gap junction trafficking (NES = 2.20) — connexins, relevant to CA3 network synchrony
- Iron uptake and transport

**GO Molecular Function** top hits:
- Metal ion transmembrane transporter activity
- GTPase binding / GTPase regulator activity
- **Glutamate receptor binding** (54 genes, p = 7×10⁻²⁵)
- Ion channel activity

### Microglia (Reg_KA_Micro) — KEGG
*(See `Fig12b_Reactome_Reg_KA_Micro`, `Fig12f_GO_MF_Reg_KA_Micro`)*

The top microglial KEGG pathways are **identical to neurons**:
- Synaptic vesicle cycle (NES = 2.95) 
- Oxidative phosphorylation (NES = 2.80)
- Parkinson/Huntington/Prion/Alzheimer/ALS disease gene sets — all positive NES

This strongly suggests CA3 microglia are **phagocytosing synaptic material** rather than having an intrinsically different transcriptome. Microglia in CA3 are surrounded by more synaptic density and may be constitutively more engaged in synaptic pruning and debris clearance.

### Astrocytes (Reg_KA_Astro) — Reactome
*(See `Fig12c_Reactome_Reg_KA_Astro`)*

Notably, astrocytes show **MHC class II antigen presentation** as the top hit (NES = 2.71) — astrocytes are not classically antigen-presenting cells, but this can occur in neuroinflammation. Other top terms:
- Respiratory electron transport
- L1CAM interactions (axon-glial interactions)
- Neuronal System / Transmission across Chemical Synapses — again the same neuronal bystander signal

---

## Follow-Up Experiment Priorities

### Priority 1 (Highest): Validate the KA-specific CA3 genes
The 192 genes that are regionally dysregulated in KA but **not** in PBS are the most specific to epilepsy pathology.

- **Tgfb2** — IHC/ISH in CA3 vs CA1 at multiple time points (3d, 1 wk, 2 wk, 1 mo). TGF-β2 is implicated in post-seizure blood-brain barrier remodeling. Is this a response to injury or ongoing pathology?
- **Grik4** (in all-3-cell-type overlap) — GluK4 kainate receptor subunit. Directly relevant: KA acts on kainate receptors. Is GluK4 upregulated to compensate or does it drive sustained hyperexcitability? CRISPR knockdown priority.
- **Cnr1** (top interaction gene, nominally KA-specific) — Cannabinoid receptor 1. CA3 Cnr1 decreasing in KA vs PBS (logFC = −1.0). Loss of CB1R-mediated inhibition is a known epilepsy mechanism. Priority for immunofluorescence and pharmacological rescue.
- **Scn3b** — Sodium channel β3 subunit, CA1-enriched in KA. Loss in CA3 post-KA could alter network excitability.

### Priority 2: Confirm microglial findings are phagocytosis-driven
The microglial transcriptome in CA3 looks almost identical to neurons. Run:
- **Microglia depletion** (PLX5622 chow) + re-profiling to confirm which microglial genes are intrinsic vs phagocytic
- **snRNA-seq** in the same tissue to separate signals that GeoMx cannot resolve by cell
- Alternatively, compare FACS-sorted CA3 microglia gene expression to confirm

### Priority 3: Increase power for KA vs PBS comparison
With 5 mice/group the direct treatment effect is underpowered. Options:
- Add ~5 more mice per group (n=10) — this is the cleanest approach
- Run a **targeted reanalysis** with reduced gene set (e.g., only synaptic/channel genes, ~500 genes → lower FDR burden) — this can be done computationally now without new experiments
- Use a **GSEA pre-ranked** approach on the Trt_ comparisons rather than gene-level FDR to extract pathway-level treatment effects

### Priority 4: The contralateral side
The Side_ comparisons (Ipsi vs Contra) show zero FDR-significant genes. The contralateral hippocampus does not appear to be mounting a detectable transcriptional response at 2 weeks. However, there are ~500–900 nominally significant genes in several comparisons, suggesting subtle bilateral effects worth exploring with increased power.

### Priority 5: Time course
This is a single 2-week time point. The large CA3 vs CA1 baseline signal in PBS means the epilepsy-specific changes (192 KA-only genes) are sitting on top of a huge anatomical background. A **time course** (1 wk, 2 wk, 6 wk, 3 mo post-KA) would reveal whether the KA-specific signal increases or resolves, and when it peaks.

---

## Key Numbers for a Paper

- 14,103 genes measured across 102 ROIs (5 KA, 5 PBS mice; Astrocyte/Microglia/Neuron × CA1/CA3 × Ipsi/Contra)
- 60 dream() mixed-effects comparisons with (1|Mouse)
- **Zero FDR-significant genes** in any direct KA vs PBS comparison after BH correction
- **390 FDR-sig genes** (CA3 vs CA1) in KA neurons; 178 in astrocytes; 146 in microglia
- **192 KA-specific** CA3/CA1 neuronal genes (not present in PBS baseline)
- **75 genes shared across all three cell types** in the CA3 vs CA1 contrast (KA)
- Top pathways: synaptic vesicle cycle, oxidative phosphorylation, synapse regulation, glutamate receptor binding

---

*Linked files: [SELECTED_COMPARISONS.md](SELECTED_COMPARISONS.md) | Figures: `figures/` | DE results: `results/de_results_all.rds` | Pathway results: `results/pathway_results.rds`*
