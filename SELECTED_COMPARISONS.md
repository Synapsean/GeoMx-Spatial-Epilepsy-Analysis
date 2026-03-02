# Selected Comparisons Summary

---

## 1. TREATMENT EFFECTS (KA vs PBS) - 18 comparisons

### A. Highly Specific (All 12)
- **CA3 Ipsi**: Astrocyte, Microglia, Neuron (3)
- **CA3 Contra**: Astrocyte, Microglia, Neuron (3)
- **CA1 Ipsi**: Astrocyte, Microglia, Neuron (3)
- **CA1 Contra**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Most sensitive to detect treatment effects in specific anatomical locations

### C. Pooled Across Regions (6)
- **Ipsi**: Astrocyte, Microglia, Neuron (3)
- **Contra**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Higher power to detect side-specific treatment effects regardless of region

---

## 2. REGIONAL EFFECTS (CA3 vs CA1) - 18 comparisons

### A. Within Treatment + Side (All 12)
- **KA Ipsi**: Astrocyte, Microglia, Neuron (3)
- **KA Contra**: Astrocyte, Microglia, Neuron (3)
- **PBS Ipsi**: Astrocyte, Microglia, Neuron (3)
- **PBS Contra**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Test CA3 vulnerability in all conditions; important for understanding if regional differences exist at baseline (PBS) vs after injury (KA), and if they differ by side

### B. Within Treatment Pooled Across Sides (6)
- **KA**: Astrocyte, Microglia, Neuron (3)
- **PBS**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Higher power to detect overall regional differences within each treatment group

---

## 3. SIDE EFFECTS (Ipsi vs Contra) - 18 comparisons

### A. Within Treatment + Region (All 12)
- **KA CA3**: Astrocyte, Microglia, Neuron (3)
- **KA CA1**: Astrocyte, Microglia, Neuron (3)
- **PBS CA3**: Astrocyte, Microglia, Neuron (3)
- **PBS CA1**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Test your hypothesis about contralateral "protective" mechanisms; detect compensatory responses in less-damaged hemisphere

### B. Within Treatment Pooled Across Regions (6)
- **KA**: Astrocyte, Microglia, Neuron (3)
- **PBS**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Higher power to detect overall side differences within each treatment, regardless of region

---

## 4. INTERACTION EFFECTS - 6 comparisons

### A. Region × Treatment (Ipsi and Contra separately)
- **Ipsi**: Astrocyte, Microglia, Neuron (3)
- **Contra**: Astrocyte, Microglia, Neuron (3)

**Rationale**: Test if CA3 vulnerability changes with KA treatment, separately for each hemisphere; addresses whether regional differences are injury-induced or baseline anatomical

---

## IMPLEMENTATION NOTES

**Dream() formulas for each comparison type:**

1. **Treatment Effects (highly specific)**: `~ Group + (1|Mouse)`
   - Example: CA3 Ipsi Astrocytes only

2. **Treatment Effects (pooled regions)**: `~ Group + Region + (1|Mouse)`
   - Example: Ipsi Astrocytes (CA1 + CA3), adjust for Region

3. **Regional Effects (within treatment+side)**: `~ Region + (1|Mouse)`
   - Example: KA Ipsi Astrocytes, compare CA3 vs CA1

4. **Regional Effects (pooled sides)**: `~ Region + Side + (1|Mouse)`
   - Example: KA Astrocytes, compare CA3 vs CA1, adjust for Side

5. **Side Effects (within treatment+region)**: `~ Side + (1|Mouse)`
   - Example: KA CA3 Astrocytes, compare Ipsi vs Contra

6. **Side Effects (pooled regions)**: `~ Side + Region + (1|Mouse)`
   - Example: KA Astrocytes, compare Ipsi vs Contra, adjust for Region

7. **Interactions (Region × Treatment)**: `~ Region * Group + (1|Mouse)`
   - Example: Ipsi Astrocytes only, test Region:Group interaction term

---

## STATISTICS CONSIDERATIONS

- **Total tests**: 60 comparisons
- **Multiple testing**: With FDR < 0.10, expect some false positives
- **Sample sizes per comparison**:
  - Highly specific: ~10-15 ROIs from 5 mice per group
  - Pooled: ~20-30 ROIs from 5 mice per group


---

### Fig 3 — DE Summary
**File:** `Fig3_DE_summary`
**Shows:** Horizontal bar chart of the number of DE genes per comparison at two thresholds (FDR < 0.05 and nominal p < 0.05). Key takeaway: **Reg_** (CA3 vs CA1) comparisons dominate; **Trt_** (KA vs PBS) comparisons do not reach FDR significance.

---

### Fig 4 — Volcano Plots (Key Comparisons)

| File | Comparison | What it shows |
|------|-----------|---------------|
| `Fig4a_volcano_KA_neuron_regional` | `Reg_KA_Neuron` | CA3 vs CA1 in KA neurons (pooled sides) — most significant regional contrast |
| `Fig4b_volcano_KA_micro_regional` | `Reg_KA_Micro` | CA3 vs CA1 in KA microglia |
| `Fig4c_volcano_KA_astro_regional` | `Reg_KA_Astro` | CA3 vs CA1 in KA astrocytes |
| `Fig4d_volcano_interaction_neuron` | `Int_Ipsi_Neuron` | Region × Treatment interaction in ipsilateral neurons — genes whose CA3/CA1 difference is **uniquely driven by KA injury** |

---

### Fig 5 — Heatmap: Top 50 Neuronal Regional DE Genes
**File:** `Fig5_heatmap_neuron_top50.pdf` *(PDF only)*
**Shows:** Z-scored expression of the top 50 FDR-significant genes from `Reg_KA_Neuron` across all neuron ROIs, annotated by Treatment, Region, and Side. Reveals gene modules associated with CA3 identity in KA.

---

### Fig 6 — Pathway Enrichment: Regional Neurons (CA3 vs CA1, KA)

| File | Method | What it shows |
|------|--------|---------------|
| `Fig6a_KEGG_dotplot_neuron` | KEGG GSEA | Top 20 KEGG pathways enriched in CA3 vs CA1 neurons (KA) |
| `Fig6b_KEGG_network_neuron` | KEGG emapplot | Network map of pathway similarity for KEGG results (generated only if ≥ 5 pathways significant) |
| `Fig6c_GO_BP_dotplot_neuron` | GO Biological Process GSEA | Top 20 GO-BP terms for CA3 vs CA1 neurons (KA) |

---

### Fig 7 — Regional Effect: KA vs PBS Scatter
**File:** `Fig7_regional_KA_vs_PBS`
**Shows:** Scatter plot comparing CA3/CA1 logFC in KA (y) vs PBS (x) for each cell type. Points above the diagonal = genes with **amplified regional differences in epilepsy**. Colored by FDR significance in KA. Includes LM trendline.

---

### Fig 8 — Primary Experimental Comparison: CA3 Ipsi KA vs PBS

| File | Comparison | What it shows |
|------|-----------|---------------|
| `Fig8a_CA3_Ipsi_KA_vs_PBS_Neuron` | `Trt_CA3_Ipsi_Neuron` | Volcano: KA vs PBS in CA3 ipsilateral neurons |
| `Fig8b_CA3_Ipsi_KA_vs_PBS_Micro` | `Trt_CA3_Ipsi_Micro` | Volcano: KA vs PBS in CA3 ipsilateral microglia |
| `Fig8c_CA3_Ipsi_KA_vs_PBS_Astro` | `Trt_CA3_Ipsi_Astro` | Volcano: KA vs PBS in CA3 ipsilateral astrocytes |
| `Fig8_CA3_Ipsi_KA_vs_PBS_combined` | All 3 above | 3-panel combined figure (patchwork). **Note: no FDR-sig genes** in these direct KA vs PBS comparisons; genes labeled by nominal p-value |

---

### Fig 9 — Pathway Enrichment: CA3 Ipsi KA vs PBS

| File | Comparison | Method |
|------|-----------|--------|
| `Fig9a_KEGG_CA3_Ipsi_Neuron_treatment` | `Trt_CA3_Ipsi_Neuron` | KEGG GSEA dotplot |
| `Fig9b_KEGG_CA3_Ipsi_Micro_treatment` | `Trt_CA3_Ipsi_Micro` | KEGG GSEA dotplot |
| `Fig9c_GO_BP_CA3_Ipsi_Neuron_treatment` | `Trt_CA3_Ipsi_Neuron` | GO-BP GSEA dotplot |

*(Only generated if pathways reach significance; check for empty/missing files)*

---

### Fig 10 — Heatmap: Top 50 CA3 Ipsi KA vs PBS Genes
**File:** `Fig10_heatmap_CA3_Ipsi_KA_vs_PBS_top50.pdf` *(PDF only)*
**Shows:** Z-scored expression of top 50 genes by nominal p-value from `Trt_CA3_Ipsi_Neuron` across CA3 ipsilateral neuron ROIs, annotated by Treatment. Uses p-value ranking (not FDR) because no FDR-significant DE genes exist in this direct treatment comparison.

---

### Fig 11 — Overlap / Venn / UpSet Plots

| File | Status | Description |
|------|--------|-------------|
| `Fig11a_euler_KA_vs_PBS_celltypes_nominal` | ✅ Present | eulerr 3-set Venn of **nominally significant** (p < 0.05, \|logFC\| ≥ 1.0) KA vs PBS genes per cell type (Ipsi, pooled regions). Note: zero FDR-significant genes exist in direct KA vs PBS comparisons — this uses uncorrected p-values and is labelled as such on the figure. Counts: Neuron=6, Microglia=1, Astrocyte=3. |
| `Fig11b_UpSetR_Regional_all_celltype` | ✅ Present | UpSetR intersection plot of FDR-sig CA3 vs CA1 genes across all 6 regional comparisons (`Reg_KA/PBS × Neuron/Micro/Astro`). 522 total FDR-sig genes; shows which genes are shared across cell types and whether the regional signal is KA-specific or also present in PBS. |
| `Fig11c_euler_Regional_KA_celltypes` | ✅ Present | eulerr 3-set Venn of FDR-sig CA3 vs CA1 genes (`Reg_KA_*`) per cell type (KA animals only). Shows overlap of regionally DE genes between Neuron, Microglia, and Astrocyte. |

---

### Fig 12 — Reactome & GO Molecular Function Pathway Plots

| File | Comparison | Method | Status |
|------|-----------|--------|--------|
| `Fig12a_Reactome_Reg_KA_Neuron` | `Reg_KA_Neuron` | Reactome GSEA dotplot | ✅ Present |
| `Fig12b_Reactome_Reg_KA_Micro` | `Reg_KA_Micro` | Reactome GSEA dotplot | ✅ Present |
| `Fig12c_Reactome_Reg_KA_Astro` | `Reg_KA_Astro` | Reactome GSEA dotplot | ✅ Present |
| `Fig12d_Reactome_CA3_Ipsi_Neuron_treatment` | `Trt_CA3_Ipsi_Neuron` | Reactome GSEA dotplot | ✅ Present |
| `Fig12e_GO_MF_Reg_KA_Neuron` | `Reg_KA_Neuron` | GO-MF GSEA dotplot | ✅ Present |
| `Fig12f_GO_MF_Reg_KA_Micro` | `Reg_KA_Micro` | GO-MF GSEA dotplot | ✅ Present |
| `Fig12g_GO_MF_CA3_Ipsi_Neuron_treatment` | `Trt_CA3_Ipsi_Neuron` | GO-MF GSEA dotplot | ✅ Present |

---

## Supplementary Figures (~70 volcano plots)

One volcano plot per DE comparison (all 60), saved as `SuppFig_volcano_{comparison_name}`. Labels top genes by nominal p-value. Organized by comparison prefix:

| Prefix | N | Description |
|--------|---|-------------|
| `Trt_CA3_Ipsi/Contra_*` | 6 | KA vs PBS in CA3 × Side × Cell Type |
| `Trt_CA1_Ipsi/Contra_*` | 6 | KA vs PBS in CA1 × Side × Cell Type |
| `Trt_Ipsi/Contra_*` | 6 | KA vs PBS pooled regions, by Side × Cell Type |
| `Reg_KA/PBS_Ipsi/Contra_*` | 12 | CA3 vs CA1, within Treatment × Side × Cell Type |
| `Reg_KA/PBS_*` | 6 | CA3 vs CA1, pooled sides, within Treatment × Cell Type |
| `Side_KA/PBS_CA3/CA1_*` | 12 | Ipsi vs Contra, within Treatment × Region × Cell Type |
| `Side_KA/PBS_*` | 6 | Ipsi vs Contra, pooled regions, within Treatment × Cell Type |
| `Int_Ipsi/Contra_*` | 6 | Region × Treatment interaction, by Side × Cell Type |

---

## Outstanding / To Do

- [x] **Fig11a** — Generated as nominal p-value version (`_nominal` suffix); no FDR-sig genes in direct KA vs PBS comparison is a biological result, not a failure
- [x] **Fig11b** — Generated; 522 FDR-sig genes across 6 regional comparisons
- [ ] Review Fig3 (DE summary) — confirm category labels parsed correctly from comparison names
- [ ] Review Fig8 volcanos — no FDR-sig genes expected; confirm gene labels look reasonable at nominal threshold
- [ ] Consider: cross-comparison pathway analysis (pathways shared across Neuron/Micro/Astro in KA response)
- [ ] Consider: update Shiny dashboard to include pathway results

