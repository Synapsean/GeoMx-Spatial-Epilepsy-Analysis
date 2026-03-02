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

