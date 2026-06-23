# Analysis of non-coding genome parallel evolution of Hymenoptera eusociality

A comparative genomic analysis of repeat element content across 186 Hymenoptera genomes in relation to social behaviour, using phylogenetic comparative methods.

---

## Project Overview

This project investigates whether the evolution of sociality in Hymenoptera (bees, wasps, ants, and sawflies) is associated with changes in transposable element (TE) content. We integrate RepeatMasker annotation data from 186 genomes with sociality classifications and a taxonomic phylogeny to test for associations using both standard and phylogenetically informed statistical approaches.

**Key question:** Do social and non-social Hymenoptera differ in their repeat element composition, and if so, does this reflect genuine convergent genomic evolution or phylogenetic confounding?

---

## Repository Structure

```
├── scripts/
│   ├── 01_hymenoptera_tree.R             # Tree + heatmap visualisation
│   ├── 02_hymenoptera_wilcoxon.R         # Wilcoxon test + ACE transitions
│   ├── 03_hymenoptera_wilcoxon_plot.R    # Boxplots + results table
│   ├── 04_hymenoptera_pgls.R             # pGLS + BRUNCH + phylo.d
│   ├── 05_hymenoptera_blomberg.R         # Blomberg's K phylogenetic signal
│   └── 06_hymenoptera_final_table.R      # Summary table + bubble plot
│
└── README.md
```

> **Note:** Scripts must be run in numbered order (01 → 06), as later scripts depend on objects created by earlier ones within the same R session.

---

## Sociality Categories

Raw annotations were harmonised into six categories:

| Category | Binary group | n |
|---|---|---|
| Eusocial | Social | 45 |
| Primitive (primitively eusocial) | Social | 20 |
| Partially social | Social | 2 |
| Solitary | Non-social | 97 |
| Kleptoparasite | Non-social | 10 |
| Unknown | excluded | — |

---

## Methods Summary

### Step 1 — Tree Construction and Visualisation (`01_hymenoptera_tree.R`)

- Taxonomic tree built from a six-level hierarchy (suborder → superfamily → family → subfamily → genus → species) using `ape::as.phylo()`
- Degree-2 nodes collapsed with `ape::collapse.singles()`; polytomies resolved with `ape::multi2di(random = TRUE, seed = 42)`
- Tree visualised with `ggtree` with tip points coloured by sociality category
- Heatmap of 17 TE classes (Unclassified excluded; values min-max normalised per column) attached via `gheatmap()`

---

### Step 2 — Ancestral State Reconstruction + Wilcoxon Test (`02_hymenoptera_wilcoxon.R`)

- Ancestral states of Social/Non-social estimated by Maximum Likelihood ACE (equal-rates model) using `ape::ace()`
- Evolutionary transitions identified by comparing each internal node to its parent node
- Two-sample Wilcoxon rank-sum test for each of 17 TE classes (Social vs Non-social)
- Multiple testing correction: Benjamini–Hochberg FDR across all 17 tests

**Key result:** 6 TE classes significant at p_adj < 0.05 (Small RNA, L2/CR1/Rex, Gypsy/DIRS1, Rolling-circles higher in Non-social; Low complexity, Simple repeats higher in Social). 8 independent Non-social → Social transitions inferred.

---

### Step 3 — Wilcoxon Visualisation (`03_hymenoptera_wilcoxon_plot.R`)

- Boxplots for significant TE classes only, with jittered individual points
- Summary results table rendered as PDF
- Depends on objects from script 02

---

### Step 4 — Phylogenetic Comparative Methods (`04_hymenoptera_pgls.R`)

Three methods implemented via the `caper` package:

**phylo.d** — phylogenetic signal of the binary Social/Non-social trait  
- D = −0.327 (p < 0.001 vs random; p = 0.955 vs Brownian)  
- Sociality is strongly phylogenetically conserved (more clustered than Brownian motion)

**pGLS** — phylogenetic generalised least squares regression  
- Model: `log1p(TE) ~ Social_binary` with Pagel's lambda estimated by ML  
- Pagel's lambda 0.47–1.00 for most TE classes → strong phylogenetic structuring  
- No TE class significant after BH correction

**BRUNCH** — independent clade contrasts at Social/Non-social transition nodes  
- 8 independent contrasts available  
- No TE class significant after BH correction

---

### Step 5 — Blomberg's K (`05_hymenoptera_blomberg.R`)

- Blomberg's K estimated for each TE class using `phytools::phylosig()` with 1,000 permutations on log1p-transformed values
- K ranges from 0.13 to 0.65 across TE classes — all below 1 (faster than Brownian)
- No TE class significant after BH correction
- Results integrated with pGLS lambda and Wilcoxon p-values into a unified interpretation table

**Key finding:** The discordance between non-significant K and high pGLS lambda is methodologically expected — K has lower power at this sample size, while lambda captures residual covariance conservatively.

---

### Step 6 — Final Summary (`06_hymenoptera_final_table.R`)

- Bubble plot: Wilcoxon significance (x-axis) vs Pagel's lambda (y-axis), point size = Blomberg's K
  - Red zone (upper right): Wilcoxon significant + high lambda → phylogenetic confounding
  - Blue zone (lower right): Wilcoxon significant + low lambda → potentially real signal
- All three method results combined in a single exportable table

---

## Dependencies

All analyses were performed in **R**. Required packages:

```r
install.packages(c(
  "readxl",    # reading Excel input
  "dplyr",     # data manipulation
  "stringr",   # string operations
  "tidyr",     # reshaping data
  "ape",       # phylogenetic tree operations
  "ggtree",    # tree visualisation (Bioconductor)
  "ggplot2",   # plotting
  "ggnewscale",# multiple colour scales
  "cowplot",   # combining plots
  "caper",     # pGLS, BRUNCH, phylo.d
  "phytools",  # Blomberg's K, ACE utilities
  "ggrepel",   # non-overlapping labels
  "gridExtra", # table rendering
  "grid",      # grid.draw()
  "scales"     # axis formatting
))

```

---

## Notes

- **Unclassified TE category** is set to zero before analysis. In the Extended results sheet, Unclassified values reach up to 75% of some genomes — an artefact of the RepeatMasker run compared to near-zero values in the Partial results sheet. This category is retained as a column in the heatmap (rendered white) but excluded from all statistical tests.
- **Branch lengths** are set uniformly to 1 (topological tree). A time-calibrated molecular phylogeny would improve accuracy of all phylogenetic methods.
- **Polytomies** within genera (e.g., *Bombus*, *Polistes*, *Apis*) are resolved randomly with `set.seed(42)` for reproducibility.
- **Kleptoparasites** are classified as Non-social: they lack cooperative brood care and colony structure despite being evolutionarily derived from social lineages.
