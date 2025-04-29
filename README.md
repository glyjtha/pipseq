# pipseq

This repository contains code, data, and analyses
---

## Repository Structure

### `deg_mol/`
- **Purpose:** Differential expression analysis for MOL subtypes.
- **Contents:** 
  - `.csv` files with DEG results
  - `.Rmd` script for running DEG
  - Volcano and violin plots per subtype
  - `README.md` inside the folder for more detail

### `data/`
- Contains all datas need to use from Ellie.

### `pseudotime_tSNE/`
- Pseudotime and dimensionality reduction (t-SNE) analyses on selected clusters.


### `relevant_code/`
- Contains Ellie's scripts for manipulating Seurat objects, adding layers, and clustering.

### `differential_expression_astro_micro/`
- **Purpose:** DEG analysis for astrocyte and microglia subsets (organized separately).


R packages installed:
```r
install.packages(c("Seurat", "ggplot2", "ggrepel", "plyr", "Matrix", "patchwork"))
