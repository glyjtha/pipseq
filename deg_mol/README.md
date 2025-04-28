# DEG Analysis for MOL Subtypes

This folder contains the results and code for the differential expression (DEG) analysis between Thrb+ and SafeTarget+ cells across three MOL subtypes:

- **MOL_12**
- **MOL_34**
- **MOL_56**

### Contents
- **CSV files**:
  - `DEG_MOL_12_Thrb_vs_ST.csv`
  - `DEG_MOL_34_Thrb_vs_ST.csv`
  - `DEG_MOL_56_Thrb_vs_ST.csv`
  - Each file contains DEG results, ranked by adjusted p-value (FDR) and log2 fold change.
- **R Markdown file**:
  - `DEG_mol.Rmd`
  - Script for running DEG analysis and generating visualizations (volcano plots, violin plots).

### Notes
- Only cells with **total guide UMI > 2** were included in the analysis.
- MOL subtypes were **reassigned** based on updated clustering:
  - Cluster 0, 4 → MOL_12
  - Cluster 5 → MOL_34
  - Cluster 1, 2, 3 → MOL_56
- Differential expression was performed using normalized RNA counts, with:
  - `logfc.threshold = 0.25`
  - `min.pct = 0.05`

### How to Interpret
- **Positive log2 fold change**: gene is upregulated in **Thrb+** cells.
- **Negative log2 fold change**: gene is upregulated in **SafeTarget+** cells.
- Top DE genes are visualized in violin plots per subtype in the `.Rmd` output.

