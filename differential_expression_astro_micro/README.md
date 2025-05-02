### Differential Gene Expression on ASTROCYTES AND MICROGLIA

## Data
obj.RData from /pipseq/data which is from line 0 to line 212 from Ellie's notebook

## Steps
- Used Ellie's code for assign cluster
- Subset microglia and astrocytes to only include those with fewer than 3000 detected genes
- Normalize microglia and astrocytes
- Compares astrocytes vs. microglia to identify differentially expressed genes

## Output
result csv
`differential_expression_astro_micro/astro_vs_microglia_DE.csv`
`differential_expression_astro_micro/astro_micro_de_sorted.csv`
