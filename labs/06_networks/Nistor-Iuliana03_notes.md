# Lab 06 - Gene Co-Expression Networks (GCEs)

## Correlation method and threshold used
- Correlation method: **Pearson**
- Absolute correlation used: **Yes**
- Variance threshold: **0.0** (kept all genes)
- Adjacency threshold: **0.2**

## Short reflection
Compared to classic clustering (Lab 5), a gene co-expression network focuses on pairwise relationships between genes.  
Instead of grouping genes only by global similarity, GCE networks keep the edges (correlation strengths) between pairs of genes.  
Modules are detected as graph communities, which reveal connectivity patterns and hub genes that classical clustering cannot capture.
