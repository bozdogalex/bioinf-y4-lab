## Note
Aceste fișiere contin date **reale**, dar de dimensiune redusă pentru demonstrații.
- FASTA: secvențe TP53 (om + ortologi)
- Proteină TP53 (UniProt)
- Subset de expresie genică (GEO)

**Pentru laboratoare:** descărcați date reale extinse din NCBI/GEO/Ensembl/TCGA, conform instrucțiunilor din README-urile laboratoarelor.

### Expression Matrix + Metadata

- **`expr_small_tp53_geo.csv`** – subset of expression values for 8 genes across 8 samples.  
- **`expr_small_tp53_geo_metadata.csv`** – companion file giving sample-level labels.

Example:

```csv
SampleID,Condition
Sample1,Tumor
Sample2,Tumor
Sample3,Tumor
Sample4,Normal
```
How to use:

Merge the expression table with the metadata by SampleID.

Use Condition as the label column for:

Clustering → compare Tumor vs Normal grouping.

Supervised ML (RF, SVM, FL labs) → train/predict tumor status.

This file is intentionally small and balanced (toy example).
For assignments, students must fetch and annotate real GEO/TCGA metadata.
