"""
Exercise 10.2 — Identify top SNP–Gene correlations

TODO:
- încărcați matricea integrată multi-omics
- împărțiți rândurile în SNPs vs gene (după indice sau după nume)
- calculați corelații între fiecare SNP și fiecare genă
- filtrați |r| > 0.5
- exportați snp_gene_pairs_<handle>.csv
"""

from pathlib import Path
import pandas as pd
from scipy import stats

HANDLE = "numipasaaa"
JOINT_CSV = Path(f"data/work/{HANDLE}/lab10/multiomics_concat_{HANDLE}.csv")

OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

# Load matrix and identify SNP/Gene columns (features); rows are samples (e.g., S010)
df_joint = pd.read_csv(JOINT_CSV, index_col=0)

# Detect columns
snp_cols = [c for c in df_joint.columns if str(c).startswith('rs')]
gene_cols = [c for c in df_joint.columns if str(c).startswith('Gene')]

correlations = []

for snp_col in snp_cols:
    x = pd.to_numeric(df_joint[snp_col], errors='coerce')
    for gene_col in gene_cols:
        y = pd.to_numeric(df_joint[gene_col], errors='coerce')

        mask = x.notna() & y.notna()
        x_valid = x[mask]
        y_valid = y[mask]

        if len(x_valid) < 2 or len(y_valid) < 2:
            continue
        if x_valid.nunique() < 2 or y_valid.nunique() < 2:
            continue

        r, p_value = stats.pearsonr(x_valid.values, y_valid.values)

        if abs(r) > 0.5:
            correlations.append({
                'SNP': snp_col,
                'Gene': gene_col,
                'Correlation': r,
                'P_value': p_value
            })

# Create DataFrame with SNP-gene pairs
snp_gene_pairs = pd.DataFrame(correlations)

# Check if any correlations were found
if len(snp_gene_pairs) > 0:
    # Sort by absolute correlation value
    snp_gene_pairs = snp_gene_pairs.reindex(
        snp_gene_pairs['Correlation'].abs().sort_values(ascending=False).index
    )
    
    # Save the result
    snp_gene_pairs.to_csv(OUT_CSV, index=False)
    
    print(f"Found {len(snp_gene_pairs)} SNP-gene pairs with |r| > 0.5")
    print(f"\nTop 10 correlations:")
    print(snp_gene_pairs.head(10))
else:
    print("No SNP-gene pairs found with |r| > 0.5")
    print("Consider lowering the correlation threshold.")       