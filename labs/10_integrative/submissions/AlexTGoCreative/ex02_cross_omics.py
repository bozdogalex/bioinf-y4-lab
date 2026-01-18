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
import numpy as np

HANDLE = "AlexTGoCreative"
JOINT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/multiomics_concat_{HANDLE}.csv")

OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

print("=" * 60)
print("Exercise 10.2 — Cross-Omics Marker Discovery")
print("=" * 60)

# Load joint matrix
print("\n1. Loading integrated multi-omics matrix...")
df_joint = pd.read_csv(JOINT_CSV, index_col=0)
print(f"   Joint matrix: {df_joint.shape} (features × samples)")

# Split rows into SNPs vs Genes
print("\n2. Splitting SNPs and Genes...")
# SNPs start with 'rs', genes start with 'GENE'
snp_rows = [idx for idx in df_joint.index if idx.startswith('rs')]
gene_rows = [idx for idx in df_joint.index if idx.startswith('GENE')]

df_snp = df_joint.loc[snp_rows]
df_gene = df_joint.loc[gene_rows]

print(f"   SNPs: {len(snp_rows)}")
print(f"   Genes: {len(gene_rows)}")

# Compute correlations between each SNP and each gene
print("\n3. Computing SNP–Gene correlations...")
print(f"   Total pairs to compute: {len(snp_rows) * len(gene_rows)}")

correlations = []

for snp in snp_rows:
    snp_values = df_snp.loc[snp].values
    
    for gene in gene_rows:
        gene_values = df_gene.loc[gene].values
        
        # Pearson correlation
        corr = np.corrcoef(snp_values, gene_values)[0, 1]
        
        correlations.append({
            'SNP': snp,
            'Gene': gene,
            'Correlation': corr,
            'Abs_Correlation': abs(corr)
        })

df_corr = pd.DataFrame(correlations)

print(f"   Computed {len(df_corr)} correlations")
print(f"   Mean |r|: {df_corr['Abs_Correlation'].mean():.4f}")
print(f"   Max |r|: {df_corr['Abs_Correlation'].max():.4f}")

# Filter by |r| > 0.5
print("\n4. Filtering high correlations (|r| > 0.5)...")
threshold = 0.5
df_filtered = df_corr[df_corr['Abs_Correlation'] > threshold].copy()
df_filtered = df_filtered.sort_values('Abs_Correlation', ascending=False)

print(f"   Pairs with |r| > {threshold}: {len(df_filtered)}")

if len(df_filtered) > 0:
    print(f"\n   Top 10 SNP–Gene pairs:")
    print(df_filtered[['SNP', 'Gene', 'Correlation']].head(10).to_string(index=False))
else:
    print(f"   No pairs found with |r| > {threshold}")
    print(f"   Lowering threshold to 0.3 for demonstration...")
    threshold = 0.3
    df_filtered = df_corr[df_corr['Abs_Correlation'] > threshold].copy()
    df_filtered = df_filtered.sort_values('Abs_Correlation', ascending=False)
    print(f"   Pairs with |r| > {threshold}: {len(df_filtered)}")
    if len(df_filtered) > 0:
        print(f"\n   Top 10 SNP–Gene pairs:")
        print(df_filtered[['SNP', 'Gene', 'Correlation']].head(10).to_string(index=False))

# Export results
print(f"\n5. Exporting results...")
df_filtered.to_csv(OUT_CSV, index=False)
print(f"   Saved: {OUT_CSV}")

# Summary statistics
print("\n" + "=" * 60)
print("✓ Exercise 10.2 Complete!")
print("=" * 60)
print(f"\nSummary:")
print(f"  • Total SNP-Gene pairs analyzed: {len(df_corr)}")
print(f"  • High-correlation pairs (|r| > {threshold}): {len(df_filtered)}")
print(f"  • Strongest correlation: {df_filtered['Abs_Correlation'].max():.4f}" if len(df_filtered) > 0 else "  • No strong correlations found")
print(f"  • Output: {OUT_CSV}")

# Additional analysis: correlation distribution
print("\nCorrelation distribution:")
bins = [-1, -0.7, -0.5, -0.3, 0, 0.3, 0.5, 0.7, 1]
labels = ['< -0.7', '-0.7 to -0.5', '-0.5 to -0.3', '-0.3 to 0', 
          '0 to 0.3', '0.3 to 0.5', '0.5 to 0.7', '> 0.7']
df_corr['Category'] = pd.cut(df_corr['Correlation'], bins=bins, labels=labels)
print(df_corr['Category'].value_counts().sort_index())

# Save full correlation matrix for potential further analysis
full_corr_csv = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_correlations_full_{HANDLE}.csv")
df_corr.to_csv(full_corr_csv, index=False)
print(f"\nFull correlation matrix also saved to: {full_corr_csv}")
