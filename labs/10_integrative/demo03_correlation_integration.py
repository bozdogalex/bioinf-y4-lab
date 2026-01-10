"""
Demo 3 — SNP–Gene correlation
"""

import pandas as pd

df_joint = pd.read_csv("multiomics_concat_demo.csv", index_col=0)

# split synthetic dataset (first 50 = SNPs, rest = genes)
snp_rows = df_joint.index[:50]
gene_rows = df_joint.index[50:]

df_snp = df_joint.loc[snp_rows]
df_gene = df_joint.loc[gene_rows]

corr = df_snp.T.corrwith(df_gene.T, axis=0)

print(corr.head())
