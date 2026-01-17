"""
Exercise 10.2 — Identify top SNP–Gene correlations
"""

from pathlib import Path
import pandas as pd
import numpy as np

HANDLE = "banmepls"
JOINT_CSV = Path(f"data/work/{HANDLE}/lab10/multiomics_concat_{HANDLE}.csv")

OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

CORR_THRESHOLD = 0.5
TOP_N = 1000


def main():
    if not JOINT_CSV.exists():
        raise FileNotFoundError(
            f"Fișierul joint lipsește: {JOINT_CSV}\n"
            "Rulează întâi Exercise 10.1 pentru a genera multiomics_concat_{HANDLE}.csv"
        )
    
    df_joint = pd.read_csv(JOINT_CSV, index_col=0)
    
    snp_rows = [idx for idx in df_joint.index if idx.startswith(('rs', 'SNP'))]
    gene_rows = [idx for idx in df_joint.index if idx not in snp_rows]
    
    if not snp_rows or not gene_rows:
        raise ValueError("Nu am putut separa SNP-uri și gene după nume. Verifică index-ul!")
    
    
    df_snp = df_joint.loc[snp_rows]
    df_gene = df_joint.loc[gene_rows]
    
    corr_matrix = df_snp.T.corrwith(df_gene.T, axis=0, method='pearson')
    
    df_corr = pd.DataFrame()
    for snp in df_snp.index:
        for gene in df_gene.index:
            r = df_snp.loc[snp].corr(df_gene.loc[gene])
            df_corr = pd.concat([df_corr, pd.DataFrame({
                'snp': [snp],
                'gene': [gene],
                'correlation': [r]
            })], ignore_index=True)
    
    df_filtered = df_corr[np.abs(df_corr['correlation']) > CORR_THRESHOLD]
    
    df_filtered = df_filtered.assign(abs_corr=np.abs(df_filtered['correlation']))
    df_filtered = df_filtered.sort_values('abs_corr', ascending=False).drop(columns='abs_corr')
    
    if len(df_filtered) > TOP_N:
        df_filtered = df_filtered.head(TOP_N)
    
    df_filtered.to_csv(OUT_CSV, index=False)
    
    print(df_filtered.head(10))


if __name__ == "__main__":
    main()
