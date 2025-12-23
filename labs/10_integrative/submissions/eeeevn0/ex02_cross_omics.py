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
from scipy.stats import pearsonr

HANDLE = "eeeevn0"
JOINT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/multiomics_concat_{HANDLE}.csv")
OUT_CSV = Path(f"labs/10_integrative/submissions/{HANDLE}/snp_gene_pairs_{HANDLE}.csv")

THRESHOLD = 0.5 # prag pentru |r|

#incarcare date
df = pd.read_csv(JOINT_CSV, index_col=0)

#separare SNPs vs genes
snp_cols = [c for c in df.columns if "SNP" in c]
gene_cols = [c for c in df.columns if c not in snp_cols]

print(f"Detected {len(snp_cols)} SNPs and {len(gene_cols)} genes")

#calculare corelatii
results = []

for snp in snp_cols:
    for gene in gene_cols:
        r, p = pearsonr(df[snp], df[gene])

        if abs(r) > THRESHOLD:
            results.append({
                "SNP": snp,
                "Gene": gene,
                "Pearson_r": r,
                "p_value": p
            })
#export rezultate
df_out = pd.DataFrame(results)

if not df_out.empty:
    df_out = df_out.sort_values("Pearson_r", key=abs, ascending=False)

df_out.to_csv(OUT_CSV, index=False)
print(f"Exercise completed. Results saved to {OUT_CSV}")
