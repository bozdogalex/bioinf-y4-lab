import pandas as pd
import numpy as np
from pathlib import Path

HANDLE = "eeeevn0"
DATA_DIR = Path(f"data/work/{HANDLE}/lab10")
DATA_DIR.mkdir(parents=True, exist_ok=True)

np.random.seed(123)   # pentru reproductibilitate

n_samples = 60 # numar de probe
n_snps = 25 # numar de SNPs
n_genes = 40 # numar de gene

samples = [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]

# snp maxtrix
snp_names = [f"SNP_{i:02d}" for i in range(n_snps)]

# genotip realist: 0,1,2 cu probabilitate
snp_data = np.random.binomial(2, 0.35, size=(n_samples, n_snps))
df_snp = pd.DataFrame(snp_data, index=samples, columns=snp_names)

#matrice expresie genica
gene_names = [f"GENE_{i:02d}" for i in range(n_genes)]
expr_noise = np.random.normal(0, 1.5, size=(n_samples, n_genes))
df_exp = pd.DataFrame(expr_noise, index=samples, columns=gene_names)

# adaugare corelatii intre unele SNPs si gene
for i in range(5):
    snp = snp_names[i]

    gene_a = gene_names[2 * i]
    gene_b = gene_names[2 * i + 1]

    df_exp[gene_a] = (
        0.6 * df_snp[snp] +
        np.random.normal(0, 0.4, n_samples)
    )

    df_exp[gene_b] = (
        -0.5 * df_snp[snp] +
        np.random.normal(0, 0.4, n_samples)
    )

# salvare matrice
df_snp.to_csv(DATA_DIR / f"snp_matrix_{HANDLE}.csv")
df_exp.to_csv(DATA_DIR / f"expression_matrix_{HANDLE}.csv")

print("Synthetic multi-omics dataset generated successfully.")
