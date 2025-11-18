import pandas as pd
from pathlib import Path

in_path = Path("data/work/eeeevn0/lab06/GSE115469_Data.csv")
out_path = Path("data/work/eeeevn0/lab06/expression_matrix1.csv")

print("[INFO] Reading full matrix...")
df = pd.read_csv(in_path, index_col=0)

print("[INFO] Full shape:", df.shape)

num_genes = 2000     
df_subset = df.sample(n=num_genes, random_state=42)
df_subset.to_csv(out_path)

print("[OK] Subset salvat :", out_path)
print("[INFO] Subset shape:", df_subset.shape)
