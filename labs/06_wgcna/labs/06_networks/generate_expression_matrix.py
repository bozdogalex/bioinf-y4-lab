import numpy as np
import pandas as pd
from pathlib import Path

# Config
output_dir = Path("data/work/MadalinaDenisa/lab06")
output_dir.mkdir(parents=True, exist_ok=True)
output_csv = output_dir / "expression_matrix.csv"

# Parametri matrice
n_genes = 200
n_samples = 10
np.random.seed(42)

# Simulare valori RNA-Seq (ex: numar de citiri)
data = np.random.poisson(lam=50, size=(n_genes, n_samples))  # medie 50 counts

# Creare index (gene) și coloane (samples)
genes = [f"Gene{i+1}" for i in range(n_genes)]
samples = [f"Sample{i+1}" for i in range(n_samples)]

# Construire DataFrame
df = pd.DataFrame(data, index=genes, columns=samples)

# Salvare CSV
df.to_csv(output_csv)

print(f"CSV generat cu {n_genes} gene × {n_samples} probe → {output_csv}")
