"""
Demo — Calculul corelației și pragului pe date toy

Acest demo arată cum se calculează o matrice de corelație între gene
și cum se aplică un prag pentru a obține o matrice de adiacență simplă.
"""

import pandas as pd
import numpy as np

# Toy expression matrix (3 gene x 5 samples)
data = {
    "Sample1": [5, 3, 8],
    "Sample2": [4, 3, 9],
    "Sample3": [6, 2, 7],
    "Sample4": [5, 4, 10],
    "Sample5": [4, 3, 8],
}
df = pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC"])
print("Matrice expresie (toy):")
print(df)

# Corelație Spearman între gene
corr = df.T.corr(method="spearman")
print("\nMatrice corelație (Spearman):")
print(corr)

# Prag pe corelație (ex: 0.7)
threshold = 0.7
adjacency = (corr >= threshold).astype(int)
np.fill_diagonal(adjacency.values, 0)  # fără self-loops

print(f"\nMatrice adiacență cu prag {threshold}:")
print(adjacency)
