"""
Demo 2 — PCA pe date integrate
"""

import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

df = pd.read_csv("multiomics_concat_demo.csv", index_col=0)

pca = PCA(n_components=2)
proj = pca.fit_transform(df.T)

plt.scatter(proj[:, 0], proj[:, 1], c="blue")
plt.title("PCA on Integrated Multi-Omics")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.show()
