Saved dendrogram to: labs/05_clustering/submissions/StanaAndrei/hierarchical_StanaAndrei.png
Saved KMeans plot to: labs/05_clustering/submissions/StanaAndrei/kmeans_StanaAndrei.png
/workspaces/bioinf-y4-lab/labs/05_clustering/submissions/StanaAndrei/ex01_clustering.py:91: UserWarning: You passed a edgecolor/edgecolors ('k') for an unfilled marker ('x').  Matplotlib is ignoring the edgecolor in favor of the facecolor.  This behavior may change in the future.
  plt.scatter(X_pca[mask, 0], X_pca[mask, 1], label=label_name, alpha=0.7, marker=marker, edgecolors="k", linewidths=0.3)
Saved DBSCAN plot to: labs/05_clustering/submissions/StanaAndrei/dbscan_StanaAndrei.png
Saved clusters CSV to: labs/05_clustering/submissions/StanaAndrei/clusters_StanaAndrei.csv

Summary:
Diagnosis distribution:
Diagnosis
0    357
1    212
Name: count, dtype: int64

KMeans cluster counts:
KMeans_Cluster
0    375
1    194
Name: count, dtype: int64

DBSCAN cluster counts (includes -1 for noise):
DBSCAN_Cluster
-1    550
 0     19
Name: count, dtype: int64




------------------------------------------
1.
 K-Means este cea mai potrivită metodă pentru acest set de date — oferă separarea cea mai coerentă, este rapidă și se aliniază bine cu natura datelor (două clase clare, bine delimitate).

2.
Clustering-ul (ex. K-Means, DBSCAN, hierarchical) grupează datele în funcție de similarități matematice între trăsături (feature vectors). Nu are nevoie de informații prealabile despre origini sau succesiuni evolutive; pur și simplu identifică grupuri similare în spațiul datelor.

Arborii filogenetici, în schimb, sunt o formă specializată de clustering ierarhic, dar cu semnificație biologică temporală: arată descendența comună și distanțele evolutive între specii sau gene, bazate pe mutații și moștenire, nu doar pe similarități numerice.