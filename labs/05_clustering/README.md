# Săptămâna 5 — Clustering în Bioinformatică

## Scopuri
- Înțelegerea metodelor de bază pentru clustering: **Hierarchical, K-means, DBSCAN**.  
- Aplicarea clustering-ului pe date biologice reale.  
- Vizualizarea și interpretarea clusterelor.  
- Conectarea clustering-ului la filogenetică și la aplicații avansate (drug repurposing, gene co-expression networks).  

---

## Context
După ce în săptămâna 4 am construit arbori filogenetici, acum aplicăm metode de **clustering** pentru a descoperi grupări ascunse în date biologice.  
În prezentare vom discuta atât metode clasice (Hierarchical, K-means, DBSCAN), cât și metode avansate (dimensionality reduction, cluster validity indices, probabilistic și fuzzy clustering).  
În laborator vom implementa algoritmii de bază pe un dataset de cancer mamar.  

---

## Hands-on
**Rulați**  
- `demo01_k_means.py` — k means cu vizualizare PCA
**Rulați și completați**  
- `ex01_clustering.py` — clustering pe date de expresie genică din cancer mamar (toy dataset).  
  - Standardizați datele.  
  - Aplicați **Hierarchical Clustering** și vizualizați dendrograma.  
  - Aplicați **K-means (K=2)** și vizualizați rezultatele cu PCA.  
  - Aplicați **DBSCAN** și comparați clusterele.  
  - Salvați fișierele generate (CSV cu etichete, imagini cu ploturi).  

---

## Livrabile
În PR trebuie să apară:
1. Fișierul `labs/05_clustering/submissions/<github_handle>_notes.md` cu:  
   - ce metodă ați considerat cea mai potrivită pentru datele analizate,  
   - o scurtă reflecție: **Cum se compară clustering-ul cu arborii filogenetici în descoperirea relațiilor biologice?**  
2. Scriptul completat `ex06_clustering.py`.  
3. Fișierele generate:  
   ```bash
   labs/05_clustering/submissions/<handle>/clusters_<handle>.csv
   labs/05_clustering/submissions/<handle>/hierarchical_<handle>.png
   labs/05_clustering/submissions/<handle>/kmeans_<handle>.png
   labs/05_clustering/submissions/<handle>/dbscan_<handle>.png
   ```
4. Completarea checklist-ului din șablonul PR.

---

## Săptămâna următoare

- Gene Co-expression Networks.
- De la clustering clasic la module biologice și integrare multi-omics.
- [Vezi Săptămâna 6 — Gene Co-expression Networks](./06_wgcna/README.md))

---

## Competențe

- aplicarea Hierarchical, K-means și DBSCAN pe date biologice.
- Reducerea dimensionalității și vizualizarea clusterelor (PCA).
- Compararea și interpretarea rezultatelor între metode.
- Înțelegerea limitelor clustering-ului și legătura cu analiza filogenetică. 

---

## Resurse

-[Fișa laborator](../../docs/lab_onepagers/05_clustering.md)
-[Scikit-learn Clustering](https://scikit-learn.org/stable/modules/clustering.html)
-[Scipy Hierarchical](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html)
-[PCA în scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)
-[Drug repurposing cu clustering pe rețele de gene.](../../docs/papers/clustering_paper.pdf)