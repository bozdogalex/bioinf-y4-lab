# Lab 8 — Machine Learning aplicat pe date omice (Supervised, Unsupervised, Semi-Supervised)

## Obiective
- Construirea unui pipeline complet de clasificare.  
- Analiza performanței unui model ML pe date biologice.  
- Explorarea clusterelor nesupravegheate.  
- Simularea unui scenariu semi-supervised.  
- Interpretarea biologică a rezultatelor.

---

## Task 1 — Pregătirea datelor (1p)
- Folosiți `expression_matrix_<handle>.csv`.  
- Separați **X (gene)** și **y (Label)**.  
- Encodați etichetele cu `LabelEncoder`.  
- Împărțiți în train/test stratificat.

---

## Task 2 — Supervised ML: Random Forest (3p)
Implementați un pipeline complet:
1. antrenarea modelului,  
2. generarea `classification_report`,  
3. calculul matricei de confuzie + salvare `.png`,  
4. extragerea `feature_importances_` + salvare `.csv`.

**Întrebări pentru raport**:
- Ce gene sunt cele mai importante?  
- Modelul greșește mai des între anumite clase? De ce?  

---

## Task 3 — Logistic Regression (opțional +1p)
- antrenați Logistic Regression cu scaling,  
- comparați performanța cu RF,  
- discutați diferențele dintre modele liniare vs non-liniare.

---

## Task 4 — Unsupervised ML: PCA + KMeans (2p)
1. rulați PCA pentru vizualizare,  
2. aplicați KMeans cu 2–4 clustere,  
3. salvați scatter PCA colorat după cluster,  
4. generați crosstab: Label × Cluster.

**Întrebări pentru raport**:
- Clusterizarea recuperează parțial structura de clase?  
- Există probe care sar în clusterul „greșit”? De ce?

---

## Task 5 — Semi-Supervised (mini-experiment) (2p)
Simulați lipsa etichetelor:
- marcați 30–50% dintre etichete drept „unknown”,  
- antrenați RF doar pe datele etichetate,  
- generați pseudo-etichete pentru cele lipsă,  
- reantrenați modelul pe setul complet.

**Întrebări**:
- Performanța crește sau scade?  
- De ce ajută pseudo-labeling în bioinformatică?

---

## Livrabile
Livrați în `labs/08_ml/submissions/<handle>/`:

1. `classification_report_<handle>.txt`  
2. `confusion_rf_<handle>.png`  
3. `feature_importance_<handle>.csv`  
4. `cluster_crosstab_<handle>.csv` *(opțional)*  
5. `sup_vs_unsup_scatter_<handle>.png`  
6. `report_lab08_<handle>.pdf` (max 3 pagini)  

---

## Structura raportului (max 3 pagini)
1. **Introducere** — tipul datelor, obiectivele analizei.  
2. **Supervised ML** — RF (și LogReg dacă ați făcut opționalul).  
3. **Unsupervised ML** — PCA + KMeans, crosstab.  
4. **Semi-Supervised** — pseudo-labels + discuție.  
5. **Interpretare biologică** — ce gene sunt informative și de ce.  
6. **Limitări** — bias, small sample size, variabilitate biologică.  

---

## Bonus (+1p)
Comparați vizual PCA înainte și după eliminarea a 10 gene cu varianță mică.  
Discută dacă vizibilitatea clusterelor devine mai bună.

