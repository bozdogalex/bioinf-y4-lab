# Săptămâna 8 — Machine Learning în Bioinformatică  
## Supervised, Unsupervised & Semi-Supervised Learning

## Scopuri
- Înțelegerea diferențelor dintre învățarea supravegheată, nesupravegheată și semi-supravegheată.  
- Aplicarea ML pe date de expresie genică.  
- Construirea unui pipeline complet de clasificare (Random Forest + Logistic Regression).  
- Explorarea clusterelor nesupravegheate (PCA + KMeans).  
- Simularea unui scenariu semi-supervised cu etichete lipsă.  
- Analiză critică a performanței modelelor ML în omică.

---

## Context
După ce în săptămânile 5–7 am explorat clustering-ul, rețelele de co-expresie și modulele, acum introducem fundația analizelor predictive:  
**Machine Learning aplicat pe date omice**.

ML este esențial pentru:
- clasificarea tipurilor de țesut,  
- subtipuri tumorale,  
- descoperirea biomarkerilor,  
- modele clinice predictive,  
- analize multicentru.

Laboratorul de azi include:
1. Supervised Learning (Random Forest, Logistic Regression)  
2. Unsupervised Learning (PCA + KMeans)  
3. Semi-Supervised (pseudo-labels, scenariu simplificat)

---

## Hands-on

### 🔹 `demo01_sup_unsup.py`
Demonstrație conceptuală:
- generare date sintetice 2D,  
- Random Forest vs KMeans,  
- vizualizare diferențe SL/UL.

### 🔹 `demo02_expr.py`
Pipeline complet pe date reale:
- încărcare subset expresie genică,  
- clasificare cu Random Forest,  
- vizualizare PCA,  
- clustering cu KMeans,  
- crosstab Label × Cluster.

### 🔹 `demo03_logreg.py`
Model liniar fundamental:
- Logistic Regression cu scaling,  
- comparație vizuală cu RF.

---

## Completați:
### `ex01_sup_pipeline.py`
Task-uri:
1. încărcați `expression_matrix_<handle>.csv`,  
2. separați X și y,  
3. encodați etichetele,  
4. împărțire stratificată train/test,  
5. antrenare Random Forest,  
6. generarea classification_report,  
7. salvarea matricei de confuzie și a feature importance.

### `ex02_logreg_vs_ref.py`
- antrenare Logistic Regression + Random Forest,  
- comparație rapoarte,  
- observații: modele liniare vs non-liniare.

---

## Livrabile
În PR trebuie să apară:

1. `classification_report_<handle>.txt`  
2. `confusion_rf_<handle>.png`  
3. `feature_importance_<handle>.csv`  
4. *(opțional)* `cluster_crosstab_<handle>.csv`  
5. checklist PR completat

---

## Săptămâna următoare
- Network-based Drug Repositioning  
- Similarity networks, bipartite drug–gene  
- Integrarea expresiei diferențiale și modulelor  
- Analiza topologică și prioritizarea țintelor

---

## Competențe
- Construirea unui pipeline ML end-to-end.  
- Evaluarea modelelor predictive.  
- Analiza comparației supervised vs unsupervised.  
- Înțelegerea aplicabilității semi-supervised în omică.  
- Interpretarea rezultatelor în context biologic.

---

## Resurse
- Scikit-learn Documentation  
- Breiman — Random Forests (2001)  
- Jolliffe — PCA (2002)  
- GTEx Portal  
- van Dam et al., 2018 — ML în omică  
- Ng & Jordan — Semi-Supervised Learning
