# Săptămâna 9 — Network-Based Drug Repurposing  
## Similarity Networks, Drug–Gene Interactions & Network Proximity

## Scopuri
- Înțelegerea conceptului de **drug repurposing** prin modelare de rețele.  
- Construirea rețelelor de tip **drug–gene** și **drug similarity networks**.  
- Aplicarea metodelor de **proximity**, **shortest paths** și **centrality** pentru prioritizare de medicamente.  
- Rularea unui demo complet pe un subset de date (drug–target interactions).  
- Dezvoltarea unei perspective critice asupra rețelelor folosite în descoperirea de medicamente.

---

## Context
După ce în săptămânile 6 și 7 am construit și interpretat rețele de gene (GCEs), iar în săptămâna 8 am introdus machine learning-ul, în această săptămână:

**aplicăm modelele de rețele pentru a sugera medicamente candidate în repurposing.**

Laboratorul urmărește logica folosită în cercetare și industrie:
- Open Targets  
- BenevolentAI (Baricitinib pentru COVID-19)  
- Atomwise  
- Pfizer / Novartis (polypharmacology)  

Ne concentrăm pe concepte fundamentale:
- rețele drug–gene (bipartite),  
- rețele drug–drug (similarity),  
- centralitate, proximitate, shortest paths.

---

## Hands-on

### 🔹 `demo01_construct_drug_gene_network.py`
Demonstrație simplificată:
- încărcarea unui subset de interacțiuni drug–gene (toy dataset),  
- construirea rețelei bipartite,  
- proiectarea pe layer-ul **drug similarity**,  
- vizualizare + colorare după numărul de gene țintă.

### 🔹 `demo02_network_proximity.py`
Demonstrație conceptuală:
- definirea unui set de *disease genes*,  
- calcularea distanței medii drug ↔ disease genes (network proximity),  
- prioritizarea medicamentelor.

### 🔹 `demo03_random_walk_basics.py`
Demonstrație educațională:
- implementare simplificată de **Random Walk with Restart (RWR)** pe o rețea mică,  
- ranking-ul nodurilor cele mai apropiate de un set țintă.

---

## Exerciții

### Exercițiul 1 — Drug Similarity Network
Completați:
- încărcarea unui fișier CSV cu drug–gene,  
- construirea rețelei bipartite,  
- proiectarea pe layer-ul de medicamente,  
- calcularea ponderilor (Jaccard similarity între seturi de gene),  
- vizualizare PCA / force layout.

### Exercițiul 2 — Disease Proximity
Pentru 3–5 gene asociate cu o boală:
- calculați distanța medie la fiecare medicament,  
- ordonați medicamentele,  
- discutați interpretabilitatea rezultatului.

---

## Livrabile
În PR trebuie să apară:

1. `labs/09_repurposing/submissions/<handle>/REPURPOSING_<handle>.csv`  
2. `network_drug_gene_<handle>.png`  
3. `drug_priority_<handle>.csv` (ranking după network proximity)  
4. `report_repurposing_<handle>.pdf` (max 3 pagini)  
5. checklist PR completat

---

## Săptămâna următoare
- Multi-layer networks  
- Integrarea omics (genomic, transcriptomic, proteomic) în rețele  
- Graph embeddings pentru predictivitate reală în drug discovery  

---

## Competențe
- Construirea rețelelor bipartite drug–gene.  
- Definirea metodelor de similaritate și distanță.  
- Calculul centralității și al proximității în rețele complexe.  
- Aplicarea metodelor moderne folosite în repurposing.  
- Interpretarea rezultatelor în context biologic și translational.

---

## Resurse
- Open Targets Platform  
- Gysi et al., *Network medicine framework for drug repurposing*  
- Cheng et al., *Drug–target interaction prediction*  
- Yıldırım et al., *Drug–gene network organization*  
- scikit-learn, networkx documentație  
