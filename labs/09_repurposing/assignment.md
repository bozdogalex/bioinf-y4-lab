# Lab 9 — Drug Repurposing Using Network-Based Approaches

## Obiective
- Construirea unei rețele bipartite drug–gene (subset din Open Targets / DrugBank-like).  
- Construirea unui layer de similaritate între medicamente.  
- Aplicarea metodelor de proximitate pentru ranking-ul medicamentelor.  
- Vizualizare și interpretare biologică.

---

## Task 1 — Datele și rețeaua bipartită (2p)
- Încărcați dataset-ul `drug_gene_<handle>.csv` (oferit).  
- Construți un graf bipartit:  
  - layer 1: drug nodes  
  - layer 2: gene nodes  
- Calculați pentru fiecare drug numărul de gene țintă.  
- Exportați `drug_summary_<handle>.csv`

---

## Task 2 — Proiecția pe layer-ul de medicamente (3p)
- Pentru fiecare pereche de medicamente, calculați similaritatea dintre seturile de gene.  
  Recomandat:
  - **Jaccard similarity**  
  - alternativ: Overlap Coefficient, Cosine on binaries
- Construiți rețeaua **drug similarity** (weighted graph).  
- Exportați: `drug_similarity_<handle>.csv`

---

## Task 3 — Disease Proximity (3p)
Folosind lista `disease_genes_<handle>.txt`:
- calculați distanța medie a fiecărui medicament față de genele bolii,  
- distanțele se calculează în rețeaua **drug–gene** (nu în similarity),  
- ordonați medicamentele după proximitate,  
- exportați: `drug_priority_<handle>.csv`

---

## Task 4 — Vizualizare (1p)
Generează:
- figura `network_drug_gene_<handle>.png`  
  - medicamente: albastru  
  - gene: roșu  
  - mărimea nodului = număr gene target  

---

## Task 5 — Raport (1p)
Raport PDF (max 3 pagini) cu:
- metodologia,  
- cele mai apropiate medicamente,  
- interpretare biologică,  
- limitări (bias adnotări, incompletitudine, topologie parțială).

---

## Bonus (+1p)
Implementați o versiune simplificată de **Random Walk with Restart (RWR)** pe rețeaua drug–gene și comparați ranking-ul cu proximitatea.

---

## Livrabile
- `drug_summary_<handle>.csv`  
- `drug_similarity_<handle>.csv`  
- `drug_priority_<handle>.csv`  
- `network_drug_gene_<handle>.png`  
- `report_repurposing_<handle>.pdf`

