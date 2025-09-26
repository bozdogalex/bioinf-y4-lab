# Assignment Gene Co-Expression Networks, Visualization & Diseasome

## Obiective
- construirea rețelelor de co-expresie genică (GCEs),
- detectarea modulelor,
- vizualizarea și identificarea genelor hub,
- validarea biologică prin analiză de îmbogățire,

---

## Task-uri

### 1. Date și preprocesare (2p)
- Refolosiți dataset-ul TP53 asociat din **Lab 3**.  
- Preprocesați datele:
  - log2(x+1),
  - filtrați genele cu varianță scăzută.  

### 2. Rețea și module (3p)
- Construiți matricea de corelație (Pearson sau Spearman).  
- Aplicați un prag pentru a obține matricea de adiacență.  
- Construiți graful și detectați module cu algoritmul **Louvain**.  
- Exportați `modules_tp53_<handle>.csv` (gene → module).  

### 3. Vizualizare și hub genes (3p)
- Vizualizați rețeaua cu **NetworkX** sau un tool extern (Cytoscape/Gephi).  
- Colorați nodurile după modul.  
- Evidențiați genele hub 
- Exportați `network_tp53_<handle>.png`.  
- Exportați `hubs_tp53_<handle>.csv`.  

### 4. Interpretare biologică și diseasome (2p)
- Alegeți **un modul** și faceți analiză de îmbogățire (GO/KEGG) cu g:Profiler sau DAVID.  
- Scrieți un scurt raport (max 2 pagini PDF) care să includă:
  - descrierea modulului și genelor hub,  
  - rezultatele de îmbogățire,  
  - reflecție: *Cum ar putea acest modul să se integreze într-un context de tip diseasome?*  

### Bonus (+1p)
- Comparați vizualizările dintre NetworkX și Cytoscape/Gephi (screenshot + observații).  

---

## Livrabile
Încărcați într-un PR un fișier `.zip` cu:
- `modules_tp53_<handle>.csv`  
- `network_tp53_<handle>.png`  
- `hubs_tp53_<handle>.csv`  
- Codurile folosite (`.py` sau `.ipynb`)  
- `report_<handle>.pdf` (max 2 pagini, cu interpretarea și referința la diseasome)  

