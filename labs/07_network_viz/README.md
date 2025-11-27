# Săptămâna 7 — Vizualizarea și interpretarea rețelelor de co-expresie (GCEs) + Diseasome

## Scopuri
- Vizualizarea rețelelor de co-expresie genică și a modulelor identificate.  
- Identificarea genelor hub (cele mai conectate noduri din module).  
- Validarea modulelor prin analiză de îmbogățire funcțională (GO, KEGG).  
- Înțelegerea conceptului de **Diseasome** și legătura cu rețelele de co-expresie.  
- Interpretarea rezultatelor în context biomedical (cancer, comorbidități, drug repurposing).  

---

## Context
În săptămâna 6 am construit o rețea de co-expresie și am detectat module.  
Acum mergem mai departe: **vizualizare și interpretare**.  

Rețelele sunt greu de înțeles doar ca matrici — o vizualizare bună scoate la suprafață structura, modulele și genele hub.  
Mai mult, putem valida biologic modulele prin **analiză de îmbogățire funcțională** și putem conecta aceste rezultate la **diseasome** — harta bolilor umane legate prin gene și module comune.  

---

## Hands-on
**Rulați și completați**  
- `ex08_network_viz.py` — vizualizați rețeaua construită în Lab 6.  
  - Încărcați `modules_<handle>.csv` și matricea de adiacență.  
  - Colorați nodurile în funcție de modul.  
  - Evidențiați genele hub (cele cu cel mai mare grad).  
  - Exportați figura în `network_<handle>.png`.  

**La alegere (bonus)**  
- Încercați vizualizarea într-un tool extern (Cytoscape sau Gephi).  
- Comparați aspectul și ușurința interpretării față de NetworkX.  

---

## Livrabile
În PR trebuie să apară:
1. Fișierul `labs/07_networkviz/submission/<github_handle>_notes.md` cu:  
   - ce metodă de layout ați folosit (ex: spring, kamada-kawai),  
   - o scurtă reflecție: **Ce avantaje aduce vizualizarea față de analiza numerică din Lab 6?**  
2. Scriptul completat `ex01_network_viz.py`.  
3. Fișierul generat:  
   ```bash
   labs/07_networkviz/submissions/<github_handle>/network_<handle>.png
   ```
4. Completarea checklist-ului din șablonul PR.

---

## Săptămâna următoare
- Machine Learning în analiza datelor biomedicale.
- Clasificarea stărilor de boală și evaluarea performanței modelelor.
- [Vezi Săptămâna 8 — Machine Learning.](/labs/08_ML_flower)

---

## Competențe
- Vizualizarea rețelelor și interpretarea modulelor.
- Identificarea și analiza genelor hub.
- Aplicarea analizelor de îmbogățire funcțională pe module.
- Înțelegerea și explicarea conceptului de Diseasome.
- Conectarea rezultatelor la context biomedical (cancer, comorbidități, drug repurposing).

---

## Resurse

- [Fișa laborator](../../docs/lab_onepagers/07_network_viz.md)
- [NetworkX Drawing](https://networkx.org/documentation/stable/reference/drawing.html)
- [Cytoscape — vizualizare rețele biologice](https://cytoscape.org/)
- [Gephi — vizualizare rețele generale](https://gephi.org/)
- [g:Profiler — analiză de îmbogățire GO/KEGG](https://biit.cs.ut.ee/gprofiler/)
- [Barabási et al., Nature Genetics 2007 — The Human Diseasome](https://www.nature.com/articles/nrg2918)