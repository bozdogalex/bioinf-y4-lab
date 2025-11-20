# Săptămâna 6 — Gene Co-Expression Networks (GCEs) — Building Modules

## Scopuri
- Înțelegerea conceptului de rețele de co-expresie genică (GCEs).  
- Preprocesarea datelor RNA-Seq (normalizare, log-transformare, filtrare).  
- Calcularea matricilor de corelație și construirea matricii de adiacență.  
- Detectarea modulelor folosind algoritmi de tip community detection (ex: Louvain).  
- Analiza critică a limitărilor și capcanelor interpretative.  

---

## Context
După ce în săptămâna 5 am studiat clustering-ul, acum trecem la **rețele de co-expresie genică**.  
Clustering-ul grupează probe sau gene global, în timp ce rețelele de co-expresie surprind **relațiile pereche** dintre gene.  
În această primă parte vom construi rețeaua și vom detecta modulele. Interpretarea și vizualizarea vor fi abordate în săptămâna 7.  

---

## Hands-on
**Rulați și completați**  
- `ex07_gce_networks.py` — construcția unei rețele de co-expresie folosind date RNA-Seq din cancer mamar (GSE115469).  
  - Încărcați și normalizați datele.  
  - Aplicați log-transformare și filtrați genele cu varianță scăzută.  
  - Calculați corelația (Pearson sau Spearman).  
  - Aplicați un prag pentru a construi matricea de adiacență.  
  - Construiți graful și detectați module cu algoritmul **Louvain**.  
  - Exportați un fișier `.csv` cu mapping gene → module.  

---

## Livrabile
În PR trebuie să apară:
1. Fișierul `labs/06_networks/submissions/<github_handle>_notes.md` cu:  
   - ce metrică de corelație și ce prag ați folosit,  
   - o scurtă reflecție: **Cum diferă o rețea de co-expresie față de clustering-ul clasic?**  
2. Scriptul completat `ex07_gce_networks.py`.  
3. Fișierul generat:  
   ```bash
   labs/06_networks/submissions/<handle>/modules_<handle>.csv
   ```
4. Completarea checklist-ului din șablonul PR.

---

## Săptămâna următoare
- Vizualizarea și interpretarea rețelelor.
- Identificarea hub genes și analiza funcțională (enrichment).
- Introducerea conceptului de Diseasome.
- [Vezi Săptămâna 7 — Visualization & Diseasome](./07_network_viz/README.md)

---

## Competențe

- Preprocesarea datelor RNA-Seq pentru analize de rețea.
- Construirea matricilor de corelație și adiacență.
- Detectarea modulelor cu algoritmi de tip community detection.
- Evaluarea critică a rețelelor.

---

## Resurse
- [Fișa laborator](../../docs/lab_onepagers/06_WGCNA.md)
- [NetworkX Documentation](https://networkx.org/documentation/stable/)
- [GEO Accession GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469)
- [van Dam et al., Brief Bioinform 2018](https://doi.org/10.1093/bib/bbw139)
- [Langfelder & Horvath, BMC Bioinformatics 2008](https://doi.org/10.1186/1471-2105-9-559)
