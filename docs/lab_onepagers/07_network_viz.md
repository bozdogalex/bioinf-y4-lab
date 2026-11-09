# Lab 7 — Vizualizarea și interpretarea aceleiași rețele

Încărcăm `adjacency.csv`, `modules.csv` și `parameters.json` din `data/work/student/lab06/`, fără reconstruirea adiacenței. Astfel, o schimbare de layout nu modifică rețeaua analizată.

```bash
python labs/07_network_viz/ex01_network_viz.py --handle student
```

## Ce urmărim

- Poziția în layout este o alegere de afișare. Culorile indică modulele; izolatele sunt păstrate și desenate gri.
- Gradul numără muchii, strength însumează ponderi, iar conectivitatea intramodulară însumează ponderile spre genele din același modul.
- Eigengene-ul rezumă modulul prin prima componentă principală; poate fi pus în relație cu caracteristici ale probelor.
- Stabilitatea la reeșantionare și verificarea în date independente completează interpretarea.
- Genele sintetice nu au adnotări biologice. Enrichment GO/KEGG se demonstrează separat pe gene reale, cu background și corecție pentru testări multiple.
- Diseasome folosește alte noduri și alte relații: de exemplu, boli legate prin gene asociate comune. Nu este o conversie automată a corelațiilor în comorbidități.

Scriptul exportă `network.png`, `hubs.csv` și `visualization_parameters.json` în `data/work/student/lab07/`, ignorat de Git. Aceste rezultate rămân local.

## Exercițiul obligatoriu

Rulați același graf cu `--seed 42` și `--seed 7`, folosind `--output-dir` separat pentru fiecare rezultat. Comparați figurile și tabelele de hub-uri. În programul propriu `submissions/<handle>/analyze_hub.py`, calculați scorurile unei gene pe graful complet, apoi extrageți gena și vecinii direcți pentru afișare. Explicați ce susține simularea și ce date ar fi necesare pentru interpretarea biologică.

## Extensie opțională

```bash
python labs/07_network_viz/ex02_module_profiles.py --handle student
```

Calculează profilurile PC1 din expresia standardizată, proporția variației explicate și corelația descriptivă cu `control`/`tratament`. Folosește expresia și metadatele salvate în Lab 6, aliniate după probe. Rezultatele rămân în `data/work/student/lab07/profiles/`.

În PR-ul privat includeți programul propriu și `labs/07_network_viz/submissions/student/notes.md`, cu comenzile de reproducere, parametrii, hash-ul intrării și un rezumat al rezultatelor. Păstrați scripturile didactice nemodificate. Trimiteți linkul PR și commit SHA în LMS. Opțional, includeți cel mult 10 linii de metadate; nu adăugați matrici, liste complete sau PNG-uri generate. [Politica de date](../GDPR_and_DataPolicy.md).

- [Instrucțiuni complete](../../labs/07_network_viz/README.md)
- [Enunț și punctaj](../../labs/07_network_viz/assignment.md)
- [Prezentare în română](../presentations/07_network_viz.pdf)
- [NetworkX Drawing](https://networkx.org/documentation/stable/reference/drawing.html)
- [Goh et al., The human disease network (2007)](https://doi.org/10.1073/pnas.0701361104)

Material actualizat cu asistență AI, sub îndrumarea cadrului didactic.
