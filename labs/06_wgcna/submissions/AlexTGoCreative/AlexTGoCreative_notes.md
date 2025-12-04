# Lab 6 — Gene Co-Expression Networks — Notițe

**Autor:** AlexTGoCreative  
**Data:** 4 decembrie 2025

---

## Parametri utilizați

### Metrică de corelație
- **Metoda:** Spearman
- **Justificare:** Spearman este mai robustă la outlieri și captează relații monotone non-liniare între gene, fiind preferată pentru datele RNA-Seq care pot avea distribuții non-normale.

### Praguri aplicate
- **Prag varianță:** 0.1
  - Filtrează genele cu varianță scăzută care nu contribuie la co-expresie
- **Prag adiacență:** 0.6
  - Se folosesc valorile absolute ale corelațiilor (|cor| ≥ 0.6)
  - Selectează doar relațiile puternice dintre gene

### Configurație rețea
- **Tip rețea:** Neorientată (undirected)
- **Ponderare:** Binară (1 dacă cor ≥ threshold, 0 altfel)
- **Algoritm module:** Louvain community detection

---

## Reflecție: Cum diferă o rețea de co-expresie față de clustering-ul clasic?

### Clustering clasic
- **Perspectivă globală:** Grupează entități (gene sau sample) bazat pe similaritatea globală în spațiul multidimensional
- **Distanță uniformă:** Folosește o singură metrică de distanță pentru toate perechile
- **Hierarchie fixă:** Rezultatul este adesea un dendrogram cu grupuri distincte
- **Nu capturează relații:** Nu păstrează informații despre relațiile specifice dintre elemente individuale

### Rețele de co-expresie
- **Perspectivă locală:** Capturează **relații pereche (pairwise)** dintre gene
- **Structură de graf:** Fiecare genă este un nod, fiecare corelație puternică este o muchie
- **Informație bogată:** 
  - Identifică hub genes (gene centrale cu multe conexiuni)
  - Permite analiza gradului de conectivitate
  - Păstrează informația despre intensitatea relațiilor (dacă este ponderată)
- **Module dinamice:** Comunitățile detectate reflectă grupuri de gene cu pattern-uri coordonate de expresie
- **Interpretare biologică:** 
  - Genele din același modul pot fi co-reglate sau implicate în aceleași căi biologice
  - Hub genes pot fi regulatori cheie sau markeri de boală

### Avantaje rețele vs. clustering
1. **Granularitate:** Rețelele păstrează detalii despre fiecare relație, nu doar apartenența la un grup
2. **Flexibilitate:** Permit analiza la multiple niveluri (noduri individuale, module, întreaga rețea)
3. **Context biologic:** Reflectă mai bine interacțiunile moleculare reale și căile de semnalizare
4. **Vizualizare:** Oferă reprezentări vizuale intuitive ale relațiilor complexe

### Limitări
- Pragul de corelație poate fi arbitrar și influențează dramatic structura rețelei
- Computațional mai costisitor pentru seturi mari de date
- Interpretarea necesită cunoștințe de biologie sistemică
- Corelația nu implică cauzalitate

---

## Observații suplimentare

- Preprocesarea (log2 transformare + filtrare varianță) este crucială pentru a reduce zgomotul
- Algoritmul Louvain maximizează modularitatea și identifică comunități în mod eficient
- Pentru analize viitoare: identificarea hub genes, enrichment analysis pe module, validare experimentală
