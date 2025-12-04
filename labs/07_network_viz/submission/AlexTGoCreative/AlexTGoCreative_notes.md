# Lab 7 — Vizualizarea rețelelor de co-expresie — Notițe

**Autor:** AlexTGoCreative  
**Data:** 4 decembrie 2025

---

## Metodă de layout utilizată

### Spring Layout (Fruchterman-Reingold)
- **Algoritm:** `nx.spring_layout()` cu seed=42 pentru reproducibilitate
- **Principiu:** Algoritmul Fruchterman-Reingold simulează un sistem fizic în care:
  - Muchiile acționează ca **arcuri** (springs) care atrag nodurile conectate
  - Nodurile se resping reciproc ca **particule încărcate**
  - Sistemul ajunge la un echilibru energetic minim
- **Avantaje:**
  - Nodurile puternic conectate (hub genes) tind să fie poziționate central
  - Modulele (comunități) se grupează natural în spațiu
  - Vizualizarea este estetică și intuitivă
- **Parametri:** Seed fix pentru consistență între rulări

### Alternative considerate
- **Kamada-Kawai:** Bun pentru rețele mai mici, optimizează distanțele geodezice
- **Circular layout:** Util pentru a evidenția simetria, dar pierde structura modulară
- **Force-Atlas (Gephi):** Similar cu spring, dar disponibil doar în tool-uri externe

---

## Reflecție: Ce avantaje aduce vizualizarea față de analiza numerică din Lab 6?

### Analiza numerică (Lab 6)
**Ce am obținut:**
- Matricea de corelație: valori numerice precise ale relațiilor
- Matricea de adiacență: reprezentare binară/ponderată a conexiunilor
- Module detectate: liste de gene grupate prin algoritmi (Louvain)
- Metrici: număr de noduri, muchii, modularitate

**Limitări:**
- Informația este abstractă și greu de interpretat rapid
- Nu putem vedea **structura globală** a rețelei
- Greu de identificat visual pattern-uri sau anomalii
- Nu evidențiază rolul central al unor gene specifice

### Vizualizarea (Lab 7)
**Avantaje clare:**

#### 1. **Intuiție instantanee**
- Vedem imediat câte module există și cât de separate sunt
- Identificăm rapid zonele dense vs. rare în rețea
- Observăm dacă rețeaua este conectată sau fragmentată

#### 2. **Identificarea hub genes**
- Nodurile mari/centrale sunt hub genes (grad mare)
- Vizual evident care gene sunt conectoare critice între module
- În context biologic: **hub genes = potențiali regulatori master sau markeri de boală**

#### 3. **Structura modulară**
- Culorile diferite pentru module fac grupurile evidente
- Putem vedea dacă modulele sunt:
  - Complet separate (module funcționale distincte)
  - Parțial suprapuse (funcții partajate)
  - Conectate prin hub genes (gene de interconectare biologică)

#### 4. **Pattern recognition**
- **Topologie scale-free:** Observăm dacă câteva gene dominante controlează rețeaua
- **Small-world properties:** Vedem dacă există "scurtături" între module distante
- **Outlieri:** Identificăm gene izolate sau conectate neașteptat

#### 5. **Comunicare științifică**
- O imagine valorează cât 1000 de rânduri de date
- Esențială în publicații și prezentări
- Permite non-experților să înțeleagă complexitatea biologică

#### 6. **Generare de ipoteze biologice**
- **Hub genes în cancer:** Pot fi oncogene sau tumor supresori
- **Module specifice:** Pot corespunde căilor de semnalizare sau proceselor celulare
- **Gene bridges:** Conectează module diferite → potențial pentru drug repurposing

### Exemplu concret (date toy)
Chiar și pentru 2 gene (GENE_A, GENE_B):
- **Numeric:** `cor=1.0, module=0` → abstract
- **Vizual:** Două noduri conectate, aceeași culoare → **imediat clar că sunt co-exprimate**

---

## Observații suplimentare

### Hub genes în context biologic
- În cancer mamar (GSE115469), hub genes ar putea fi:
  - **TP53, BRCA1, MYC:** Cunoscuți oncogene/tumor supresori
  - **Gene metabolice:** Dacă modulul este legat de metabolism
  - **Receptori hormonali:** ESR1, PGR în cancere hormono-dependente

### Limitări ale vizualizării
- **Scalabilitate:** Pentru >1000 noduri devine aglomerată (necesită filtrare)
- **2D vs 3D:** Pierdem informație spațială în proiecția 2D
- **Layout dependence:** Algoritmi diferiți → vizualizări diferite
- **Nu înlocuiește analiza statistică:** Vizualizarea completează, nu înlocuiește metricile numerice

### Integrare cu Diseasome
Vizualizarea rețelelor de co-expresie poate fi extinsă la:
- **Diseasome networks:** Noduri = boli, muchii = gene partajate
- **Drug-target networks:** Identificarea medicamentelor care țintesc hub genes
- **Multi-omics networks:** Integrarea cu proteomică, metabolomică

---

## Concluzie

Vizualizarea transformă **date abstracte** în **knowledge acționabil**. Combinația Lab 6 (rigoare numerică) + Lab 7 (intuiție vizuală) oferă o înțelegere completă a rețelelor de co-expresie genică.

**"Seeing is understanding"** — vizualizarea este esențială pentru descoperirea biologică.
