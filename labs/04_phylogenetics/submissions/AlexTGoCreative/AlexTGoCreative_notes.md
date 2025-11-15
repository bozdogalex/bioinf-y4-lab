# Săptămâna 4 — Note și Reflecții asupra Filogeneticii

## Secvențe FASTA utilizate

**Fișier sursă:** `data/sample/tp53_protein_multi.fasta`

**Descriere:**
- Fișier multi-FASTA conținând secvențe proteice TP53 din multiple organisme
- Include variante de TP53 din diferite specii (Homo sapiens, Mus musculus, Danio rerio)
- Secvențele au fost aliniate cu MUSCLE înainte de construirea arborelui
- Salvat local în: `data/work/AlexTGoCreative/lab04/your_sequences.fasta`

**Accession-uri incluse:**
- TP53 uman (Homo sapiens) - P04637
- TP53 șoarece (Mus musculus) 
- TP53 pește zebră (Danio rerio)

---

## Reflecție: Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?

### Context
O **matrice de distanțe** conține doar valori numerice care descriu cât de diferite sunt secvențele între ele, două câte două. Un **arbore filogenetic**, în schimb, organizează aceste relații într-o reprezentare ierarhică și evolutivă.

### Avantajele arborelui filogenetic

#### 1. **Relații evolutive și ancestralitate**
- **Matricea de distanțe:** Oferă doar "cât de departe" sunt două secvențe
- **Arborele filogenetic:** Arată **drumul evolutiv** dintre secvențe și indică potențiali **ancești comuni**
- Exemplu: Putem vedea că TP53 uman și șoarece sunt mai apropiate între ele decât față de TP53 peștele zebră, sugerând un ancestor comun mai recent între mamifere

#### 2. **Ierarhia relațiilor**
- **Matricea de distanțe:** Date "plate" - toate comparațiile au aceeași importanță
- **Arborele filogenetic:** Structură ierarhică - grupează secvențele în **clade** (grupuri monofiletice)
- Permite identificarea rapid a grupurilor de secvențe înrudite

#### 3. **Vizualizare intuitivă**
- **Matricea de distanțe:** Tabel de numere dificil de interpretat pentru >10 secvențe
- **Arborele filogenetic:** Reprezentare grafică care arată imediat:
  - Ce secvențe sunt cel mai apropiate
  - Ordinea divergenței evolutive
  - Structura generală a relațiilor

#### 4. **Lungimea ramurilor = distanța evolutivă**
- Ramurile mai lungi → mai multe mutații acumulate
- Ramurile scurte → conservare evolutivă
- Permite identificarea **regiunilor sub selecție** (ramuri scurte) sau **regiunilor neutre** (ramuri lungi)

#### 5. **Inferență a caracterelor ancestrale**
- Arborele permite **reconstrucția stărilor ancestrale** la nodurile interne
- Putem estima cum arăta secvența TP53 la ancestrul comun al tuturor vertebratelor
- Matricea de distanțe nu oferă această informație

#### 6. **Testarea ipotezelor evolutive**
- Arborele poate fi folosit pentru a testa dacă anumite **evenimente evolutive** (duplicări genice, speciații) sunt coerente cu datele
- Permite compararea cu arbori bazați pe alte gene sau pe cronologia geologică
- Suportă analize statistice (bootstrap, likelihood ratio tests)

### Exemplu concret din lab

În exercițiul nostru cu TP53:
- **Matricea de distanțe** ne spune că distanța dintre TP53 uman și șoarece este, de exemplu, 0.15
- **Arborele NJ** ne arată că:
  - Omul și șoarecele formează un clade (mamifere)
  - Peștele zebră este outgroup (vertebrate non-mamifere)
  - Lungimea ramurilor indică rata de evoluție moleculară
  - Topologia arborelui reflectă istoria evolutivă

### Concluzie

Arborele filogenetic **nu înlocuiește** matricea de distanțe — este construit **pe baza** ei. Dar transformă date numerice brute în **ipoteze testabile despre istoria evolutivă**. 

Pentru cercetători, arborele este esențial pentru:
- Înțelegerea funcțiilor conservate vs. inovatoare ale genelor
- Identificarea ortologilor și paralogilor
- Predicția funcției proteinelor necunoscute bazată pe omologie
- Contextul evolutiv al variantelor genetice în medicină

**În rezumat:** Matricea de distanțe = "cât de diferite sunt secvențele". Arborele filogenetic = "cum au evoluat aceste diferențe și ce ne spune asta despre biologie".

---

## Metoda utilizată

- **Algoritm:** Neighbor-Joining (NJ)
- **Metric de distanță:** Identity (proporția de poziții identice în aliniere)
- **Software de aliniere:** MUSCLE (opțional, dacă este disponibil)
- **Format output:** Newick (.nwk)

---
