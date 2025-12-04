# Săptămâna 5 — Note și Reflecții asupra Clustering-ului

## Date utilizate

**Dataset:** Date de expresie genică din cancer mamar (toy dataset)
- **Sursă:** Dataset inclus în laborator pentru demonstrarea metodelor de clustering
- **Caracteristici:** Profile de expresie genică pentru eșantioane de cancer mamar
- **Scop:** Identificarea subgrupurilor de pacienți bazate pe pattern-uri de expresie genică

---

## Metodă considerată cea mai potrivită

### Analiza comparativă a metodelor

Am aplicat trei metode principale de clustering pe datele de expresie genică:

#### 1. **Hierarchical Clustering**
**Avantaje:**
- Vizualizare intuitivă prin dendrogramă
- Nu necesită specificarea prealabilă a numărului de clustere
- Arată ierarhia relațiilor între eșantioane
- Potrivit pentru date biologice unde relațiile ierarhice sunt importante

**Dezavantaje:**
- Sensibil la outlieri
- Computațional intens pentru dataset-uri mari
- Deciziile timpurii (la bază) afectează întreaga structură

#### 2. **K-means (K=2)**
**Avantaje:**
- Rapid și eficient computațional
- Simplu de implementat și interpretat
- Funcționează bine când clusterele sunt sferice și de dimensiuni similare

**Dezavantaje:**
- Necesită specificarea prealabilă a K
- Sensibil la inițializare (rezultate pot varia)
- Presupune clustere de formă sferică și densitate similară
- Nu gestionează bine noise-ul

#### 3. **DBSCAN**
**Avantaje:**
- Nu necesită specificarea numărului de clustere
- Identifică outlieri (puncte de zgomot)
- Găsește clustere de forme arbitrare
- Robust la zgomot

**Dezavantaje:**
- Sensibil la parametrii eps și min_samples
- Dificil de aplicat pe date de dimensionalitate înaltă
- Nu funcționează bine când clusterele au densități diferite

### **Metoda recomandată: Hierarchical Clustering**

Pentru datele de expresie genică din cancer mamar, consider **Hierarchical Clustering** cea mai potrivită metodă din următoarele motive:

1. **Context biologic relevant:** Subtipurile de cancer mamar (luminal A, luminal B, HER2+, basal-like) au relații ierarhice bazate pe profile moleculare

2. **Vizualizare intuitivă:** Dendrograma oferă o vedere de ansamblu asupra relațiilor între eșantioane, permițând identificarea atât a grupurilor majore cât și a subgrupurilor

3. **Flexibilitate:** Nu trebuie să stabilim априори câte subtipuri căutăm - putem "tăia" dendrograma la diferite niveluri pentru a explora grupări la diverse rezoluții

4. **Validare clinică:** Rezultatele pot fi comparate cu clasificări clinice existente (PAM50, intrinsic subtypes)

**Notă:** În practică, o abordare **combinată** este optimă:
- Hierarchical pentru explorare inițială și vizualizare
- K-means pentru validare și clustering rapid pe subseturi
- DBSCAN pentru identificarea outlierilor (eșantioane atipice care pot reprezenta cazuri rare)

---

## Reflecție: Cum se compară clustering-ul cu arborii filogenetici în descoperirea relațiilor biologice?

### Context

Atât **clustering-ul** cât și **arborii filogenetici** sunt metode de grupare și vizualizare a relațiilor, dar au scopuri și interpretări fundamentale diferite.

### Diferențe fundamentale

#### 1. **Natura relațiilor**

**Arbori filogenetici:**
- Reprezintă **relații evolutive** și **descendență comună**
- Nodurile interne = **ancești comuni** reali (existați în trecut)
- Ramurile = **timp evolutiv** și **mutații acumulate**
- Interpretare: "Secvența A și B au un ancestor comun mai recent decât A și C"

**Clustering:**
- Reprezintă **similarități** în spațiul caracteristicilor
- Nodurile interne (în hierarchical) = **puncte de agregare matematică**, nu entități biologice
- Distanțele = **diferențe în profile de expresie/caracteristici**, nu timp
- Interpretare: "Eșantionul A și B au profile de expresie similare"

#### 2. **Presupuneri și constrângeri**

**Arbori filogenetici:**
- Presupun **un singur arbore adevărat** (istoria evolutivă reală)
- Constrângere: structură de arbore (fiecare nod are exact un părinte)
- Model evolutiv: mutații, selecție, drift genetic
- Validare prin consistența cu alte gene, fosile, geografie

**Clustering:**
- **Nu există o grupare "adevărată"** unică - depinde de scop și metodă
- Flexibilitate: clustere pot fi suprapuse (fuzzy) sau ierarhice
- Nu presupune un proces generativ specific
- Validare prin stabilitate, separabilitate, validare externă (labels clinice)

#### 3. **Tipuri de date și aplicații**

**Arbori filogenetici:**
- Input: **secvențe ADN/proteină** (caractere omogoase)
- Aplicații:
  - Clasificarea speciilor
  - Tracking-ul evoluției genelor (ortologi vs. paralogi)
  - Inferența funcției prin omologie
  - Epidemiologie moleculară (tracking virusuri)

**Clustering:**
- Input: **orice tip de date** (expresie genică, imagistică, clinice)
- Aplicații:
  - Identificarea subtipurilor de boli (cancer subtypes)
  - Descoperirea gene co-expression modules
  - Stratificarea pacienților pentru tratament personalizat
  - Drug repurposing (gruparea medicamentelor cu efecte similare)

### Similarități și convergențe

#### 1. **Hierarchical clustering ≈ phylogeny superficial**

- Ambele produc **dendrograme** (arbori ierarhici)
- Ambele folosesc **distanțe între entități**
- Ambele permit **vizualizare intuitivă** a relațiilor

**Dar:** Dendrograma de clustering nu are interpretare evolutivă!

#### 2. **Metode computaționale suprapuse**

- Unii algoritmi sunt comuni:
  - **UPGMA** (Unweighted Pair Group Method with Arithmetic Mean) - folosit în ambele
  - **Neighbor-Joining** - original pentru filogenie, dar e un algoritm de clustering
- Ambele folosesc **matrice de distanțe** ca input

#### 3. **Aplicații hibride**

În practică, metodele se completează:

**Exemplu 1: Evoluția cancerului**
- **Phylogeny** pentru tracking-ul evoluției clonale a tumorilor (cancer phylogenetics)
- **Clustering** pentru gruparea pacienților în subtipuri bazate pe profile moleculare

**Exemplu 2: Analiza metagenomică**
- **Phylogeny** pentru clasificarea taxonomică a microbilor
- **Clustering** pentru identificarea comunităților microbiene funcționale

**Exemplu 3: Gene families**
- **Phylogeny** pentru a reconstrui istoria duplicărilor genice
- **Clustering** pentru a grupa gene cu funcții similare (gene ontology)

### Când să folosim fiecare?

| Criteriu | Arbori filogenetici | Clustering |
|----------|---------------------|------------|
| **Tip de date** | Secvențe (ADN/proteină) | Orice (expresie, clinice, imagistică) |
| **Întrebare** | "Cum au evoluat aceste entități?" | "Ce grupări naturale există în date?" |
| **Interpretare** | Istorie evolutivă | Pattern-uri de similaritate |
| **Validare** | Consistență cu alte surse evolutive | Stabilitate, separabilitate |
| **Output** | Un arbore (istorie unică) | Multiple grupări posibile |

### Concluzie

**Arborii filogenetici** și **clustering-ul** sunt **complementare**, nu interschimbabile:

- **Phylogeny** răspunde la **"de unde venim?"** - reconstituie istoria evolutivă
- **Clustering** răspunde la **"ce grupuri există?"** - descoperă pattern-uri în date

În **bioinformatica modernă**, ambele sunt esențiale:
- Phylogeny pentru înțelegerea evoluției și funcției genelor
- Clustering pentru medicina personalizată și descoperirea de biomarkeri

**Metafora perfectă:**
- **Arborele filogenetic** = arborele genealogic al familiei (relații de rudenie)
- **Clustering-ul** = gruparea persoanelor după ocupație, stil de viață, preferințe (similarități, nu neapărat rudenie)

---

## Metode aplicate în exercițiu

### Preprocessing
- **Standardizare:** StandardScaler pentru normalizarea datelor de expresie

### Algoritmi
1. **Hierarchical Clustering:** Ward linkage cu dendrogramă
2. **K-means:** K=2, vizualizare cu PCA
3. **DBSCAN:** eps=0.5, min_samples=5

### Reducerea dimensionalității
- **PCA (Principal Component Analysis):** Pentru vizualizarea clusterelor în 2D

---

## Referințe

- Sørlie T, et al. (2001). "Gene expression patterns of breast carcinomas distinguish tumor subclasses with clinical implications". PNAS. 98(19):10869-74.
- Scikit-learn Clustering: https://scikit-learn.org/stable/modules/clustering.html
- Scipy Hierarchical Clustering: https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
- Huang S, et al. (2018). "Applications of clustering methods in biomedical image analysis". Genomics, Proteomics & Bioinformatics. 16(5):322-337.
