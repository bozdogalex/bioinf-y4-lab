# Raport: Analiza Rețelei de Co-expresie Genică TP53
**Student:** AlexTGoCreative  
**Data:** 6 Ianuarie 2026  
**Laborator:** 07 - Network Visualization & Diseasome

---

## 1. Descrierea Rețelei și Modulelor Detectate

### 1.1 Caracteristici generale
Am construit o rețea de co-expresie genică folosind 30 de gene relevante pentru cancerul asociat cu TP53, analizate pe 12 probe TCGA. Rețeaua finală conține:
- **23 noduri** (gene)
- **93 muchii** (conexiuni)
- **3 module distincte** (detectate cu algoritmul Louvain)

### 1.2 Modulele identificate

**Modul 0 (8 gene) - P53 Pathway & Apoptosis:**
- Gene: TP53, MDM2, CDKN1A, BAX, PUMA, NOXA, GADD45A, PTEN
- Funcție biologică: Răspuns la stres celular, apoptoză, oprirea ciclului celular
- Pattern de expresie: Expresie înaltă în probe 1-4, scăzută în restul

**Modul 1 (7 gene) - Cell Cycle Regulation:**
- Gene: CCND1, CCNE1, CDK2, CDK4, E2F1, RB1, CCNA2
- Funcție biologică: Reglarea progresiei prin ciclul celular, tranziții G1/S și G2/M
- Pattern de expresie: Expresie înaltă în probe 5-8, scăzută în restul

**Modul 2 (8 gene) - DNA Repair & Genome Stability:**
- Gene: BRCA1, BRCA2, ATM, ATR, CHEK1, CHEK2, RAD51, XRCC1
- Funcție biologică: Repararea daunelor ADN, menținerea stabilității genomului
- Pattern de expresie: Expresie înaltă în probe 9-12, scăzută în restul

### 1.3 Hub Genes Identificate

Top 10 hub genes bazat pe grad și betweenness centrality:

| Gene | Degree | Betweenness | Modul | Rol Biologic |
|------|--------|-------------|-------|--------------|
| CCND1 | 14 | 0.087 | 1 | Reglator ciclu celular G1/S |
| CCNE1 | 14 | 0.087 | 1 | Reglator ciclu celular G1/S |
| TP53 | 9 | 0.000 | 0 | Supresor tumoral central |
| MDM2 | 9 | 0.000 | 0 | Inhibitor TP53, ubiquitin ligase |
| BAX | 9 | 0.000 | 0 | Pro-apoptotic, activat de TP53 |
| CDKN1A | 9 | 0.000 | 0 | Inhibitor CDK (p21), arrest ciclu |

**Observație:** CCND1 și CCNE1 sunt hub genes dominante (betweenness > 0), acționând ca punți între module, ceea ce sugerează rol în coordonarea ciclului celular cu răspunsul la stres.

---

## 2. Analiza Funcțională - Modul 0 (P53 Pathway)

Am selectat **Modulul 0** pentru analiză de îmbogățire deoarece conține gena centrală TP53 și reprezintă calea canonică de răspuns la stres celular.

### 2.1 Gene analizate (Modul 0)
TP53, MDM2, CDKN1A, BAX, PUMA, NOXA, GADD45A, PTEN

### 2.2 Îmbogățire funcțională GO/KEGG (rezultate g:Profiler)

Analiza de îmbogățire a fost realizată pe 8 gene din Modulul 0 (TP53, MDM2, CDKN1A, BAX, PUMA, NOXA, GADD45A, PTEN) folosind g:Profiler. Rezultatele confirmă implicarea directă în calea p53 și apoptoză.

**Top 10 GO Molecular Function (GO:MF):**

| Rank | GO Term | Description | p-value |
|------|---------|-------------|---------|
| 1 | GO:0019904 | protein domain specific binding | 2.603×10⁻³ |
| 2 | GO:0031625 | ubiquitin protein ligase binding | 1.113×10⁻² |
| 3 | GO:0042802 | identical protein binding | 1.332×10⁻² |
| 4 | GO:0046982 | protein heterodimerization activity | 1.849×10⁻² |
| 5 | GO:1990841 | promoter-specific chromatin binding | 2.783×10⁻² |
| 6 | GO:0002039 | **p53 binding** ✓ | 2.865×10⁻² |
| 7 | GO:0044877 | protein-containing complex binding | 4.012×10⁻² |

**Interpretare MF:** Genele interacționează direct cu p53 (GO:0002039) și sunt implicate în ubiquitinare (MDM2 → degradare TP53) și heterodimerize (BAX/BCL2).

**Top 15 GO Biological Process (GO:BP) - Cele mai semnificative:**

| Rank | GO Term | Description | p-value |
|------|---------|-------------|---------|
| 1 | GO:0071214 | cellular response to abiotic stimulus | **1.501×10⁻⁸** |
| 2 | GO:0051248 | negative regulation of protein metabolic process | **9.737×10⁻⁷** |
| 3 | GO:0072332 | **intrinsic apoptotic signaling pathway by p53 class mediator** ✓ | **3.005×10⁻⁶** |
| 4 | GO:0051726 | **regulation of cell cycle** ✓ | **1.682×10⁻⁵** |
| 5 | GO:0072359 | circulatory system development | 2.300×10⁻⁵ |
| 6 | GO:0042770 | **signal transduction in response to DNA damage** ✓ | **5.853×10⁻⁵** |
| 7 | GO:0072717 | cellular response to actinomycin D | 3.635×10⁻⁴ |
| 8 | GO:1904705 | regulation of vascular associated smooth muscle cell proliferation | 4.123×10⁻⁴ |
| 9 | GO:2000379 | positive regulation of reactive oxygen species metabolic process | 6.265×10⁻⁴ |
| 10 | GO:0000079 | **regulation of cyclin-dependent protein serine/threonine kinase activity** ✓ | **9.394×10⁻⁴** |
| 11 | GO:0048145 | regulation of fibroblast proliferation | 1.342×10⁻³ |
| 12 | GO:0090398 | **cellular senescence** ✓ | **1.901×10⁻³** |
| 13 | GO:0050793 | regulation of developmental process | 2.261×10⁻³ |
| 14 | GO:0031667 | response to nutrient levels | 4.408×10⁻³ |
| 15 | GO:0070230 | positive regulation of lymphocyte apoptotic process | 9.252×10⁻³ |

**Procese cheie identificate:**
- ✅ **Apoptoză mediată de p53** (GO:0072332, p=3.0×10⁻⁶) - FOARTE SEMNIFICATIV
- ✅ **Reglarea ciclului celular** (GO:0051726, p=1.7×10⁻⁵) - confirmat
- ✅ **Răspuns la daune ADN** (GO:0042770, p=5.9×10⁻⁵) - rol central p53
- ✅ **Senescență celulară** (GO:0090398, p=1.9×10⁻³) - alternativă la apoptoză
- ✅ **Reglare CDK** (GO:0000079, p=9.4×10⁻⁴) - prin CDKN1A (p21)

**GO Cellular Component (GO:CC):**

| Rank | GO Term | Description | p-value |
|------|---------|-------------|---------|
| 1 | GO:0016604 | nuclear body | 1.134×10⁻⁴ |
| 2 | GO:0017053 | transcription repressor complex | 3.444×10⁻² |
| 3 | GO:0097144 | **BAX complex** ✓ | **4.998×10⁻²** |

**Human Phenotype (HP):**

| Rank | HP Term | Description | p-value |
|------|---------|-------------|---------|
| 1 | HP:0030448 | **Soft tissue sarcoma** | **6.640×10⁻...** |

**Interpretare HP:** Asociere cu sarcom (cancer țesut moale), consistent cu mutații TP53 în Li-Fraumeni syndrome.

### 2.3 Interpretare biologică

**TP53** este centrul acestui modul, funcționând ca "gardianul genomului":
1. **Răspuns la stres:** Activare în caz de daune ADN, stres oncogenic
2. **Activare transcrițională:** 
   - CDKN1A (p21) → arrest ciclu celular G1/S
   - BAX, PUMA, NOXA → inițierea apoptozei
   - GADD45A → reparare ADN
3. **Feedback negativ:** MDM2 ubiquitinează TP53 → degradare (menține homeostazia)
4. **PTEN:** Inhibă PI3K/AKT → sensibilizează la apoptoză

**Semnificație clinică:** Mutații în TP53 (>50% cancere) → pierderea checkpoint-urilor → proliferare necontrolată.

---

## 3. Integrare în Contextul Diseasome

### 3.1 Conceptul de Diseasome
Diseasome-ul reprezintă rețeaua de conexiuni între boli și genele asociate. Genele din același modul tind să fie implicate în boli similare sau inter-relate.

### 3.2 Modulul 0 în Diseasome

**Conexiuni centrale TP53:**
```
TP53 Mutations → Multiple Cancers:
├─ Breast Cancer (BRCA connection via Modul 2)
├─ Colorectal Cancer
├─ Lung Cancer
├─ Ovarian Cancer
├─ Li-Fraumeni Syndrome (germline mutations)
└─ Increased risk pentru ~15+ cancer types
```

**Cross-talk între module și boli:**

1. **Modul 0 ↔ Modul 1:**
   - TP53 → CDKN1A (p21) → inhibă CDK2/CCND1
   - Mutații TP53 → pierderea inhibiției → hiperproliferare
   - **Conexiune în diseasome:** Cancer cu ciclu celular deregulat

2. **Modul 0 ↔ Modul 2:**
   - TP53 activează ATM/ATR în răspuns la daune ADN
   - BRCA1/2 (Modul 2) colaborează cu TP53 în reparare
   - **Conexiune în diseasome:** 
     - Pacienți BRCA1/2+ și TP53 mutant → prognostic foarte slab
     - Cancer mamar triplu negativ (basal-like)

3. **Hub genes inter-module:**
   - CCND1/CCNE1 (betweenness ridicat) → legături între apoptoză și proliferare
   - **Implicații:** Target terapeutic pentru CDK4/6 inhibitori (palbociclib)

### 3.3 Perspectiva Diseasome - Network Medicine

**Principii aplicabile:**

1. **Disease modules:** Genele din Modul 0 formează un "disease module" pentru cancere TP53-mutante

2. **Disease-disease relationships:**
   - Pacienți cu Li-Fraumeni (TP53 germline) → risc pentru multiple cancere
   - Overlap genic → overlap fenotipic

3. **Comorbidități:**
   - Sindroame cancer + imunodepresie (ex: p53 și reglare imună)
   - Aging accelerat (senescență celulară controlată de p53/p21)

4. **Drug repurposing opportunities:**
   - MDM2 inhibitori (Nutlin-3) → reactivare p53 în tumori wild-type
   - PARP inhibitori (Modul 2: BRCA) + activatori p53
   - Combinații sinergice bazate pe topologia rețelei

### 3.4 Predicții clinice bazate pe modul

- **Biomarkeri prognostici:** Expresia coordonată în Modul 0 → răspuns la chimioterapie
- **Stratificare pacienți:** 
  - TP53 wild-type + expresie înaltă Modul 0 → candidați pentru MDM2 inhibitori
  - TP53 mutant → terapii alternative (restaurare funcție p53)
- **Vulnerabilități sintetice:** Pierderea p53 + inhibiție CHK1/CHK2 → synthetic lethality

---

## 4. Concluzii

1. **Modularitatea clară:** 3 module distincte reflectă funcții biologice complementare esențiale pentru homeostazia celulară

2. **TP53 ca hub în diseasome:** Modulul 0 este central în patogeneza cancerului, conectând multiple căi de semnalizare

3. **Cross-talk funcțional:** Hub genes inter-module (CCND1, CCNE1) sugerează coordonare fină între apoptoză și proliferare

4. **Aplicații clinice:** Înțelegerea rețelei permite:
   - Identificare biomarkeri
   - Design terapii combinate
   - Predicție răspuns la tratament
   - Drug repurposing bazat pe topologie

5. **Limitări:** 
   - Dataset sintetic (30 gene, 12 probe) → validare necesară pe date reale TCGA
   - Corelație ≠ cauzalitate → necesită validare experimentală
   - Lipsa datelor clinice asociate probelor

---

## Referințe

1. Barabási, A.L., et al. (2011). "Network medicine: a network-based approach to human disease." *Nature Reviews Genetics*, 12(1), 56-68.

2. Vogelstein, B., et al. (2013). "Cancer genome landscapes." *Science*, 339(6127), 1546-1558.

3. Menche, J., et al. (2015). "Uncovering disease-disease relationships through the incomplete interactome." *Science*, 347(6224), 1257601.

4. Lane, D.P. (1992). "Cancer. p53, guardian of the genome." *Nature*, 358(6381), 15-16.

5. Zhang, B., & Horvath, S. (2005). "A general framework for weighted gene co-expression network analysis." *Statistical Applications in Genetics and Molecular Biology*, 4(1).

---
