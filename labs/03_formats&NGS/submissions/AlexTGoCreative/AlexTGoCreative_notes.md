# Săptămâna 3 — Note și Reflecții asupra NGS și QC

## Date utilizate

**FASTQ utilizat:**
- **Accession:** SRR684066
- **Sursă:** ENA (European Nucleotide Archive)
- **Locație locală:** `data/work/AlexTGoCreative/lab03/SRR684066.fastq.gz`
- **Descriere:** Secvențe NGS TP53-related descărcate din arhiva ENA

---

## Reflecție: De ce este esențială verificarea calității datelor înainte de analiza variantelor?

### Context
Analiza variantelor genetice depinde critic de **calitatea datelor de secvențiere**. Fișierele FASTQ conțin milioane de citiri scurte (reads), fiecare cu un scor de calitate Phred asociat. Înainte de a identifica și valida variante, trebuie să asigurăm că datele sunt fiabile.

### Motive principale pentru QC

#### 1. **Erorile de secvențiere sunt sistematice**
- Primele și ultimele poziții din citiri au adesea scor Phred mai scăzut
- Bazele 'N' (ambigue) indică incertitudine în apelul bazei
- O rată mare de N → citiri de calitate scăzută care nu vor alinia corect

#### 2. **Impactul asupra alinierii și detecției variantelor**
- Citiri cu scor Phred scăzut vor alinia incorect pe genomul de referință
- Variante false pozitive (SNP-uri care sunt de fapt erori de secvențiere)
- Variante false negative (variante adevărate mascate de zgomot)

#### 3. **Costuri computaționale și timp**
- QC permite filtrarea citirilor de slabă calitate **înainte** de mapare
- Evită procesarea datelor problematice prin toată conducta
- Reducere semnificativă a timpului și resurselor

#### 4. **Validare și reproductibilitate**
- Documentarea parametrilor QC (read count, lungime medie, N-rate, Phred mediu)
- Comparație între diferite runde de secvențiere
- Detectarea problemelor în laboratorul de secvențiere (chip defect, probleme chimice)

### Metrici monitorizate în QC

Din `qc_report_AlexTGoCreative.txt` am calculat:

- **Reads:** Numărul total de citiri din eșantion
- **Mean length:** Lungimea medie a citirilor (pentru estimări de acoperire)
- **N rate:** Proporția de baze ambigue (N) — valori mari indică probleme
- **Mean Phred:** Scorul de calitate mediu — cel puțin 30 (Q30) pe bază este considerat "high quality"

### Concluzie

Verificarea calității nu este doar o etapă administrativă — **este o cerință pentru analiza fiabilă a variantelor**. O conductă NGS fără QC este ca o cercetare științifică fără control de calitate: rezultatele sunt nesigure și nereprodubile.

---

## Referințe

- Illumina FASTQ Format: https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html
- ENA Portal API: https://www.ebi.ac.uk/ena/portal/api/
- Biopython SeqIO: https://biopython.org/wiki/SeqIO
