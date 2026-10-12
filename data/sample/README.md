## Note
Fișierele FASTA de referință sunt descrise în [proveniența Lab 1/2](lab12-provenance.md). `toy_alignment.fasta` este un exemplu artificial, nu o secvență biologică. Setul proteic TP53 a fost corectat la 5 octombrie 2026.
- FASTA: secvențe TP53 (om + ortologi)
- Proteină TP53 (UniProt)
- Subset de expresie genică (GEO)

**Pentru laboratoare:** descărcați date reale extinse din NCBI/GEO/Ensembl/TCGA, conform instrucțiunilor din README-urile laboratoarelor.

### FASTQ pentru demonstrația Lab 3

`sample.fastq` păstrează cele două înregistrări complete din eșantionul existent, cu secvențele și calitățile originale nemodificate. A treia înregistrare era trunchiată, fără separator și scoruri de calitate, și a fost eliminată. Nu au fost reconstruite sau inventate calități.

Fișierul este un eșantion minimal pentru citire și calculul statisticilor, nu un set suficient pentru concluzii biologice. Sursa externă și licența datelor originale nu sunt documentate în repository; nu îl atribuiți unei gene sau unui organism pe baza numelui laboratorului.

### Expression Matrix + Metadata

- **`expr_small_tp53_geo.csv`** – subset of expression values for 8 genes across 8 samples.  
- **`expr_small_tp53_geo_metadata.csv`** – companion file giving sample-level labels.

Example:

```csv
SampleID,Condition
Sample1,Tumor
Sample2,Tumor
Sample3,Tumor
Sample4,Normal
```
How to use:

Merge the expression table with the metadata by SampleID.

Use Condition as the label column for:

Clustering → compare Tumor vs Normal grouping.

Supervised ML (RF, SVM, FL labs) → train/predict tumor status.

This file is intentionally small and balanced (toy example).
For assignments, students must fetch and annotate real GEO/TCGA metadata.
