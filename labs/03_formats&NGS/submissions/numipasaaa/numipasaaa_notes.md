# Demo 02 - Mapping Toy
```
Read ATGCTAGC mapped at position 0
Read GATCGATC mapped at position 19
Read TACGATCG mapped at position 28
Read GGGGGGGG did not map
```

# Exercitiu 1 - Fetch FASTQ
```
Descarcare din: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR305/058/SRR30503258/SRR30503258_1.fastq.gz
Downloaded: data/work/numipasaaa/lab03/SRR30503258_1.fastq.gz
```

# Exercitiu 2 - FASTQ Stats
```
[OK] QC report -> /workspaces/bioinf-y4-lab/labs/03_formats&NGS/submissions/numipasaaa/qc_report_numipasaaa.txt
```

# Livrabile
- ce FASTQ ați folosit (link sau accession SRA)
  - **Accession**: SSR30503258
- o scurtă reflecție: **De ce este esențială verificarea calității datelor înainte de analiza variantelor?**
  - **Verificarea calitatii** (QC) este esentiala din cauza principiului **"Garbage In, Garbage Out"** (GIGO). Daca datele de secventiere sunt de proasta calitate, analiza variantelor va produce rezultate fundamental gresite.
  - Fara QC, risti trei probleme majore:
        
        1. Fals-Pozitive: Erori de citire care pot fi interpretate gresit.

        2. Fals-Negative: Datele de calitate slaba sunt adesea respinse in timpul analizei. Daca o mutatie reala se afla doar pe acele citiri, ea va fi ratata.

        3. Bias: De exemplu, regiuni duplicate.