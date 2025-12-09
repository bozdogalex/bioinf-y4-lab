# Săptămâna 3 — Formate și NGS  
Autor: IrisDanila

## Date folosite
Am folosit fișierul FASTQ real descărcat de la ENA (European Nucleotide Archive):
- **Accession**: SRR000001
- **Sursă**: http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000001/SRR000001.fastq.gz
- **Dimensiune**: 53.21 MB
- **Locație**: `data/work/IrisDanila/lab03/your_reads.fastq.gz`

Am rulat următoarele scripturi:
```bash
python ex01_fetch_fastq.py
python ex02_fastq_stats.py

Rezultate și observații
Scriptul ex02_fastq_stats.py a analizat fișierul FASTQ și a generat următoarele statistici:

Rezultate QC:

Total citiri: 236,003 reads
Lungime medie: 254.62 bp
N-rate: 0.0005 (0.05%)
Scor Phred mediu: 25.96
Evaluare calitate:

PASS - Low N-rate (< 1%)



Observații:

Structura FASTQ: 4 linii per citire (ID @, secvență, +, calitate)
Scoruri Phred encodate ASCII (Q = ord(char) - 33)
Dataset real NGS cu ~236K citiri
Calitate moderată dar acceptabilă pentru analiză (Phred ~26)
N-rate foarte scăzut (0.05%) indică date de bună calitate
Reflecție: De ce este esențială verificarea calității datelor înainte de analiza variantelor?
Verificarea calității datelor NGS este crucială din următoarele motive:

1. Prevenirea variantelor false pozitive

Scorurile Phred scăzute (<20) indică probabilitate mare de eroare de secvențiere
Un Phred de 20 = 1% probabilitate de eroare; Phred 30 = 0.1% eroare
În cazul meu (Phred mediu 25.96), există ~0.25% șansă de eroare per bază
Fără filtrare, aceste erori par variante reale în analiză
2. Detectarea problemelor tehnice timpurii

N-rate ridicat (>1-2%) semnalează degradare ADN sau contaminare
Dataset-ul meu (0.05% N-rate) este foarte bun
Distribuții anormale ale Phred pot indica probleme cu chimia de secvențiere
QC permite re-secvențierea înainte de cheltuieli mari de analiză
3. Optimizarea parametrilor bioinformatici

Lungimea medie (254.62 bp) ghidează strategia de aliniere
Calitatea influențează toleranța la mismatch-uri în mapare
Permite setarea corectă a threshold-urilor pentru variant calling
Ex: cu Phred 26, pot tolera mai puține mismatch-uri decât cu Phred 35
4. Eficiență computațională

Eliminarea citirii de calitate scăzută reduce timpul de procesare
Trimming-ul capetelor cu Phred scăzut îmbunătățește alinierea
În cazul meu: 236K citiri × 255 bp = ~60 Mb date de procesat
Filtrarea economisește storage și costuri cloud
5. Reproductibilitatea și transparența

QC standardizat permite comparații între experimente
Documentarea calității este obligatorie în publicații științifice
Alți cercetători pot evalua fiabilitatea datelor
Criterii QC clare asigură rezultate reproductibile
6. Impact clinic și diagnostic

În genomica clinică, variantele false duc la diagnostic greșit
QC-ul este obligatoriu în pipeline-uri certificate (ISO, CAP, CLIA)
Salvează vieți evitând decizii medicale bazate pe artefacte
Reduce costuri evitând validări inutile prin Sanger sequencing
7. Exemplu concret din datele mele:

Cu 236,003 citiri × 254.62 bp = ~60 million de baze
La Phred 26 (0.25% eroare) → ~150,000 baze potențial greșite
Fără QC, aceste 150K erori pot părea variante reale!
Filtrarea previne zeci de mii de false positive



Concluzie:

QC-ul este prima și cea mai importantă barieră împotriva datelor de calitate scăzută. Fără el, întregul pipeline bioinformatic devine nesigur, iar concluziile biologice sunt compromise.

În NGS, regula de aur este: "garbage in, garbage out" - QC-ul asigură că intrarea este de înaltă calitate, rezultând în variante fiabile și reproducibile.