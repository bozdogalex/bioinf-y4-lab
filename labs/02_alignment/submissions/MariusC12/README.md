# Lab 2 – Sequence Alignment

**Student:** MariusC12

## 1. Descriere

Proiectul analizează secvențe nucleotidice ale genei TP53 de la trei organisme: Homo sapiens, Mus musculus și Danio rerio.

Au fost realizate:
- Calculul distanțelor p între secvențe.
- Implementarea algoritmilor Needleman–Wunsch și Smith–Waterman.
- Compararea rezultatelor implementărilor proprii cu Biopython.
- Alinierea multiplă a celor trei secvențe cu Clustal Omega.
- Identificarea unei regiuni conservate.

## 2. Date utilizate

| Accession | Organism | Tip moleculă |
|---|---|---|
| NM_000546.6 | Homo sapiens | mRNA |
| NM_011640.3 | Mus musculus | mRNA |
| NM_131327.2 | Danio rerio | mRNA |

Secvențele au fost descărcate din baza de date NCBI Nucleotide la data de 06.11.2026.

Fișierele FASTA sunt salvate în `data/work/MariusC12/lab01/`.

### Descărcarea secvențelor

Secvențele pot fi descărcate din nou folosind NCBI E-utilities:

```bash
mkdir -p data/work/MariusC12/lab01

curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_000546.6&rettype=fasta&retmode=text" -o data/work/MariusC12/lab01/seq1.fasta

curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_011640.3&rettype=fasta&retmode=text" -o data/work/MariusC12/lab01/seq2.fasta

curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_131327.2&rettype=fasta&retmode=text" -o data/work/MariusC12/lab01/seq3.fasta

cat data/work/MariusC12/lab01/seq1.fasta data/work/MariusC12/lab01/seq2.fasta data/work/MariusC12/lab01/seq3.fasta > data/work/MariusC12/lab01/sequences.fasta
```

## 3. Structura proiectului

- `task1_distance.py` – calculul distanțelor p.
- `task1_distances.csv` – rezultatele distanțelor.
- `ex01_global_nw.py` – implementarea Needleman–Wunsch.
- `ex02_local_sw.py` – implementarea Smith–Waterman.
- `task2_regions.py` – extragerea regiunilor CDS și selectarea fragmentelor.
- `task2_regions.fasta` – cele două fragmente utilizate.
- `task2_compare.py` – compararea implementărilor proprii cu Biopython.
- `task3_msa.fasta` – alinierea multiplă completă.
- `task3_msa_excerpt.fasta` – fragmentul conservat selectat.
- `task3_analyze_msa.py` – analiza regiunii conservate.

## 4. Rulare

Calcularea distanțelor:

```bash
python "labs/02_alignment/submissions/MariusC12/task1_distance.py" --fasta "data/work/MariusC12/lab01/sequences.fasta" > "labs/02_alignment/submissions/MariusC12/task1_distances.csv"
```

Extragerea regiunilor:

```bash
python "labs/02_alignment/submissions/MariusC12/task2_regions.py"
```

Compararea algoritmilor:

```bash
python "labs/02_alignment/submissions/MariusC12/task2_compare.py"
```

Alinierea multiplă:

```bash
clustalo -i "data/work/MariusC12/lab01/sequences.fasta" -o "labs/02_alignment/submissions/MariusC12/task3_msa.fasta" --outfmt=fasta --force
```

Analiza regiunii conservate:

```bash
python "labs/02_alignment/submissions/MariusC12/task3_analyze_msa.py"
```

## 5. Parametri de aliniere

Pentru comparația algoritmilor au fost folosite:
- Match: +1
- Mismatch: -1
- Gap: -2

## 6. Versiuni software

- Python: 3.11.18
- Biopython: 1.85
- Clustal Omega: 1.2.4

## 7. Surse și utilizarea AI

Sursa secvențelor: NCBI Nucleotide.

Biopython și Clustal Omega au fost utilizate pentru aliniere și verificarea rezultatelor. ChatGPT a fost utilizat pentru explicații, depanarea codului și sprijin în interpretarea rezultatelor.