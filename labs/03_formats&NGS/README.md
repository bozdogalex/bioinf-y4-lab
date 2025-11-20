# Săptămâna 3 — Formate și Next-Generation Sequencing (NGS)

## Scopuri
- Înțelegerea formatelor utilizate în bioinformatică: **FASTA, FASTQ, SAM, VCF**.  
- Exersarea verificării calității datelor NGS (QC: read count, lungime, N-rate, Phred).  
- Familiarizarea cu pașii dintr-un flux de analiză NGS.  
- Consolidarea legăturii dintre date biologice și literatura științifică.  

---

## Context
După ce în săptămâna 2 am studiat alinierea secvențelor (globală și locală), acum trecem la **date NGS**, unde lucrăm cu milioane de citiri scurte stocate în format FASTQ. Vom învăța să verificăm calitatea datelor și să conectăm rezultatele la articole științifice.

---

## Partea 1 — File formats & NGS
**Rulați**  
- `demo01_fastq_qc.py` — citirea unui fișier FASTQ și calcularea statisticilor de bază (număr citiri, lungime medie, N-rate, scor Phred mediu).  
- `demo02_mapping_toy.py` — exemplu simplificat de mapare a citirilor pe o secvență de referință (potrivire exactă de șiruri).  

**Completați și rulați**  
- `ex01_fetch_fastq.py` — exercițiu scurt: descărcați **un FASTQ propriu** (TP53-related) folosind API-ul ENA și salvați-l în `data/work/<handle>/lab03/your_reads.fastq.gz`.  
- `ex02_fastq_stats.py` — calculați QC pe fișierul vostru FASTQ (acceptă `.fastq` sau `.fastq.gz`) și salvați raportul în `labs/03_formats&NGS/submissions/<handle>/qc_report_<handle>.txt`.

## Livrabile
În PR trebuie să apară:
1. Fișierul `labs/03_formats&NGS/submissions/<github_handle>_notes.md` cu:  
   - ce FASTQ ați folosit (link sau accession SRA),  
   - o scurtă reflecție: **De ce este esențială verificarea calității datelor înainte de analiza variantelor?**  
2. Exercițiul completat, salvat în:  
   ```bash
   labs/03_formats&NGS/submissions/<github_handle>/ex04_fastq_stats.py
   labs/03_formats&NGS/submissions/<github_handle>/qc_report_<handle>.txt
   ```
3. Completarea checklist-ului din șablonul PR.

## Săptămâna următoare

- Arbori filogenetici folosind distanțe între secvențe și aliniamente multiple.
- De la analiza variantelor la relațiile evolutive.
- [Vezi Săptămâna 4 — Filogenetică](../04_phylogenetics/README.md)

## Competențe

- Înțelegerea diferenței între formatele FASTA, FASTQ, SAM, VCF.
- Utilizarea Biopython pentru citirea și procesarea fișierelor FASTQ.
- Implementarea verificărilor de bază pentru QC.
- Conectarea datelor NGS la informații din literatura științifică.

## Resurse

- [Fișa laborator](../../docs/lab_onepagers/02_alignment.md)
- [FASTQ (Illumina)](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)
- [Biopython SeqIO](https://biopython.org/wiki/SeqIO)
- [Formatul SAM/VCF (htslib)](http://samtools.github.io/hts-specs/)
- [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) — căutare și filtrare rulari SRA.  
- [Entrez Programming Utilities (E-utilities)](https://www.ncbi.nlm.nih.gov/books/NBK25501/) — interfață pentru metadate și linkuri.  
- [NCBI Datasets API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) — API modern pentru secvențe și rulari.
### **ENA (European Nucleotide Archive)**
- [ENA Browser](https://www.ebi.ac.uk/ena/browser/home) — interfață web pentru genomuri și rulari.  
- [ENA Portal API (filereport)](https://www.ebi.ac.uk/ena/portal/api/) — acces programatic la linkurile FASTQ.
- [Biopython Entrez](https://biopython.org/docs/1.75/api/Bio.Entrez.html) — pentru acces la NCBI 