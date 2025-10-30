# Săptămâna 3 — Assignment (Formate și NGS)

## Obiectiv
Simulați un flux simplificat de analiză NGS folosind instrumente programatice (Biopython, requests sau API-uri publice) și fișiere proprii.  
Tema are rolul de a integra cunoștințele despre baze de date, formate și QC, și de a face legătura dintre datele brute și literatura științifică.

---

## Cerințe 

### 1) PubMed parsing — **2p**
- Scrieți `pubmed_query.py` care:
  - Interoghează PubMed cu `TP53 AND cancer` (max. 5 articole).  
  - Salvează titlul, autorii și rezumatul într-un fișier `pubmed_<handle>.txt`.

### 2) FASTA (din Lab 1) — **2p**
- Refolosiți fișierul FASTA cu secvența TP53 descărcat în Lab 1.  
- Salvați o copie ca `fasta_<handle>.fasta`.

### 3) FASTQ QC — **3p**
- Folositi codul din `ex04_fastq_stats.py` pentru a calcula:
  - număr citiri, lungime medie, proporție de N, scor Phred mediu.  
- Folosiți fișierul stocat în `data/work/<handle>/lab03/your_reads.fastq.gz`.  
- Salvați raportul ca `qc_report_<handle>.txt`.  
- **Dacă folosiți fișierul FASTQ de rezervă:** max **2p** (din 3).

### 4) VCF → PubMed — **3p**
- Scrieți `vcf_pubmed.py` care:
  - Citește un fișier VCF și extrage cel puțin 2 variante.  
  - Pentru fiecare variantă:
    - Dacă are ID (ex. `rs12345`) → căutați pe PubMed acel rsID.  
    - Dacă nu are ID → căutați `chr<CHROM>:<POS> AND TP53`.  
  - Salvați rezultatele ca `variants_<handle>.txt`.  
- **Dacă folosiți fișierul VCF sample:** max **2p** 

---

### Bonus (+1p)
- Generați o vizualizare simplă (matplotlib) cu distribuția lungimilor si scorurilor Phred și salvați ca `qc_plot_<handle>.png`.

---

## Livrabile (upload pe CV)
Încărcați pe Moodle un fișier `.zip` care să conțină:  
- `pubmed_query.py`, `vcf_pubmed.py`, `ex04_fastq_stats.py`  
- `pubmed_<handle>.txt`, `fasta_<handle>.fasta`, `qc_report_<handle>.txt`, `variants_<handle>.txt`  
- `qc_plot_<handle>.png` (opțional, bonus)  
- `README.txt` și `notes.pdf`

---

### `README.txt` trebuie să conțină:
- Autori (1–2) + handle-uri GitHub.  
- Disclosure AI (ce ați folosit, cum ați verificat).  
- Sursa FASTQ 

### `notes.pdf` (max 1 pagină):
- 3–5 propoziții: „De ce este important QC înainte de variant calling?”  
- 2–3 propoziții: cum ați formulat căutările PubMed pentru variante.
