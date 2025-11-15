LAB 02 - SEQUENCE ALIGNMENT
===========================

Autor: AlexTGoCreative
Data: 15 noiembrie 2025

DESCRIERE
---------
Acest pachet conține soluțiile pentru assignment-ul Lab 02 (Sequence Alignment):
- Task 1: Matrice de distanțe Hamming
- Task 2: Aliniamente pairwise (global vs local)
- Task 3: Multiple sequence alignment (MSA) cu Clustal Omega
- Bonus: Aliniament semiglobal

CERINȚE SISTEM
--------------
- Python 3.8 sau mai recent
- Biopython >= 1.79

INSTALARE DEPENDENȚE
--------------------
pip install biopython

STRUCTURA FIȘIERELOR
-------------------
lab02_solution/
├── task01_hamming_distance.py      # Task 1: Distanțe Hamming
├── task02_pairwise_alignments.py   # Task 2: Global vs Local
├── task03_msa_guide.py             # Task 3: MSA + analiză
├── bonus_semiglobal.py             # Bonus: Semiglobal alignment
├── tp53_protein_multi.fasta        # Date input (secvențe TP53)
├── sequences_for_msa.fasta         # Secvențe pentru MSA online
├── msa_result.aln                  # Rezultat MSA de la Clustal Omega
├── notes.md                        # Notițe detaliate (convertit în PDF)
└── README.txt                      # Acest fișier

RULARE
------
Toate scripturile se rulează din directorul curent (self-contained):

1. Task 1 - Hamming distance:
   python task01_hamming_distance.py

2. Task 2 - Pairwise alignments:
   python task02_pairwise_alignments.py

3. Task 3 - MSA guide și analiză:
   python task03_msa_guide.py

4. Bonus - Semiglobal alignment:
   python bonus_semiglobal.py

NOTE IMPORTANTE
--------------
- Toate fișierele necesare sunt incluse în arhivă
- Nu sunt necesare path-uri externe
- Rezultatul MSA (msa_result.aln) este deja inclus
- Pentru notes.pdf, convertește notes.md cu:
  pandoc notes.md -o notes.pdf

RESURSE EXTERNE
--------------
- GitHub Copilot: asistență în cod
- Clustal Omega (EBI): https://www.ebi.ac.uk/Tools/msa/clustalo/
- Documentație Biopython: https://biopython.org/
