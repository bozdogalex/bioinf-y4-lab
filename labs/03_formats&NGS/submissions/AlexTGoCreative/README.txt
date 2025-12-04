================================================================================
README - Lab 03: Formate și NGS Assignment
================================================================================

AUTORI:
-------
- GitHub Handle: AlexTGoCreative

DISCLOSURE AI:
--------------
Pentru realizarea acestui assignment am folosit următoarele instrumente AI:
- GitHub Copilot pentru generarea codului Python și scripturi
- ChatGPT pentru verificarea sintaxei și debugging
- Am verificat toate rezultatele manual și am testat fiecare script individual
- Am înțeles logica fiecărui script și pot explica funcționarea acestora
- Tot codul a fost revizuit și adaptat pentru cerințele specifice ale assignment-ului

SURSA FASTQ:
------------
Fișierul FASTQ utilizat: SRR684066.fastq.gz
Sursa: NCBI Sequence Read Archive (SRA)
Locație: /workspaces/bioinf-y4-lab/data/work/AlexTGoCreative/lab03/SRR684066.fastq.gz
Descriere: Secvențe NGS pentru gena TP53

STRUCTURA FIȘIERELOR:
---------------------
1. pubmed_query.py         - Script pentru interogarea PubMed
2. ex02_fastq_stats.py     - Script pentru analiza QC a fișierului FASTQ
3. vcf_pubmed.py           - Script pentru parsarea VCF și căutare PubMed
4. fasta_AlexTGoCreative.fasta - Fișier FASTA cu secvența TP53
5. pubmed_AlexTGoCreative.txt  - Rezultate interogare PubMed pentru "TP53 AND cancer"
6. qc_report_AlexTGoCreative.txt - Raport QC pentru FASTQ
7. qc_plot_AlexTGoCreative.png - Vizualizare distribuții (BONUS)
8. variants_AlexTGoCreative.txt - Rezultate căutare PubMed pentru variante VCF
9. README.txt              - Acest fișier
10. notes.pdf              - Notițe despre importanța QC și strategii de căutare

INSTRUCȚIUNI DE RULARE:
------------------------
1. PubMed Query:
   python pubmed_query.py
   
2. FASTQ QC Analysis:
   python ex02_fastq_stats.py
   
3. VCF to PubMed:
   python vcf_pubmed.py

REZULTATE OBȚINUTE:
-------------------
- 5 articole despre TP53 și cancer din PubMed
- Statistici QC complete pentru fișierul FASTQ (număr citiri, lungime medie, 
  proporție N, scor Phred mediu)
- Visualizări pentru distribuția lungimilor de citiri și scorurilor Phred (BONUS)
- Căutări PubMed pentru variante din fișierul VCF (rsID sau poziție cromozomială)
