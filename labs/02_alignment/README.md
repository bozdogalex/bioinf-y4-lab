# Săptămâna 2 — Sequence Alignment

## Scopuri
- Înțelegerea tipurilor de aliniere: global (Needleman–Wunsch), local (Smith–Waterman), semiglobal.  
- Folosirea matricilor de substituție (PAM, BLOSUM).  
- Exersarea cu Biopython și instrumente externe (BLAST, Clustal Omega).  
- Construirea competențelor de implementare a algoritmilor de bază pentru aliniere.

---

## Partea 1 — Demo / Exerciții
**Rulați**  
- `demo01_pairwise_biopython.py` — aliniere globală și locală cu Biopython (pairwise2).  
- `demo02_distance_matrix.py` — calcul distanțe (p-distance, Hamming) pe FASTA.  

**Completați și rulați**  
- `ex02_global_nw.py` — schelet de implementare pentru aliniere globală (TODO).  
- `ex03_local_sw.py` — schelet de implementare pentru aliniere locală (TODO).  

Notă: folosiți datele descărcate în Lab 1 (din `data/work/<handle>/lab01/`)

---

## Livrabile
În PR trebuie să apară:
1. Fișierul `labs/02_alignment/submissions/<github_handle>_notes.md` cu:  
   - ce date ați folosit (ex. TP53 vs. BRCA1),  
   - o scurtă reflecție: **Când este de preferat alinierea globală vs. locală?”**  
2. Exercițiile completate, salvate în:  
   ```bash
   labs/02_alignment/submissions/<github_handle>/ex02_global_nw.py
   labs/02_alignment/submissions/<github_handle>/ex03_local_sw.py
   ```
3. Completarea checklist-ului din șablonul PR.

---

### Săptămâna următoare
- Vom extinde analiza la citiri NGS (FASTQ → mapare → variant calling).
- Vom folosi alinierile obținute pentru a valida maparea și analizele NGS
- [Vezi Săptămâna 3 — NGS Analysis](../03_formats&NGS/README.md)

---

## Competențe
- Înțelegerea diferenței dintre aliniere globală și locală.
- Utilizarea Biopython pentru aliniere simplă.
- Implementarea de bază a algoritmilor NW și SW.
- Interpretarea rezultatelor și compararea cu BLAST/Clustal.

---

## Resurse 
- [Fișa laborator](../../docs/lab_onepagers/02_alignment.md)  
- [Aliniere globală (Needleman–Wunsch)](../../docs/presentations/alignment1.pdf)  
- [Aliniere locală (Smith–Waterman)](../../docs/presentations/alignment2.pdf)  
- [Applied Bioinformatics of Nucleic Acids — Cap. 1](../../docs/papers/Applied_Bioinformatics.pdf)  
- [Scoring Matrix Development (BLOSUM62) (pdf în /papers)](../../docs/papers/Scoring_matrix_development_BLOSUM62.pdf)  
- Substitution matrices: [BLOSUM62 (NCBI)](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62)  
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)  
- [Clustal Omega — Multiple Sequence Alignment](https://www.ebi.ac.uk/Tools/msa/clustalo/)  
- [Biopython pairwise2](https://biopython.org/docs/1.75/api/Bio.pairwise2.html)  
- [Earth bio genome project](https://www.earthbiogenome.org/)
- [Open Tree of Life](https://tree.opentreeoflife.org/)