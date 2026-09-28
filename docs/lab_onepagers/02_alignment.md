# Săptămâna 2 — Sequence Alignment

## Obiective
- Înțelegerea tipurilor de aliniere: global (Needleman–Wunsch), local (Smith–Waterman), semiglobal.  
- Aplicarea matricilor de substituție (PAM, BLOSUM).  
- Exersarea cu Biopython și compararea cu unelte externe (BLAST, Clustal Omega).  
- Implementarea de bază a algoritmilor NW și SW.

---

## Pași principali
1. **Rularea demo-urilor**: pairwise Biopython și calcul distanțe (Hamming, p-distance).  
2. **Completarea exercițiilor**:  
   - `ex01_global_nw.py` → aliniere globală (Needleman–Wunsch).  
   - `ex02_local_sw.py` → aliniere locală (Smith–Waterman).  

---

## Competențe 
- Diferențierea între aliniere globală și locală.  
- Utilizarea Biopython pentru aliniere rapidă.  
- Înțelegerea și implementarea algoritmilor NW și SW.  
- Interpretarea rezultatelor și compararea cu unelte consacrate (BLAST, Clustal).  

---

## Resurse
- [Aliniere globală (Needleman–Wunsch)](../../docs/presentations/02_alignment1.pdf)  
- [Aliniere locală (Smith–Waterman)](../../docs/presentations/02_alignment2.pdf)  
- [Applied Bioinformatics of Nucleic Acids — Cap. 1](../../docs/papers/Applied_Bioinformatics.pdf)  
- [Scoring Matrix Development (BLOSUM62) (pdf în /papers)](../../docs/papers/Scoring_matrix_development_BLOSUM62.pdf)  
- Substitution matrices: [BLOSUM62 (NCBI)](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62)  
- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)  
- [Clustal Omega — Multiple Sequence Alignment](https://www.ebi.ac.uk/Tools/msa/clustalo/)  
- [Biopython pairwise2](https://biopython.org/docs/1.75/api/Bio.pairwise2.html)  
