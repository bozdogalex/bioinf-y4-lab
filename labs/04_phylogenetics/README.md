# Săptămâna 4 — Filogenetică

## Scopuri
- Înțelegerea conceptelor de bază în filogenetică.  
- Calcularea distanțelor între secvențe și construirea arborilor filogenetici.  
- Exersarea metodei **Neighbor-Joining (NJ)** cu Biopython.  
- Consolidarea legăturii dintre aliniamente multiple și relațiile evolutive.  

---

## Context
După ce în săptămâna 3 am explorat formatele și analiza datelor NGS, acum trecem la **arbori filogenetici**.
Pornim de la un set de secvențe (multi-FASTA), calculăm o matrice de distanțe și construim un arbore NJ.
Vom compara rezultatele cu un MSA realizat online pentru a observa regiunile conservate.  

---

## Hands-on — Distanțe și arbori filogenetici
**Rulați**  
- `demo01_distance_matrix.py` — calculul unei matrice de distanțe pe un multi-FASTA (Hamming/p-distance).  

**Completați și rulați**  
- `ex05_phylo_tree.py` — exercițiu:  
  - salvați rezultatul ca `.nwk` în `labs/04_phylogenetics/submissions/<handle>/tree_<handle>.nwk`.

---

## Livrabile
În PR trebuie să apară:
1. Fișierul `labs/04_phylogenetics/submissions/<github_handle>_notes.md` cu:  
   - ce secvențe FASTA ați folosit (link/descriere),  
   - o scurtă reflecție: **Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?**  
2. Exercițiul completat, salvat în:  
   ```bash
   labs/04_phylogenetics/submissions/<github_handle>/ex05_phylo_tree.py
   labs/04_phylogenetics/submissions/<github_handle>/tree_<handle>.nwk
   ```
3. Completarea checklist-ului din șablonul PR.

---

## Săptămâna următoare
- Clustering de expresie genică și analiza grupurilor.
- De la arbori evolutivi la module de co-expresie.
- [Vezi Săptămâna 5 — Clustering](../04_phylogenetics/README.md)

---

## Competențe
- Calculul matricelor de distanțe dintr-un multi-FASTA.
- Construirea și vizualizarea arborilor NJ cu Biopython.
- Interpretarea clusterelor evolutive.
- Conectarea alinierilor multiple la relațiile filogenetice.

---

## Resurse
- [Fișa laborator](../../docs/lab_onepagers/04_phylogenetics.md)  
- [Biopython Phylo](https://biopython.org/wiki/Phylo)
- [Formatul Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html)
- [Clustal Omega (MSA online)](https://www.ebi.ac.uk/Tools/msa/clustalo/)
- [Multiple Sequence Alignment în Biopython](https://biopython.org/wiki/AlignIO)