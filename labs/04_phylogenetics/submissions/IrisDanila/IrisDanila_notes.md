# Săptămâna 4 — Filogenetică  
Autor: IrisDanila

## Date folosite
Am folosit secvențele TP53 mRNA descărcate în Lab 1:
- **Fișier**: `data/work/IrisDanila/lab04/your_sequences.fasta`
- **Sursa originală**: NCBI Nucleotide Database
- **Secvențe**: 5 variante de transcript TP53 (Homo sapiens)
  - NM_001126117.2 (variant 7) - 2063 bp
  - NM_001126115.2 (variant 5) - 2003 bp
  - NM_001126118.2 (variant 8) - 2629 bp
  - NM_001126112.3 (variant 2) - 2509 bp
  - NM_001126116.2 (variant 6) - 2136 bp

Am rulat următoarele scripturi:
```bash
python demo01_distance_matrix.py
python ex01_phylo_NJ.py


Rezultate și observații
Matrice de distanțe (rezultat cod):

Sequences: ['NM_001126117.2', 'NM_001126115.2', 'NM_001126118.2', 'NM_001126112.3', 'NM_001126116.2']
Distance Matrix (p-distance):
[[0.         0.51023465 0.7353969  0.74288567 0.50374438]
 [0.51023465 0.         0.74438342 0.73040439 0.51522716]
 [0.7353969  0.74438342 0.         0.70394408 0.74088867]
 [0.74288567 0.73040439 0.70394408 0.         0.7324014 ]
 [0.50374438 0.51522716 0.74088867 0.7324014  0.        ]]

Calculată folosind p-distance (proporția pozițiilor diferite)
Secvențele truncate la lungimea minimă (2003 bp) pentru comparație uniformă
Majoritatea variantelor TP53 arată similaritate mare (distanțe mici)
Arbore filogenetic (Neighbor-Joining):

Salvat în format Newick: tree_IrisDanila.nwk
Metoda NJ grupează secvențele după similaritate
Variantele apropiate evolutiv formează clustere


Observații:

Toate secvențele sunt din aceeași genă (TP53) și organism (H. sapiens)
Diferențele reflectă splicing alternativ, nu evoluție divergentă
Arborele arată relațiile structurale între isoformele TP53
Reflecție: Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?
Arborele filogenetic oferă mai multe avantaje semnificative față de matricea de distanțe:

1. Vizualizare ierarhică a relațiilor

Matricea de distanțe arată doar valori numerice perechi
Arborele prezintă vizual structura ierarhică și grupările
Permite identificarea rapidă a clusterelor și outlier-ilor
Exemplu: În cazul meu, variantele 5 și 7 pot forma un cluster distinct
2. Reconstrucția istoriei evolutive

Arborele arată ordinea ramificărilor (topologie)
Lungimile ramurilor indică gradul de divergență
Nodurile interne reprezintă strămoși comuni inferați
Matricea nu oferă context temporal sau direcționalitate
3. Identificarea grupurilor funcționale

Secvențele înrudite tind să aibă funcții similare
Clustering-ul din arbore sugerează grupuri funcționale
Pentru TP53: variante cu structură sau localizare similară se grupează
Utilitate în prezicerea funcției pentru secvențe noi
4. Reducerea complexității informației

Matricea n×n devine greu de interpretat pentru n mare
Arborele condensează n(n-1)/2 distanțe într-o structură simplă
Mai ușor de comunicat și publicat în articole științifice
Exemplu: 5 secvențe → 10 distanțe vs. 1 arbore clar
5. Detectarea erorilor și outlier-ilor

Secvențe contaminante sau greșit adnotate ies în evidență
Poziții neașteptate în arbore semnalează probleme
Matricea nu oferă acest context vizual imediat
6. Suport statistic pentru grupări

Arborele permite bootstrap values sau posterior probabilities
Indică încrederea în fiecare ramură
Matricea nu oferă validare statistică
7. Compararea cu cunoștințe biologice

Arborele poate fi comparat cu taxonomie sau funcții cunoscute
Discrepanțe sugerează transfer genetic orizontal sau convergență
În cazul TP53: arborele trebuie să reflecte relațiile de splicing cunoscute
8. Predicții biologice

Caractere ancestrale pot fi reconstruite la noduri interne
Permite datarea divergențelor (cu ceas molecular)
Identifică evenimente de duplicare genică sau pierdere de gene
Aplicație în cazul meu (TP53 variants):

Deși lucrez cu variante de splicing (nu evoluție divergentă clasică), arborele NJ oferă totuși informații valoroase:

Gruparea variantelor similare structural: Variante care includ aceiași exoni se vor grupa
Identificarea variantelor predominante: Varianta "centrală" în arbore poate fi isoforma majoră
Relații cu funcționalitate: Variante cu domenii conservate similare vor fi apropiate
Validare a datelor: Dacă o variantă e foarte îndepărtată, poate indica eroare de anotare
Concluzie:

În timp ce matricea de distanțe oferă date brute cantitative, arborele filogenetic transformă aceste date într-o ipoteză vizuală, testabilă despre relațiile evolutive. Este instrumentul fundamental pentru:

Taxonomie și sistematică
Genomica comparativă
Identificarea genelor ortoloage și paraloage
Studii de evoluție moleculară
Bioinformatică structurală