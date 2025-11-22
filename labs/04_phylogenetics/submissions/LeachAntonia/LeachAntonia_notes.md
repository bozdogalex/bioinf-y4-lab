1. Secventele descărcate de la ENA:
Accession	Specie / Descriere
AF209148	Human TP53
AF151353	Mouse TP53
AF209133	Human TP53 variant

Am folosit create_multifasta.py pentru a crea un multi-FASTA, iar înainte de a lucra cu fișierul FASTA, am făcut padding pentru a alinia cele 3 secvențe la aceeași lungime.

2. Reflectie
Matricea de distanțe ne arată doar cât de diferite sunt secvențele între ele, ca valori numerice.
Distance matrix:
ENA|AF209148|AF209148.1 0.000000
ENA|AF151353|AF151353.1 0.966643    0.000000
ENA|AF209133|AF209133.1 0.144074    0.972321    0.000000
    ENA|AF209148|AF209148.1 ENA|AF151353|AF151353.1 ENA|AF209133|AF209133.1

Arborele filogenetic (NJ) organizează aceste distanțe într-o structură ierarhică care arată relațiile evolutive: cine este mai apropiat de cine, care sunt grupurile de secvențe înrudite, și oferă o vizualizare clară a legăturilor între specii sau variante.
ASCII view of NJ tree:
  ___ ENA|AF209148|AF209148.1
 |
_|_____________________________________________________ ENA|AF151353|AF151353.1
 |
 |___ ENA|AF209133|AF209133.1


Arborele permite să observăm clade, divergențe și posibile strămoși comuni, ceea ce nu poate fi intuit dintr-o simplă matrice numerică.