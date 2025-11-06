**Daniela Ispas**
**Eugenia Vacarciuc**

**demo01**
python3 labs/04_phylogenetics/demo01_distance_matrix.py
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]

 **ex01**
python labs/04_phylogenetics/submissions/eeeevn0/ex01_phylo_NJ.py
Am încărcat 3 secvențe din data/work/eeeevn0/lab04/my_sequences.fasta
ID-uri secvențe: SRR000001.1, SRR000001.4, SRR000001.5

Secvențele AU lungimi diferite: [266, 240, 276]
Tăiem toate secvențele la lungimea minimă: 240 baze

Matricea de distanțe (identity):
SRR000001.1 0.000000
SRR000001.4 0.775000    0.000000
SRR000001.5 0.700000    0.758333    0.000000
    SRR000001.1 SRR000001.4 SRR000001.5

Arbore NJ salvat în: labs/04_phylogenetics/submissions/eeeevn0/tree_eeeevn0.nwk

Arbore (ASCII):
  _______________________________________________________ SRR000001.1
 |
_|_________________________________________________________________ SRR000001.4
 |
 |_____________________________________________________ SRR000001.5


Arbore (Newick):
(SRR000001.1:0.35833,SRR000001.4:0.41667,SRR000001.5:0.34167)Inner1:0.00000;


Am folosit trei secvențe extrase din fișierul `SRR000001.fastq.gz`, pe care l-am folosit la laboratorul anterior(lab3).
Le-am convertit în format FASTA și le-am salvat într-un fișier numit `my_sequences.fasta`.  

Secvențele sunt:
- SRR000001.1 – 266 baze  
- SRR000001.4 – 240 baze  
- SRR000001.5 – 276 baze  

**Reflecție**
Matricea de distanțe arată doar cât de diferite sunt secvențele,pe când arborele arată și relația dintre ele, ca un fel de „arbore de familie” al ADN-ului, fapt, ce conturează clar cum se construiește și interpretează un arbore filogenetic.