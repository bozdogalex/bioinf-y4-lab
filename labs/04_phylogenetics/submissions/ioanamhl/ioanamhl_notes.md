
**demo01**
root@codespaces-ad53c2:/workspaces/bioinf-y4-lab# python 'labs/04_phylogenetics/demo01_distance_matrix.py'
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]


**Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?**
O matrice de distanțe arată doar cât de diferite sunt secvențele între ele — valori numerice izolate, fără context evolutiv.În schimb, un arbore filogenetic organizează aceste relații într-o structură ierarhică, arătând: modul în care secvențele se grupează (care sunt mai apropiate și formează clade comune), ordinea divergenței — adică cine a avut un strămoș comun mai recent, și uneori chiar o estimare a lungimii ramurilor proporțională cu numărul de mutații acumulate.

**fisier folosit**
Am folosit fisierul your_sequences.fasta, generat din datele NGS proprii (your_reads.fastq.gz) obtinute în laboratorul 3.
Am extras 10 citiri din secventierea asociata genei TP53, folosind un script Python care a convertit primele citiri din FASTQ în FASTA.
Acest fisier conține 10 secvențe scurte de ADN, fiecare reprezentând fragmente ale aceleiasi regiuni genomice. Ele au fost folosite pentru calculul matricei de distanțe si construirea arborelui filogenetic prin metoda NJ.

