** Demo1 rulat**
Sequences: ['NM_000546.6', 'NM_011640.3', 'NM_131327.2']
Distance matrix:
 [[0.         0.51194268 0.6691879 ]
 [0.51194268 0.         0.72599663]
 [0.6691879  0.72599663 0.        ]]

""Ex05""
Calculăm matricea de distanțe...
Matricea de distanțe:
SRR390728.1 0.000000
SRR390728.2 0.777778    0.000000
SRR390728.3 0.638889    0.694444    0.000000
SRR390728.4 0.833333    0.722222    0.861111    0.000000
SRR390728.5 0.833333    0.805556    0.833333    0.777778    0.000000
SRR390728.6 0.750000    0.722222    0.888889    0.750000    0.666667    0.000000
SRR390728.7 0.694444    0.777778    0.666667    0.722222    0.833333    0.861111    0.000000
SRR390728.8 0.583333    0.805556    0.722222    0.750000    0.722222    0.611111    0.750000    0.000000
SRR390728.9 0.777778    0.722222    0.777778    0.888889    0.666667    0.638889    0.972222    0.777778    0.000000
SRR390728.10    0.777778    0.888889    0.722222    0.805556    0.805556    0.722222    0.694444    0.833333    0.916667    0.000000
    SRR390728.1 SRR390728.2 SRR390728.3 SRR390728.4 SRR390728.5 SRR390728.6 SRR390728.7 SRR390728.8 SRR390728.9 SRR390728.10

Construim arborele Neighbor-Joining...
Arbore construit!
Număr de clade/ramuri: 10

Arbore salvat în: labs/04_phylogenetics/submissions/Ana-Maria-Bojan/tree_Ana-Maria-Bojan.nwk

Vizualizare arbore (text):
               _________________________________ SRR390728.8
         _____|
        |     |_________________________________ SRR390728.1
        |
    ____|            __________________________________ SRR390728.6
   |    |    _______|
   |    |   |       |   ________________________________________ SRR390728.9
   |    |   |       |__|
   |    |___|          |_____________________________________ SRR390728.5
  _|        |
 | |        |  __________________________________________ SRR390728.4
 | |        |_|
 | |          |________________________________________ SRR390728.2
_| |
 | |_____________________________________ SRR390728.3
 |
 |____________________________________ SRR390728.7
 |
 |__________________________________________ SRR390728.10


Structura arborelui:
  - Inner8: branch_length = 0
  - Inner7: branch_length = 0.02213541666666674
  - Inner6: branch_length = 0.04123263888888884
  - Inner3: branch_length = 0.048611111111111216
  - SRR390728.8: branch_length = 0.29513888888888884
  - SRR390728.1: branch_length = 0.2881944444444444
  - Inner5: branch_length = 0.036458333333333315
  - Inner2: branch_length = 0.06770833333333329
  - SRR390728.6: branch_length = 0.29662698412698413
  - Inner1: branch_length = 0.022817460317460292
  - SRR390728.9: branch_length = 0.34548611111111105
  - SRR390728.5: branch_length = 0.3211805555555557
  - Inner4: branch_length = 0.012152777777777818
  - SRR390728.4: branch_length = 0.36875
  - SRR390728.2: branch_length = 0.3534722222222222
  - SRR390728.3: branch_length = 0.3250868055555556
  - SRR390728.7: branch_length = 0.32248263888888884
  - SRR390728.10: branch_length = 0.37196180555555547

**Ce informații suplimentare oferă arborele filogenetic față de o simplă matrice de distanțe?**

Arborele filogenetic transformă datele cantitative dintr-o matrice de distanțe într-o poveste evolutivă, arătând nu doar cât de asemănători sunt organismii, ci și cum sunt înrudiți ierarhic. El revelează relațiile de înrudire, grupele monofiletice și ordinea în care speciile s-au separat de la strămoșii comuni. Spre deosebire de matricea de distanțe care oferă doar similarități perechi, arborele filogenetic reconstituie istoria evolutivă și direcția schimbărilor, oferind o perspectivă dinamică asupra proceselor de speciație.