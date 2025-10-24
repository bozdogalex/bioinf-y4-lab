## Date folosite

*Am folosit fișierul FASTA descărcat în Lab 1, cu genele TP53 (Homo sapiens).  
*Deoarece fișierul complet conținea secvențe foarte lungi și algoritmul Needleman–Wunsch necesita prea multă memorie 
(procesul se termina cu "Terminated"), am creat o copie a fișierului care conține doar primele 500 baze din fiecare 
secvență.  
*Această versiune scurtată a permis rularea și testarea corectă a algoritmilor. 


## Reflecție

*Alinierea globala (Needleman–Wunsch) este preferată atunci când cele două secvențe sunt similare pe aproape 
toată lungimea (ex: două izoforme sau gene ortoloage).  
*Alinierea locala (Smith–Waterman) este preferată atunci când doar o parte a secvențelor este comună (ex: domenii 
proteice similare în proteine diferite).

## Rezultat demo01
[INPUT]
A: CTCCTTG
B: CCTAACC

[GLOBAL] top alignment:
C-T--CCTTG
| |  ||   
CCTAACC---
  Score=-2

[LOCAL] top alignment:
3 CCT
  |||
1 CCT
  Score=3


## Rezultat demo02

pair,hamming,p_distance,len_used
NG_017013.2-NC_060941.1,24498,0.7475,32772
NG_017013.2-NC_000017.11,32772,1.0000,32772
NC_060941.1-NC_000017.11,62200462,0.7471,83257441
