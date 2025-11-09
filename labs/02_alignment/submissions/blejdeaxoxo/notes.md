**demo1**
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

**demo2**
pair,hamming,p_distance,len_used
NG_017013.2-NC_060941.1,24498,0.7475,32772
NG_017013.2-NC_000017.11,32772,1.0000,32772
NC_060941.1-NC_000017.11,62200462,0.7471,83257441

Alinierea globala este de preferat atunci cand vrem sa comparam doua secvente cap la cap, adica acestea au lungimi similare si ne asteptam sa fie asemanatoare pe intreaga lungime.
Alinierea locala este de preferat atunci cand avem secvente foarte lungi sau foarte diferite si vrem sa gasim bucati similare intre acestea fara a fi nevoia ca intreaga secventa sa fie aliniata.