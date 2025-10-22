### Demo01

Rezultatul rularii scriptului:
 python /workspaces/bioinf-y4-lab/labs/02_alignment/demo01_pairwise.py --fasta /workspaces/bioinf-y4-lab/data/work/filippaulescu/lab01/my_tp53.fa
/usr/local/lib/python3.11/site-packages/Bio/pairwise2.py:278: BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module.
  warnings.warn(
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

### Date utilizate pentru demo01
Am folosit fisierul FASTA cu secvente ADN:
`/workspaces/bioinf-y4-lab/data/work/filippaulescu/lab01/my_tp53.fa`
Acest fisier contine secvente extrase din gena **TP53**.

### Scurta reflectie: 
Preferam alinierea globala, bazata pe algoritmul Needleman-Wunsch, atunci cand se compara secvente inrudite, de dimensiuni similare si gaseste similaritatea pe intreaga lungime. Scopul este intelegerea evolutiei acelei gene.

Preferam alinierea locala, bazata pe algoritmul Smith-Waterman,  atunci cand comparam secvente inrudite distant sau cautam similaritati, cum ar fi domenii importante ale unei proteine. In contextul acesta se observa ca e mai potrivita alinierea locala, deoarece verifica similitudinea pe o secventa mica, igonrand capetele care nu se potrivesc. 


### Demo02
Rezultatul rularii acestui script:
python labs/02_alignment/demo02_distance_matrix.py --fasta data/sample/tp53_dna_multi.fasta
pair,hamming,p_distance,len_used
NM_000546.6-NM_011640.3,1286,0.7221,1781
NM_000546.6-NM_131327.2,1681,0.7528,2233
NM_011640.3-NM_131327.2,1293,0.7260,1781

### Date utilizate pentru demo02
Am rulat scriptul pe fisierul FASTA: `data/sample/tp53_dna_multi.fasta`.
