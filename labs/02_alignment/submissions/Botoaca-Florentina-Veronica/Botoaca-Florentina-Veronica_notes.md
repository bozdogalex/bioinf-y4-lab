## Date utilizate
- **Fișier FASTA**: `/workspaces/bioinf-y4-lab/data/sample/tp53_dna_multi.fasta` 
- **Secvențe comparate**: TP53 din organismele uman și șoarece
- **Tip secvențe**: ADN (gene TP53)


### Cod rulat

///demo01_pairwise.py
root@codespaces-f06464:/workspaces/bioinf-y4-lab/labs/02_alignment# python demo01_pairwise.py --fasta /workspaces/bioinf-y4-lab/data/sample/tp53_dna_multi.fasta

[INPUT]
A: CTCAAAA
B: TTTCCCC

[GLOBAL] top alignment:
CT-CAAAA
.| | ...
TTTC-CCC
  Score=-4

[LOCAL] top alignment:
2 TC
  ||
3 TC
  Score=2


///demo2_distance_matrix.py

root@codespaces-f06464:/workspaces/bioinf-y4-lab/labs/02_alignment# python demo02_distance_matrix.
py --fasta /workspaces/bioinf-y4-lab/data/sample/tp53_dna_multi.fasta
pair,hamming,p_distance,len_used
NM_000546.6-NM_011640.3,1286,0.7221,1781
NM_000546.6-NM_131327.2,1681,0.7528,2233
NM_011640.3-NM_131327.2,1293,0.7260,1781