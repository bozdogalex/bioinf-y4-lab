## Rezultate observate

### Cod rulat

root@codespaces-74a762:/workspaces/bioinf-y4-lab/labs/02_alignment# python demo01_pairwise.py --fasta /workspaces/bioinf-y4-lab/data/work/samycoroiu/lab01/my_tp53.fa 
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

root@codespaces-74a762:/workspaces/bioinf-y4-lab/labs/02_alignment# python demo02_distance_matrix.py --fasta /workspaces/bioinf-y4-lab/data/work/samycoroiu/lab0
1/my_tp53.fa 
pair,hamming,p_distance,len_used
NG_017013.2-NC_060941.1,24498,0.7475,32772
NG_017013.2-NC_000017.11,32772,1.0000,32772
NC_060941.1-NC_000017.11,62200462,0.7471,83257441
root@codespaces-74a762:/workspaces/bioinf-y4-lab/labs/02_alignment#

### Local alignment
root@codespaces-74a762:/workspaces/bioinf-y4-lab# python labs/02_alignment/submissions/coroiusamy/ex02_local_sw.py --fasta data/work/samycoroiu/lab01/my_tp53.fa --i1 0 --i2 1 --limit 250
usage: ex02_local_sw.py [-h] --fasta FASTA [--i1 I1] [--i2 I2]
ex02_local_sw.py: error: unrecognized arguments: --limit 250
root@codespaces-74a762:/workspaces/bioinf-y4-lab# python labs/02_alignment/submissions/coroiusamy/ex02_local_sw.py --fasta data/work/samycoroiu/lab01/my_tp53.fa --i1 0 --i2 1 --limit 250
[Info] Secvențele au fost limitate la primele 250 baze.
======================================================================
Aliniere locală Smith-Waterman
======================================================================
Secvență 1: NG_017013.2 (lungime: 250 bp)
Secvență 2: NC_060941.1 (lungime: 250 bp)
Se calculează alinierea...

=== Rezultat Aliniere ===
NG_017013.2 vs NC_060941.1
CTCGAACTCCTGACCTCAGGTGATCC--A-CCT-GCCTCAGCCT--CCCAAAGTGCTGGGATTA--C-AGGA-GTCAGCC--ACCGC--ACCC-AGCCCCAA--CTAA
CT--AAC-CCTAACC-CA--TAACCCTAACCCTAACCT-ACCCTAACCCTAA-CCCT---A--ACCCTA--ACCTAACCCTAACC-CTAACCCTAACCCTAACCCTAA
Score: 81

### Global alignment
root@codespaces-74a762:/workspaces/bioinf-y4-lab# python labs/02_alignment/submissions/coroiusamy/ex01_global_nw.py --fasta data/work/samycoroiu/lab01/my_tp53.f
a --i1 0 --i2 1 --limit 250
[Info] Secvențele au fost limitate la primele 250 baze.

======================================================================
Aliniere globală Needleman-Wunsch
======================================================================
Secvență 1: NG_017013.2 (lungime: 250 bp)
Secvență 2: NC_060941.1 (lungime: 250 bp)
Se calculează alinierea...

======================================================================
=== Rezultat Aliniere Globală (NW) ===
======================================================================
Scor final: -76
Lungime aliniere: 262

Statistici:
  Matches: 105
  Mismatches: 133
  Gaps: 24
  Identitate: 40.08%

Primele 200 caractere din aliniere:
NG_017013. CTCCTTGGTTCAAGTAATTCTCCTGCCTCAGACTCCAGAGTAGCTGGGATTACAGGCGCCCGCC-ACCACGCCCAGCTAA--TTTTTTGTATTTTTAATAGAGATGGGGTTTC-ATCATGTTGGCC-AGGCTGGTCTCGAACTCCTGACCTCAGGTGATCCACCTGCCT-CAGCCT--CCCAAAGTGCT-GGGATTACAG
             |||      ||   ||   |||  | |  || | |   || |    |   |   | |   || | | |   |  ||||   |     ||    |||     |        | | |   |   || |  |    | | ||| ||| |||| |    | ||     |||  | |||  ||| ||   ||     | ||  
NC_060941. --CCTAACCCTAACCCATAACCCTAACCCTAAC-CTACCCTAAC-CCTAACCCTAACCCTAACCTAACCCTAAC-CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC-CCTAACCCTAACCCTAAC-CCAAC-CCTAACCTAACCCTAACCCTAACCCTAAACCCTAACCCTAA-CCCTAACCCTAACCC
======================================================================

root@codespaces-74a762:/workspaces/bioinf-y4-lab# 

### Context secvențe
Cele 3 secvențe TP53 din fișier:
- **NG_017013.2**
- **NC_060941.1**
- **NC_000017.11**

### Observații din demo02
Distance matrix arată:
- p-distance de 0.7475 între NG_017013.2 și NC_060941.1
- p-distance de 1.0000 între NG_017013.2 și NC_000017.11

**Interpretare:** Distanțele mari (>0.7) indică fie:
1. Compararea unor secvențe de lungimi foarte diferite (genomic vs coding)
2. Regiuni foarte divergente
3. Secvențe din contexte genomice diferite

Acest lucru face alinierea locală mai potrivită pentru găsirea  regiunilor conservate