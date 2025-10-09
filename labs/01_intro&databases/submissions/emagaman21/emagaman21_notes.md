#demo1
Găsite 1 rezultate.
ID: 2194972897
Titlu: Homo sapiens isolate CHM13 chromosome 17, alternate assembly T2T-CHM13v2.0
Length: 84276897 bp
GC fraction: 0.453
First 50 nt: CCTAACCCTAACCCATAACCCTAACCCTAACCTACCCTAACCCTAACCCT

#demo 2
root@codespaces-381dab:/workspaces/bioinf-y4-lab/labs/01_intro&databases# python demo02_seq_ops.py
DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]

#demo3
root@codespaces-381dab:/workspaces/bioinf-y4-lab/labs/01_intro&databases# python demo03_dbsnp.py
Am găsit 5 SNP IDs.
SNP_ID: 2552282559 | CHRPOS: 17:43127349 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281972 | CHRPOS: 17:43127232 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281880 | CHRPOS: 17:43127053 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281808 | CHRPOS: 17:43126996 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281780 | CHRPOS: 17:43126982 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant


#Exercitiu 
root@codespaces-381dab:/workspaces/bioinf-y4-lab/labs/01_intro&databases# python ex01_multifasta_gc.py --email emanuela.gaman@student.upt.ro --accession NM_000546 --out data/work/emagaman21/lab01/nm000546.fa
[ok] Am scris 1 înregistrări în: data/work/emagaman21/lab01/nm000546.fa
NM_000546.6     GC=0.534

root@codespaces-381dab:/workspaces/bioinf-y4-lab/labs/01_intro&databases# python ex01_multifasta_gc.py --email emanuela.gaman@student.upt.ro --query "TP53[Gene] AND Homo sapiens[Organism]" --retmax 3 --
out data/work/emagaman21/lab01/my_tp53.fa
[ok] Am scris 3 înregistrări în: data/work/emagaman21/lab01/my_tp53.fa
NG_017013.2     GC=0.490
NC_060941.1     GC=0.453
NC_000017.11    GC=0.453