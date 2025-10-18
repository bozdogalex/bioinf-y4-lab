# demo01_entrez_brca1

root@codespaces-642e05:/workspaces/bioinf-y4-lab# /usr/local/bin/python "/workspaces/bioinf-y4-lab/labs/01_intro&databases/demo01_entrez_brca1.py"
Găsite 1 rezultate.
ID: 2194972897
Titlu: Homo sapiens isolate CHM13 chromosome 17, alternate assembly T2T-CHM13v2.0
Length: 84276897 bp
GC fraction: 0.453
First 50 nt: CCTAACCCTAACCCATAACCCTAACCCTAACCTACCCTAACCCTAACCCT



# demo02_seq_ops.py

root@codespaces-642e05:/workspaces/bioinf-y4-lab# /usr/local/bin/python "/workspaces/bioinf-y4-lab/labs/01_intro&databases/demo02_seq_ops.py"
DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]



# demo03_dbsnp.py

root@codespaces-642e05:/workspaces/bioinf-y4-lab# /usr/local/bin/python "/workspaces/bioinf-y4-lab/labs/01_intro&databases/demo03_dbsnp.py"
Am găsit 5 SNP IDs.
SNP_ID: 2552282559 | CHRPOS: 17:43127349 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281972 | CHRPOS: 17:43127232 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281880 | CHRPOS: 17:43127053 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281808 | CHRPOS: 17:43126996 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281780 | CHRPOS: 17:43126982 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant



# ex01

# query
root@codespaces-642e05:/workspaces/bioinf-y4-lab/labs/01_intro&databases#  python ex01_multifasta_gc.py --email mark.inisconi@yahoo.ro \
        --query "TP53[Gene] AND Homo sapiens[Organism]" \
        --retmax 3 \
        --out data/work/MarkLucasInisconi/lab01/my_tp53.fa
[ok] Am scris 3 înregistrări în: data/work/MarkLucasInisconi/lab01/my_tp53.fa
NG_017013.2     GC=0.490
NC_060941.1     GC=0.453
NC_000017.11    GC=0.453

# accession
root@codespaces-642e05:/workspaces/bioinf-y4-lab/labs/01_intro&databases# python ex01_multifasta_gc.py --email mark.inisconi@yahoo.ro \                                                
        --accession NM_000546 \
        --out data/work/MarkLucasInisconi/lab01/nm000546.fa
[ok] Am scris 1 înregistrări în: data/work/MarkLucasInisconi/lab01/nm000546.fa
NM_000546.6     GC=0.534

# rezultat
GC fraction = 0.534