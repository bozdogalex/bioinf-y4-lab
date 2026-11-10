# Lab1

## Demo script 1
root@codespaces-fe58e9:/workspaces/bioinf-y4-lab# python labs/"01_intro&databases"/demo01_entrez_brca1.py
Source: bundled offline cache: brca1.gb
ID: NM_007294.4
Title: Homo sapiens BRCA1 DNA repair associated (BRCA1), transcript variant 1, mRNA
Length: 7088 nt
GC fraction: 0.418
First 50 nt: GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTG
This is an mRNA transcript represented using T, not the complete genomic locus.

## Demo script 2
root@codespaces-fe58e9:/workspaces/bioinf-y4-lab# python labs/"01_intro&databases"/demo02_seq_ops.py
DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]

## Demo script 3
root@codespaces-fe58e9:/workspaces/bioinf-y4-lab# python labs/"01_intro&databases"/demo03_dbsnp.py
Am găsit 5 SNP IDs.
SNP_ID: 2552282559 | CHRPOS: 17:43127349 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281972 | CHRPOS: 17:43127232 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281880 | CHRPOS: 17:43127053 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281808 | CHRPOS: 17:43126996 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281780 | CHRPOS: 17:43126982 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant

## Exercise 1: multifasta_gc
root@codespaces-fe58e9:/workspaces/bioinf-y4-lab# python labs/"01_intro&databases"/submissions/Tavi1212/ex01_multifasta_gc.py --email octaviancili@gmail.com --accession NM_000546 --out data/work/Tavi1212/lab01/nm000546.fa
[ok] Am scris 1 înregistrări în: data/work/Tavi1212/lab01/nm000546.fa
NM_000546.6     GC=0.534