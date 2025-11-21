## demo01_entrez_brca1 result

Găsite 1 rezultate.
ID: 3070263245
Titlu: Homo sapiens isolate AB17 BRCA1 protein (BRCA1) gene, partial cds
Length: 350 bp
GC fraction: 0.423
First 50 nt: CCTGATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGA

## demo02_seq_ops result

DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]

## demo03_dbsnp result

Am găsit 5 SNP IDs.
SNP_ID: 2552282559 | CHRPOS: 17:43127349 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281972 | CHRPOS: 17:43127232 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281880 | CHRPOS: 17:43127053 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281808 | CHRPOS: 17:43126996 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281780 | CHRPOS: 17:43126982 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant

## Command used

python ex01_multifasta_gc.py \
  --email oleg.garnautan@student.upt.ro \
  --query "TP53[Gene] AND Homo sapiens[Organism]" \
  --retmax 3 \
  --out data/work/olezhka0809/lab01/my_tp53_dr_protein.fasta

## Results obtained

[ok] Am scris 3 înregistrări în: data/work/olezhka0809/lab01/my_tp53_dr_protein.fasta
NG_017013.2     GC=0.490
NC_060941.1     GC=0.453
NC_000017.11    GC=0.453