**demo01_sntrez_brca1**
Găsite 1 rezultate.
ID: 3070263245
Titlu: Homo sapiens isolate AB17 BRCA1 protein (BRCA1) gene, partial cds
Length: 350 bp
GC fraction: 0.423
First 50 nt: CCTGATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGA

**demo02_seq_ops**
DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]


**demo03_dbsnp**
Am găsit 5 SNP IDs.
SNP_ID: 2552282559 | CHRPOS: 17:43127349 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281972 | CHRPOS: 17:43127232 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281880 | CHRPOS: 17:43127053 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281808 | CHRPOS: 17:43126996 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281780 | CHRPOS: 17:43126982 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant

**ex01_multifasta_gc**
python ex01_multifasta_gc.py \
  --email student@example.com \
  --accession NM_000546 \
  --out data/work/flori/lab01/nm000546.fa
[ok] Am scris 1 înregistrări în: data/work/flori/lab01/nm000546.fa
NM_000546.6     GC=0.534

python ex01_multifasta_gc.py \
  --email student@example.com \
  --query "TP53[Gene] AND Homo sapiens[Organism]" \
  --retmax 3 \
  --out data/work/flori/lab01/tp53.fa
[ok] Am scris 3 înregistrări în: data/work/flori/lab01/tp53.fa
NG_017013.2     GC=0.490
NC_060941.1     GC=0.453
NC_000017.11    GC=0.453