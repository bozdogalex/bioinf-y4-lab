demo01_entrez_brca1.py returned:
Găsite 1 rezultate.
ID: 3070263245
Titlu: Homo sapiens isolate AB17 BRCA1 protein (BRCA1) gene, partial cds
Length: 350 bp
GC fraction: 0.423
First 50 nt: CCTGATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGA

demo02_seq_ops.py returned:
DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]

demo03_dbsnp.py returned:
SNP_ID: 2552282559 | CHRPOS: 17:43127349 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281972 | CHRPOS: 17:43127232 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281880 | CHRPOS: 17:43127053 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281808 | CHRPOS: 17:43126996 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant
SNP_ID: 2552281780 | CHRPOS: 17:43126982 | Funcție: upstream_transcript_variant,2KB_upstream_variant,intron_variant

python ex01_multifasta_gc.py --email marius.jalba@student.upt.ro --query "TP53[Gene] AND Homo sapiens[Organism]" --retmax 3 --out data/work/MariusJalba/lab01/my_tp53.fa

Caut 'TP53[Gene] AND Homo sapiens[Organism]'
Am gasit 3 inregistrari, descarc FASTA
Au fost adaugate 3 inregistrari in data/work/MariusJalba/lab01/my_tp53.fa
NG_017013.2     GC=0.490
NC_060941.1     GC=0.453
NC_000017.11    GC=0.453

python ex01_multifasta_gc.py --email marius.jalba@student.upt.ro --accession NM_000546 --out data/work/MariusJalba/lab01/nm000546.fa
Am gasit 1 inregistrari, descarc FASTA
Au fost adaugate 1 inregistrari in data/work/MariusJalba/lab01/nm000546.fa
NM_000546.6     GC=0.534