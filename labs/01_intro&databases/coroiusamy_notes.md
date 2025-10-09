
python demo01_entrez_brca1.py 
Găsite 1 rezultate.
ID: 2194972897
Titlu: Homo sapiens isolate CHM13 chromosome 17, alternate assembly T2T-CHM13v2.0
Length: 84276897 bp
GC fraction: 0.453
First 50 nt: CCTAACCCTAACCCATAACCCTAACCCTAACCTACCCTAACCCTAACCCT

python demo02_seq_ops.py 
DNA: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
RNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Protein: MAIVMGR*KGAR*
Reverse complement: CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
GC fraction: 0.564
ATG positions: [0, 12]

python demo03_dbsnp.py 
Am găsit 5 SNP IDs.
Traceback (most recent call last):
  File "/workspaces/bioinf-y4-lab/labs/01_intro&databases/demo03_dbsnp.py", line 17, in <module>
    snp_id = d.get("SNP_ID")
             ^^^^^
AttributeError: 'str' object has no attribute 'get'

NG_017013.2     GC=0.490