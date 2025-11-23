
ex: demo01_fastq_qc.py
fisier folosit: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR684/SRR684066/SRR684066.fastq.gz
output:
root@codespaces-ad53c2:/workspaces/bioinf-y4-lab# python 'labs/03_formats&NGS/demo01_fastq_qc.py'
Reads: 5483467
Mean length: 46.00
N rate: 0.0002
Mean Phred: 31.37

Verificarea calității datelor NGS este crucială pentru că erorile de secvențiere, citirile scurte sau contaminările pot produce rezultate false la aliniere și detecția variantelor. QC-ul previne interpretări biologice greșite.

ex: demo02_mapping_toy.py
Read ATGCTAGC mapped at position 0
Read GATCGATC mapped at position 19
Read TACGATCG mapped at position 28
Read GGGGGGGG did not map

