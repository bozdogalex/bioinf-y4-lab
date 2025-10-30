**Demo01:**
Traceback (most recent call last):
  File "/workspaces/bioinf-y4-lab/labs/03_formats&NGS/demo01_fastq_qc.py", line 23, in <module>
    for record in SeqIO.parse(fastq_file, "fastq"):
  File "/usr/local/lib/python3.11/site-packages/Bio/SeqIO/QualityIO.py", line 1087, in __next__
    raise ValueError("Unexpected end of file")
ValueError: Unexpected end of file

Din ce am văzut, îmi apare că fișierul sample.fastq are o problemă, generând un "end of file" neașteptat. Am văzut că gena e trunchiată, așa că am înlocuit-o pe ultima cu una free de pe net. Rezultatul este acum:

Reads: 3
Mean length: 196.00
N rate: 0.0000
Mean Phred: 14.29

**Demo02:**
Read ATGCTAGC mapped at position 0
Read GATCGATC mapped at position 19
Read TACGATCG mapped at position 28
Read GGGGGGGG did not map

**Ex01:**
Downloading from: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/SRR390728/SRR390728_1.fastq.gz
Downloaded: /workspaces/bioinf-y4-lab/data/work/codmihaic/lab03/your_reads.fastq.gz

**Ex02:**
[OK] QC report -> /workspaces/bioinf-y4-lab/labs/03_formats&NGS/submissions/codmihaic/qc_report_codmihaic.txt