
---

demo01:
Reads: 4437
Mean length: 199.16
N rate: 0.0000
Mean Phred: 13.34

---

demo02:
Read ATGCTAGC mapped at position 0
Read GATCGATC mapped at position 19
Read TACGATCG mapped at position 28
Read GGGGGGGG did not map

date rulate pentru un fisier fastq de 2.4 mb

---

ex01_fetch_fastq_qc.py

root@codespaces-f06464:/workspaces/bioinf-y4-lab/labs/03_formats&NGS/submissions/Botoaca-Florentina-Veronica# python ex01_fetch_fastq.py ERR022075
🔍 Interoghez ENA API pentru ERR022075...
📡 Status code: 200
📡 URL: https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR022075&result=read_run&fields=fastq_ftp%2Csample_accession%2Cinstrument_platform&format=json&download=false
📊 Sample: SAMEA779983, Platform: ILLUMINA
🔗 Am găsit 2 fișier(e) FASTQ
   - ERR022075_1.fastq.gz
   - ERR022075_2.fastq.gz
⬇️  Descarc ERR022075_1.fastq.gz...
   📦 Progres: 100.0%
✅ Descărcat cu succes: /workspaces/bioinf-y4-lab/labs/01_intro&databases/data/work/Botoaca-Florentina-Veronica/lab01/your_reads.fastq.gz
📊 Dimensiune: 1788.81 MB (1875702420 bytes)

---
