#rulare demo 1
root@codespaces-f90576:/workspaces/bioinf-y4-lab# /usr/local/bin/python "/workspaces/bioinf-y4-lab/labs/03_formats&NGS/demo01_fastq_qc.py"
Reads: 2
Mean length: 8.50
N rate: 0.0588
Mean Phred: 40.00

#rulare demo 2
root@codespaces-f90576:/workspaces/bioinf-y4-lab# /usr/local/bin/python "/workspaces/bioinf-y4-lab/labs/03_formats&NGS/demo02_mapping_toy.py"
Read ATGCTAGC mapped at position 0
Read GATCGATC mapped at position 19
Read TACGATCG mapped at position 28
Read GGGGGGGG did not map

#rulare exercitiu
root@codespaces-f90576:/workspaces/bioinf-y4-lab# /usr/local/bin/python "labs/03_formats&NGS/ex01_fetch_fastq.py" SRR3180792
[INFO] Descărcare: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/002/SRR3180792/SRR3180792_1.fastq.gz
Downloaded: /workspaces/bioinf-y4-lab/data/work/MarioCld/lab03/SRR3180792_1.fastq.gz

root@codespaces-f90576:/workspaces/bioinf-y4-lab# /usr/local/bin/python "labs/03_formats&NGS/submissions/MarioCld/ex04_fastq_stats.py"
[OK] QC report -> /workspaces/bioinf-y4-lab/labs/03_formats&NGS/submissions/MarioCld/qc_report_MarioCld.txt

- Accession SRA: **SRR3180792**  
- Link ENA: [https://www.ebi.ac.uk/ena/browser/view/SRR3180792](https://www.ebi.ac.uk/ena/browser/view/SRR3180792)

Verificarea calității datelor (FASTQ QC) este esențială înainte de analiza variantelor din mai multe motive:

- identificăm citirile proaste sau corupte: citirile cu scoruri Phred scăzute pot induce erori false în aliniere și detecția variantelor
- controlăm lungimea citirilor: citirile foarte scurte sau incomplete pot afecta acuratețea aliniamentului și a apelului de variante
- se evită bias-ul: baze cu calitate slabă pot introduce bias în frecvența variantelor detectate
- optimizăm resurse: filtrarea citirilor de calitate slabă reduce timpul și memoria necesare pentru aliniere și analiză ulterioară