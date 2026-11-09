# Lab 3 - Bioinformatics: Formats & NGS

## Author
- Marius Cotoi
- GitHub: MariusC12

## Data sources

1. **FASTQ:** ENA (European Nucleotide Archive), accession DRR000288. The dataset is related to the experiment H1299_P53WT_input_091104_chip.
2. **FASTA:** TP53 human reference transcript NM_000546.6, reused from Lab 1.
3. **VCF:** ClinVar, downloaded as clinvar.vcf.gz. Variants were filtered for the TP53 gene.
4. **PubMed:** Articles retrieved using the NCBI Entrez API.

## AI usage disclosure

ChatGPT was used to help understand the assignment, develop the Python scripts, and debug errors. The scripts were executed in the GitHub Codespaces environment, and their outputs were checked by the student.

## Files

- `pubmed_query.py` - retrieves up to five PubMed articles for TP53 AND cancer.
- `pubmed_MariusC12.txt` - PubMed results.
- `fasta_MariusC12.fasta` - TP53 FASTA sequence.
- `ex04_fastq_stats.py` - calculates FASTQ quality statistics.
- `qc_report_MariusC12.txt` - FASTQ quality report.
- `vcf_pubmed.py` - extracts TP53 variants and searches PubMed.
- `variants_MariusC12.txt` - variant search results.

## How to run

```bash
python pubmed_query.py
python ex04_fastq_stats.py
python vcf_pubmed.py
```