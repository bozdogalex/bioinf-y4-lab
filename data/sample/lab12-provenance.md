# Lab 1 and 2 sequence provenance

Corrected on 5 October 2026. These are teaching reference records, not patient data.

| File | Record and scope | Source |
| --- | --- | --- |
| `tp53_dr_protein.fasta` | P79734, P53_DANRE, Danio rerio p53, 373 aa, sequence version 1 | [UniProt FASTA](https://rest.uniprot.org/uniprotkb/P79734.fasta) |
| `tp53_protein_multi.fasta` | Concatenation of the human P04637, mouse P02340 and zebrafish P79734 single-record files | Individual FASTA headers retain their accessions and species |
| `labs/01_intro&databases/data/brca1.gb` | NM_007294.4, human BRCA1 transcript variant 1, 7088 nt, complete sequence and annotations | [NCBI record](https://www.ncbi.nlm.nih.gov/nuccore/NM_007294.4) |
| `toy_alignment.fasta` | Artificial sequences TTTACGTAAA and GGGACGTCCC | Created for teaching; not a biological TP53 fragment |

The earlier zebrafish protein file and the corresponding combined-file entry incorrectly contained mouse Exportin-1 (Q6P5F9). Both sequence and header have been replaced. Do not compare old and corrected MSA results as if the input were unchanged.

The former BRCA1 cache was an undefined-sequence chromosome CON record (NC_060941.1), not a BRCA1 transcript. The demo now reads the bundled transcript by default and downloads exactly NM_007294.4 only with `--refresh`.

The existing nucleotide multi-FASTA contains RefSeq transcripts NM_000546.6, NM_011640.3 and NM_131327.2, represented with T. They include untranslated regions. Protein reference records were selected independently; do not assume each is the translated product of the particular transcript variant without checking CDS and isoform annotations.
