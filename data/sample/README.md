# data/sample — Real sample data (small)

This folder is populated by `tools/fetch_sample_data.py` with **real** public data:

- NCBI RefSeq mRNA (FASTA): NM_000546.6 (human TP53), NM_011640.3 (mouse Trp53), NM_131327.2 (zebrafish tp53)
- UniProt protein FASTA: P04637 (human TP53)
- GEO subset CSV: small matrix for genes [TP53, MDM2, CDKN1A, BAX, BRCA1, GADD45A, ATM, CHEK2] across ~8 samples

**Why small?** Keeps CI/Codespaces fast; students will fetch larger real datasets for their PRs.

**Reproducibility:** Files include accessions; re-run this script anytime. To switch GEO, set env var `BIOINF_GSE=GSEXXXXX`.
