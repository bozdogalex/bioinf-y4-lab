# Small debugging tasks for Labs 1 and 2

The downloading and NW/SW TODOs remain the core student exercises. Complete copies in your private submission directory. These optional tasks practise debugging without relying on corrupted data or incorrect scientific definitions.

## Input normalization

The pairwise demo currently keeps the original letter case. In a private copy, compare `acgt` with `ACGT`, explain why case should not change the biological sequence, normalize input, and add a small test. Decide what your script accepts as nucleotide symbols; remember that a FASTA file alone does not guarantee a molecule type.

## Empty and ambiguous observations

`demo02_seq_ops.py` uses a known nonempty ACGT example. Adapt a copy to user input and define its behavior for an empty sequence and ambiguous symbols. Explain your GC denominator. Do not interpret amino-acid G and C as nucleotide GC content.

## Safe optional online configuration

Before running the optional `demo03_dbsnp.py`, replace its placeholder email with a value read from `NCBI_EMAIL`, validate that it exists, and close both Entrez handles with context managers. Document how you configured the environment without committing API keys. Live dbSNP is an extension, not a prerequisite for the offline practical.

## Consistent algorithm comparisons

The NW and SW skeletons deliberately retain different default scoring values. Inspect them, configure the same values when validating against Biopython, and report the complete scoring scheme. Equal optimal scores can have different traceback alignments.

## Library maintenance

The current pairwise demo uses deprecated `pairwise2`, which still works in the pinned course dependency range. Migrate a private copy to `Bio.Align.PairwiseAligner`, set every score explicitly, and compare a normal input and a no-positive-local-alignment input. This is optional; warnings alone are not failed executions.

No points are silently added to the published assignment by this list. The instructor will identify any challenge used for assessment.
