from pathlib import Path
import shutil
import sys
from typing import List
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
    DistanceMatrix,
)

HANDLE = "RaresCalotescu"
FASTA_IN = Path(f"data/work/{HANDLE}/lab04/your_sequences.fasta")
ALIGN_OUT = Path(f"data/work/{HANDLE}/lab04/your_sequences.aln.fasta")
TREE_OUT = Path(f"labs/04_phylogenetics/submissions/{HANDLE}/tree_{HANDLE}.nwk")
TREE_OUT.parent.mkdir(parents=True, exist_ok=True)


def all_equal_length(records: List) -> bool:
    lengths = {len(r.seq) for r in records}
    return len(lengths) == 1


def build_dm_hamming(records: List) -> DistanceMatrix:
    """p-distance (Hamming) pe secvențe deja de aceeași lungime."""
    names = [r.id[:25] for r in records]  # scurtăm ID-urile dacă sunt lungi
    n = len(records)
    seqs = [str(r.seq).upper() for r in records]
    mat = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                s1, s2 = seqs[i], seqs[j]
                L = len(s1)
                # mismatches (ignorăm poziții cu N la ambele? aici le tratăm normal)
                mism = sum(1 for a, b in zip(s1, s2) if a != b)
                row.append(mism / L)
        mat.append(row)
    return DistanceMatrix(names, mat)


def try_align_with_clustalo(in_fa: Path, out_aln: Path) -> bool:
    """Încearcă Clustal Omega (clustalo). Returnează True dacă a reușit."""
    if shutil.which("clustalo") is None:
        return False
    import subprocess

    cmd = [
        "clustalo",
        "-i", str(in_fa),
        "-o", str(out_aln),
        "--force",
        "--auto",
        "--outfmt", "fasta",
    ]
    print("[INFO] Running:", " ".join(cmd))
    subprocess.check_call(cmd)
    return out_aln.exists() and out_aln.stat().st_size > 0


def main():
    if not FASTA_IN.exists():
        sys.exit(f"[ERR] Nu găsesc {FASTA_IN}. Creează multi-FASTA-ul tău (minim 3 secvențe).")

    records = list(SeqIO.parse(str(FASTA_IN), "fasta"))
    if len(records) < 3:
        sys.exit("[ERR] Trebuie cel puțin 3 secvențe în multi-FASTA.")

    # 1) Distanțe
    if all_equal_length(records):
        print("[INFO] Toate secvențele au aceeași lungime -> folosesc p-distance (Hamming).")
        dm = build_dm_hamming(records)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
    else:
        print("[INFO] Lungimi diferite -> încerc aliniere rapidă cu Clustal Omega.")
        aligned_ok = try_align_with_clustalo(FASTA_IN, ALIGN_OUT)
        if not aligned_ok:
            sys.exit(
                "[ERR] Nu am putut alinia (clustalo indisponibil). "
                "Instalează clustalo sau furnizează un FASTA deja aliniat (cu gapuri) "
                f"în {ALIGN_OUT} și rulează din nou."
            )
        aln = AlignIO.read(str(ALIGN_OUT), "fasta")
        calc = DistanceCalculator("identity")  # dist = 1 - identity
        dm = calc.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

    # 2) Salvează arborele
    Phylo.write(tree, str(TREE_OUT), "newick")
    print(f"[OK] Arbore NJ salvat -> {TREE_OUT}")

    # 3) Print ASCII (opțional, doar preview în terminal)
    try:
        print("\n=== NJ Tree (ASCII preview) ===")
        Phylo.draw_ascii(tree)
    except Exception:
        pass


if __name__ == "__main__":
    main()
