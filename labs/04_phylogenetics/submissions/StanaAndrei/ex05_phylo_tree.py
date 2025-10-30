#!/usr/bin/env python3
"""
Construire arbore Neighbor-Joining dintr-un multi-FASTA (aliniat).
Setează variabila `handle` mai jos înainte de rulare.
"""
from pathlib import Path
import sys

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def main(handle: str):
    # Paths
    fasta = Path(f"data/sample/tp53_dr_protein.fasta")
    output = Path(f"labs/04_phylogenetics/submissions/{handle}/tree_{handle}.nwk")

    if not fasta.exists():
        print(f"[ERROR] Fișierul FASTA nu există: {fasta}", file=sys.stderr)
        return 1

    # Încarcă multi-FASTA (se așteaptă un alignment; toate secvențele trebuie să aibă aceeași lungime)
    try:
        alignment = AlignIO.read(str(fasta), "fasta")
    except Exception as e:
        print(f"[ERROR] Nu am putut citi fișierul ca alignment FASTA: {e}", file=sys.stderr)
        print("  Asigură-te că fișierul conține secvențe aliniate (same length).", file=sys.stderr)
        return 1

    # Verifică dacă toate secvențele au aceeași lungime
    lengths = [len(rec.seq) for rec in alignment]
    if len(set(lengths)) != 1:
        print("[ERROR] Secvențele din FASTA nu au aceeași lungime. Trebuie să furnizezi un FASTA aliniat.", file=sys.stderr)
        print("  Recomandare: aliniere cu MAFFT/MUSCLE/Clustal înainte de a rula acest script.", file=sys.stderr)
        return 1

    # Calculează matricea de distanțe (folosesc 'identity' - proporția de site-uri diferite)
    calculator = DistanceCalculator("identity")
    try:
        dm = calculator.get_distance(alignment)
    except Exception as e:
        print(f"[ERROR] Eroare la calculul matricei de distanțe: {e}", file=sys.stderr)
        return 1

    # Construiește arborele Neighbor-Joining
    constructor = DistanceTreeConstructor()
    try:
        tree = constructor.nj(dm)
    except Exception as e:
        print(f"[ERROR] Eroare la construcția arborelui NJ: {e}", file=sys.stderr)
        return 1

    # Salvează arborele în format Newick (.nwk)
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        Phylo.write(tree, str(output), "newick")
        print(f"[OK] Arbore salvat: {output}")
    except Exception as e:
        print(f"[ERROR] Nu am putut scrie fișierul Newick: {e}", file=sys.stderr)
        return 1

    # Afișare ASCII în consolă
    try:
        print("\n===== Arbore (ASCII) =====")
        Phylo.draw_ascii(tree)
        print("==========================\n")
    except Exception:
        # draw_ascii poate eșua în unele medii; nu este critic
        print("[INFO] Afișare ASCII eșuată (mediul curent poate să nu o suporte).")

    return 0

if __name__ == "__main__":
    handle = "StanaAndrei"

    rc = main(handle)
    sys.exit(rc)
