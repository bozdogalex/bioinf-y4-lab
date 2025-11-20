"""
Exercițiul 5 — Construirea unui arbore Neighbor-Joining

Instrucțiuni (de urmat în laborator):
1. Refolosiți secvențele din laboratoarele anterioare (FASTA din Lab 2 sau FASTQ→FASTA din Lab 3).
2. Dacă aveți doar fișiere FASTA cu o singură secvență, combinați cel puțin 3 într-un fișier multi-FASTA:
3. Salvați fișierul multi-FASTA în: data/work/codmihaic/lab04/your_sequences.fasta
4. Completați pașii de mai jos:
   - încărcați multi-FASTA-ul,
   - calculați matricea de distanțe,
   - construiți arborele NJ,
   - salvați rezultatul în format Newick (.nwk).
"""

import sys, matplotlib.pyplot as plt
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def detect_seq_type(msa):
    """Întoarce 'identity' (ADN/ARN) sau 'blosum62' (proteine) pentru DistanceCalculator."""
    dna_letters = set("ACGTUN-")
    all_letters = set("".join(str(rec.seq).upper() for rec in msa))
    # dacă apar litere specifice proteinelor, folosim BLOSUM62
    proteinish = any(ch.isalpha() and ch not in dna_letters for ch in all_letters)
    return "blosum62" if proteinish else "identity"

if __name__ == "__main__":
    # TODO 1: Încărcați fișierul multi-FASTA propriu
   fasta = Path("data/work/codmihaic/lab04/codmihaic_seq.fasta")
   if not fasta.exists():
      sys.exit(f"[EROARE] Nu găsesc {fasta}. Creează un multi-FASTA cu ≥3 secvențe.")

   alignment = AlignIO.read(str(fasta), "fasta")
   lengths = {len(rec.seq) for rec in alignment}
   if len(lengths) != 1:
      print("[AVERTISMENT] Secvențele nu au aceeași lungime; pare că FASTA-ul nu e aliniat.")
      print("            Vă rog furnizați un multi-FASTA ALINIAT (cu '-').")
   print(f"[INFO] {len(alignment)} secvențe încărcate din {fasta.name}")

    # TODO 2: Calculați matricea de distanțe
   model = detect_seq_type(alignment)     # 'identity' pentru ADN/ARN, 'blosum62' pentru proteine
   calculator = DistanceCalculator(model)
   dm = calculator.get_distance(alignment)
   print("[INFO] Matricea de distanțe:")
   print(dm)

    # TODO 3: Construiți arborele NJ
   constructor = DistanceTreeConstructor()
   nj_tree = constructor.nj(dm)

    # TODO 4: Salvați arborele în format Newick
   out_dir = Path(f"labs/04_phylogenetics/submissions/codmihaic")
   out_dir.mkdir(parents=True, exist_ok=True)
   out_newick = out_dir / f"tree_codmihaic.nwk"
   Phylo.write(nj_tree, str(out_newick), "newick")
   print(f"[OK] Arbore salvat în Newick: {out_newick}")

    # TODO 5 Vizualizați arborele
   print("\n[ASCII tree]")
   Phylo.draw_ascii(nj_tree)

   out_png = out_dir / f"tree_codmihaic.png"
   fig = plt.figure(figsize=(6, 6), dpi=150)
   ax = fig.add_subplot(1, 1, 1)
   Phylo.draw(nj_tree, do_show=False, axes=ax)
   plt.tight_layout()
   fig.savefig(out_png)
   print(f"[OK] PNG salvat: {out_png}")
