import sys, matplotlib.pyplot as plt
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
	# Configurare specifică utilizatorului
	HANDLE = "filippaulescu"
	FASTA_PATH = Path(f"data/work/{HANDLE}/lab04/my_sequences.fasta")
	OUT_DIR = Path(f"labs/04_phylogenetics/submissions/{HANDLE}")
	OUT_NEWICK = OUT_DIR / f"tree_{HANDLE}.nwk"

	# --- Verificare și încărcare ---
	if not FASTA_PATH.exists():
		# Oprim execuția dacă fișierul de intrare nu există
		sys.exit(f"[EROARE] Nu găsesc {FASTA_PATH}. Creează un multi-FASTA ALINIAT cu ≥3 secvențe.")

	# TODO 1: Încărcați fișierul multi-FASTA propriu
	alignment = AlignIO.read(str(FASTA_PATH), "fasta")
	print(f"[INFO] {len(alignment)} secvențe încărcate din {FASTA_PATH.name}")

	# TODO 2: Calculați matricea de distanțe
	# Folosim 'identity' (distanța p sau distanța Hamming) ca model de bază.
	calculator = DistanceCalculator("identity")
	dm = calculator.get_distance(alignment)
	print("[INFO] Matricea de distanțe:")
	print(dm)

	# TODO 3: Construiți arborele NJ
	constructor = DistanceTreeConstructor()
	nj_tree = constructor.nj(dm)

	# TODO 4: Salvați arborele în format Newick
	OUT_DIR.mkdir(parents=True, exist_ok=True)
	Phylo.write(nj_tree, str(OUT_NEWICK), "newick")
	print(f"[OK] Arbore salvat în Newick: {OUT_NEWICK}")

	# TODO 5 Vizualizați arborele
	print("\n[ASCII tree]")
	Phylo.draw_ascii(nj_tree)

	# --- Generarea și salvarea imaginii PNG ---
	out_png = OUT_DIR / f"tree_{HANDLE}.png"
	# Setăm dimensiunea figurii pentru o vizualizare mai clară
	fig = plt.figure(figsize=(6, 6), dpi=150)
	ax = fig.add_subplot(1, 1, 1)
	# Desenăm arborele
	Phylo.draw(nj_tree, do_show=False, axes=ax)
	plt.tight_layout()
	fig.savefig(out_png)
	print(f"[OK] PNG salvat: {out_png}")
