"""
Exercițiul 5 — Construirea unui arbore Neighbor-Joining

Instrucțiuni (de urmat în laborator):
1. Refolosiți secvențele din laboratoarele anterioare (FASTA din Lab 2 sau FASTQ→FASTA din Lab 3).
2. Dacă aveți doar fișiere FASTA cu o singură secvență, combinați cel puțin 3 într-un fișier multi-FASTA:
3. Salvați fișierul multi-FASTA în: data/work/<handle>/lab04/your_sequences.fasta
4. Completați pașii de mai jos:
   - încărcați multi-FASTA-ul,
   - calculați matricea de distanțe,
   - construiți arborele NJ,
   - salvați rezultatul în format Newick (.nwk).
            """
            
from pathlib import Path
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            
if __name__ == "__main__":
    # === CONFIG — editează aici după nevoi ===
    HANDLE = "rauldavid35"  # GitHub handle pentru livrabile
    FASTA_PATH = Path("data/work/rauldavid35/lab01/my_other_tp53.aln.fa")  # multi-FASTA (ideal deja aliniat)
    PRINT_ASCII = True  # afișează arborele în consolă

    # === Verificări inițiale ===
    if not FASTA_PATH.exists():
        raise FileNotFoundError(
            f"Nu găsesc fișierul: {FASTA_PATH}\n"
            "Pune un multi-FASTA (≥3 secvențe) la calea de mai sus."
        )

    records = list(SeqIO.parse(str(FASTA_PATH), "fasta"))
    if len(records) < 3:
        raise ValueError(
            f"Fișierul are doar {len(records)} secvențe (<3). "
            "Combină cel puțin trei într-un singur multi-FASTA."
        )

    # === Verificăm dacă e deja un MSA (toate lungimi egale) ===
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        Ls = ", ".join(map(str, sorted(lengths)))
        raise ValueError(
            "Inputul nu pare aliniat (secvențe cu lungimi diferite).\n"
            f"Lungimi detectate: {Ls}\n"
            "Te rog aliniază în prealabil (Clustal Omega/MAFFT) și rulează din nou."
        )

    # === Pregătim directoarele și fișierele de livrare ===
    submit_dir = Path("labs/04_phylogenetics/submissions") / HANDLE
    submit_dir.mkdir(parents=True, exist_ok=True)
    aligned_fa = submit_dir / f"aln_{HANDLE}.fasta"    # vom scrie aici un FASTA „aliniat”
    newick_path = submit_dir / f"tree_{HANDLE}.nwk"
    dist_txt    = submit_dir / f"dist_{HANDLE}.txt"

    # === Scriem în mod explicit aliniamentul ca FASTA, apoi îl citim ca Alignment ===
    # (nu folosim Bio.Align.MultipleSeqAlignment pentru a rămâne în lista de importuri)
    SeqIO.write(records, aligned_fa, "fasta")
    alignment = AlignIO.read(str(aligned_fa), "fasta")
    ids = [rec.id for rec in alignment]

    # === Model de distanță (DNA vs protein) ===
    dna_chars = set("ACGTURYKMSWBDHVNacgturykmswbdhvn-")
    is_protein = any(any(ch not in dna_chars for ch in str(rec.seq)) for rec in alignment)
    model = "blosum62" if is_protein else "identity"

    # === Matricea de distanțe ===
    calculator = DistanceCalculator(model)
    dist_matrix = calculator.get_distance(alignment)

    # — salvăm un dump text al matricei (util pentru PR)
    with open(dist_txt, "w", encoding="utf-8") as f:
        f.write("IDs:\n")
        f.write("\t".join(ids) + "\n\n")
        f.write(f"Model distanță: {model}\n\n")
        f.write("Distance matrix (Biopython repr):\n")
        f.write(str(dist_matrix) + "\n")

    # === Arborele NJ ===
    constructor = DistanceTreeConstructor(calculator, "nj")
    nj_tree = constructor.build_tree(alignment)

    # === Scriem arborele în Newick ===
    Phylo.write(nj_tree, str(newick_path), "newick")

    # === Afișare ASCII (opțional) ===
    if PRINT_ASCII:
        try:
            Phylo.draw_ascii(nj_tree)
        except Exception as e:
            print("[WARN] ASCII tree display failed:", e)

    print(f"[OK] Matrice de distanțe: {dist_txt}")
    print(f"[OK] Arbore NJ (Newick): {newick_path}")