from Bio import AlignIO

path = "labs/02_alignment/submissions/MariusC12/task3_msa.fasta"
alignment = AlignIO.read(path, "fasta")

sequences = [str(record.seq) for record in alignment]
length = alignment.get_alignment_length()

print("Numar secvente:", len(sequences))
print("Lungime aliniere:", length)

best_start = 0
best_score = -1
window = 60

for start in range(length - window + 1):
    score = 0

    for col in range(start, start + window):
        bases = [seq[col] for seq in sequences if seq[col] != "-"]

        if len(bases) >= 2:
            if len(set(bases)) < len(bases):
                score += 1

    if score > best_score:
        best_score = score
        best_start = start

print("\nFragment conservat:")
print("Coloane:", best_start, "-", best_start + window)

for record, seq in zip(alignment, sequences):
    print(">", record.id)
    print(seq[best_start:best_start + window])

print("\nColoane cu aceeasi nucleotida la toate 3:")
conserved = []

for col in range(best_start, best_start + window):
    bases = [seq[col] for seq in sequences]

    if len(set(bases)) == 1 and bases[0] != "-":
        conserved.append(col)

print(conserved)
print("Total:", len(conserved))