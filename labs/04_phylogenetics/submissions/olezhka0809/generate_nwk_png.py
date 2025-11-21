from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read("tree_olezhka0809.nwk", "newick")

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, do_show=False, axes=ax)
plt.tight_layout()
plt.savefig("tree_olezhka0809.png", dpi=300)
plt.close(fig)
