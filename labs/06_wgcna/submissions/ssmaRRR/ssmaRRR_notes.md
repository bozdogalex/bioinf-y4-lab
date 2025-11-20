Am rulat demo codul , merge:
    Matrice expresie (toy):
       Sample1  Sample2  Sample3  Sample4  Sample5
        GeneA        5        4        6        5        4
        GeneB        3        3        2        4        3
        GeneC        8        9        7       10        8

        Matrice corelație (Spearman):
                GeneA     GeneB     GeneC
        GeneA  1.000000 -0.353553 -0.432590
        GeneB -0.353553  1.000000  0.917663
        GeneC -0.432590  0.917663  1.000000

        Matrice adiacență cu prag 0.7:
            GeneA  GeneB  GeneC
        GeneA      0      0      0
        GeneB      0      0      1
        GeneC      0      1      0


Am descarcat fisierul GSE113863_expression_raw_counts.csv de pe https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113863
si am generat fisierul modules_ssmaRRR.csv in labs/06_wgcna/submissions/ssmaRRR/modules_ssmaRRR.csv


Am folosit metoda de corelatie 'spearman', transformare: log2(x + 1) pe matricea de expresie, filtrare: am păstrat doar genele cu varianța > 0.5, rețeaua este neorientată, iar modulele au fost detectate cu Louvain (sau fallback pe `greedy_modularity_communities` din NetworkX). Prag de adiacență: 0.6 pe |cor| — muchie între două gene există doar dacă |cor(gene_i, gene_j)| ≥ 0.6.

Clustering-ul clasic grupează obiectele în funcție de similaritatea globală, în timp ce rețeaua de co-expresie transformă problema într-un graf și se uită la structura de conexiuni (muchii peste un prag), apoi extrage module ca comunități în acest graf.