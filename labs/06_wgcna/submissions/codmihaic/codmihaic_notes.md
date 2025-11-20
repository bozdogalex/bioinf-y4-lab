**Demo01:**
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

**Ex01:**
Pentru început, am folosit scriptul `generate_genes.py` pentru a genera o matrice cu 1.000 de gene si 20 de sample-uri, care au fost salvate in `data/work/codmihaic/lab06/expression_matrix.csv`, pentru a fi folosite la rezolvarea `ex01_gce_networks.py`.

Despre metrică, putem spune că Spearman este o corelație bazată pe ranguri, robustă la variații neliniare și outlieri, foarte folosită în rețelele de co-expresie a genelor pentru că surprinde bine trenduri comune între gene, chiar dacă nu sunt perfect liniare. De asemenea, pragul de adiacență folosit a fost de 0.6, cu USE_ABS_CORR = True, iar rețeaua rezultată include numai gene cu corelații suficient de puternice.

*Reflecție:*
În Lab 5, clustering-ul împarte probe sau gene în grupuri bazate pe asemănare directă (distanță / metrică), rezultatul fiind o singură partiționare a datelor (K-means), sau o structură arborescentă (hierarchical).

În schimb, rețeaua de co-expresie (Lab 6) construiește un graf, unde nodurile sunt genele, iar muchiile legături puternice de co-expresie (corelație mare). Modulele sunt, de fapt, comunități din graf, nu ca centre geometrice (K-means) sau distanțe (hierarchical), rețeaua punând accentul pe conectivitatea funcțională dintre gene, nu pe distanța între vectori.

Diferența cheie dintre ele este aceea că, în timp ce clustering-ul clasic grupează entitățile pe baza distanțelor din spațiul PCA, co-expresia genetică grupează genele pe baza corelațiilor și relațiilor lor în graf, ceea ce surprinde mai bine comportamentul „în rețea”.