Date folosite: subset cu 2 secvente extras din data/work/rauldavid35/lab01/my_tp53.fa (generat cu un script Python + Biopython)
Am fost fortat sa fac asta deoarece fisieru original are 2.4 milioane de linii
Comenzi rulate:

python labs/02_alignment/submissions/rauldavid35/ex01_global_nw.py \
  --fasta data/work/rauldavid35/lab01/my_tp53_subset.fa --i1 0 --i2 1

python labs/02_alignment/submissions/rauldavid35/ex02_local_sw.py \
  --fasta data/work/rauldavid35/lab01/my_tp53_subset.fa --i1 0 --i2 1

Am avut secventele : NG_017013.2 vs NC_060941.1

Aliniere globala (NW) : SCOR -515, scor negativ, multe gap-uri si nepotriviri
Aliniere locala (SW) : SCOR 299, scor pozitiv, identifica regiuni locale cu similaritate mare

Obs. : se foloseste gloibal atunci cand secventele sunt de lungiimi comparabile, cu obiectiv cap-la-cap, iar local atunci cand lungimile difera mult, deoarece opreste si porneste alinierea unde are sens.