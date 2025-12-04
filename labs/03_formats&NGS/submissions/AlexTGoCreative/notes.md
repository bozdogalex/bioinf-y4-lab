# Lab 03 - Notes

**Autor:** AlexTGoCreative
**Data:** 21 Noiembrie 2025

---

## De ce este important QC înainte de variant calling?

Quality Control (QC) înainte de variant calling este esențial pentru a identifica și filtra citirilе de slabă calitate care pot genera apeluri de variante false pozitive. Scorurile Phred scăzute (<Q20) indică o probabilitate mare de eroare în apelarea bazelor, ceea ce poate duce la identificarea incorectă a variantelor genetice. De asemenea, o proporție ridicată de baze N sau contaminări pot distorsiona rezultatele analizei downstream, afectând interpretarea biologică și clinică a datelor. Prin aplicarea unui QC riguros, putem elimina datele problematice înainte ca acestea să influențeze negative pipeline-ul de analiză, crescând astfel acuratețea și reproductibilitatea studiilor genomice.

---

## Cum am formulat căutările PubMed pentru variante

Pentru variantele care au un rsID (ex: rs17878362), am căutat direct utilizând acest identificator unic, deoarece rsID-urile sunt standardizate și referite frecvent în literatură. Pentru variantele fără rsID, am construit o căutare combinată folosind poziția cromozomială (ex: "chr17:7579472 AND TP53") pentru a identifica studii care discută acea poziție specifică în contextul genei TP53. Această abordare dublă maximizează șansele de a găsi literatură relevantă, fie prin identificatori standardizați, fie prin coordonate genomice, acoperind atât studiile recente care folosesc nomenclatura actuală, cât și publicațiile mai vechi.

---

## Note suplimentare

- Toate scripturile sunt documentate și testate
- Am utilizat Biopython pentru accesarea bazelor de date NCBI
- Vizualizările includ distribuții de lungimi și scoruri Phred pentru o evaluare rapidă a calității
- Fișierele de output sunt generate automat cu prefixul corespunzător handle-ului GitHub

---
