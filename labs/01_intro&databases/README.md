# Săptămâna 1 — Databases & GitHub 

## Scopuri:
1) Toți studentii rulează același mediu (Codespaces sau Docker + Jupyter).
2) Toți studentii trec o dată prin fluxul GitHub (fork → branch → PR) cu un task scurt.

> În săptămânile 2–12 vom **reutiliza datele** colectate devreme (GEO/TCGA/NCBI/Ensembl) pentru aliniere, NGS, filogenie, co-expresie, ML, etc. (vezi calendarul). 

---

## Partea 0 — Verificare mediu 

Alege **una**:

### A) Codespaces
1. Deschide un Codespace pe branch `main`.
2. Terminal:
   ```bash
   python labs/00_smoke/smoke.py
   ```
Așteptat: ok.
### B) Docker local (Windows PowerShell)
```powershell
docker pull ghcr.io/bozdogalex/bioinf-y4-lab:base
docker run -it --rm -p 8890:8888 -v "${PWD}:/work" -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base `
  bash -lc "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token='' --allow-root"
```
Deschide http://localhost:8890/lab și rulează print("ok").

**Notă (căi cu & pe Windows)**: când rulezi scripturi din acest folder, pune calea între ghilimele:
python "labs/01_intro&databases/demo_entrez_brca1.py"

## Partea 1 — Task PR 
Citeste [git-workflow.md](../../docs/git-workflow.md)
Adaugă handle-ul tău GitHub în labs/01_intro&databases/roster/handles.csv (format: Nume Prenume,github_handle) prin Pull Request.

Pași:

- fork repo → creează branch feat/roster-<handle>.
- edit labs/01_intro&databases/roster/handles.csv → adaugă o singură linie.
- git commit -m "Add <handle> to roster" → git push .
- Completează checklist-ul PR (șablonul week1_roster.md).

## Partea 2 - Demo / Exercitii
**Rulati**
- `demo01_entrez_brca1.py` — căutare + descărcare BRCA1 (GenBank) și sumar GC.
- `demo02_seq_ops.py` — operații de bază pe secvență (transcriere, traducere, GC, motif).
- `dem03_dbsnp.py` — dbSNP: interogare rapidă și sumar.

**Completati si rulati**
- `ex01_multifasta_gc.py` — **completați TODO-urile de descărcare FASTA (Entrez) și rulați calculul GC pe fișierul DESCĂRCAT**. 
  - Exemplu (interogare):  
    `python labs/01_intro&databases/ex01_multifasta_gc.py --email student@example.com --query "TP53[Gene] AND Homo sapiens[Organism]" --retmax 3 --out data/work/<handle>/lab01/my_tp53.fa`
  - Exemplu (accession):  
    `python labs/01_intro&databases/ex01_multifasta_gc.py --email student@example.com --accession NM_000546 --out data/work/<handle>/lab01/nm000546.fa`

_Notă: salvați fișierele proprii în `data/work/<handle>/lab01/`. Acest folder este ignorat de git — NU încărcați datele în repo._

## Livrabile
În PR trebuie să apară:

1. O linie nouă în `labs/01_intro&databases/roster/handles.csv`.
2. Un fișier `labs/01_intro&databases/<github_handle>_notes.md` care conține:
   - Confirmarea că ai rulat toate scripturile demo și exercițiul.
   - Un rezultat simplu observat (ex.: “GC fraction = 0.47”).
3. Exercitiul completat, salvat in:
```bash
labs/01_intro&databases/submissions/<github_handle>/ex01_multifasta_gc.py
```
4. Completarea checklist-ului din șablonul PR.

## Săptămâna următoare
- Vom folosi fișierele FASTA descarcate pentru a realiza alinieri de secvențe (global și local, NW/SW).
- Vezi [Săptămâna 2 — Sequence Alignment](../02_alignment/README.md)

## Competențe: 
- Rularea mediului reproducibil (Codespaces/Docker).
- Deschiderea și completarea corectă a unui PR (fork → branch → PR).
- Primele interogări și operații de bază pe secvențe biologice.


## Resurse : 
- [Fișa laborator](../../docs/lab_onepagers/01_intro&databases.md) 
- [NCBI](https://www.ncbi.nlm.nih.gov/)  
- [Ensembl Genome Browser](https://www.ensembl.org/)  
- [dbSNP (NCBI)](https://www.ncbi.nlm.nih.gov/snp/)  
- [TCGA (The Cancer Genome Atlas) Portal](https://portal.gdc.cancer.gov/)  
- Carte: [Pevsner, *Bioinformatics and Functional Genomics*, 3rd ed., Wiley Blackwell, 2015](https://genetics.elte.hu/oktatasi_anyag/archivum/bioinfo/Bioinformatika_2018-2019/book.pdf)  
- Carte: [Lesk, *Introduction to Bioinformatics*, 5th ed., Oxford University Press, 2019](https://edscl.in/pluginfile.php/3340/mod_folder/content/0/Introduction%20To%20Bioinformatics.pdf?forcedownload=1)  



