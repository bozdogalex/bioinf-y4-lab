# BIOINF-Y4 — Bioinformatica și Genomică Funcțională (2026–2027, Bachelor, Year 4) 

[![Open in Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/bozdogalex/bioinf-y4-lab?quickstart=1)
![CI](https://github.com/bozdogalex/bioinf-y4-lab/actions/workflows/ci.yml/badge.svg?branch=main)


> Laboratoare la nivel de licență (anul IV), combinând bioinformatica clasică cu metode moderne de învățare automată, rețele și GenAI.  
> Mediul este CPU-only și identic între Codespaces și Docker prin imaginea preconstruită `ghcr.io/bozdogalex/bioinf-y4-lab:base`.


## Student submissions — 2026–2027

**Teaching materials are public; assessed submissions are private.** Create your own private repository from a current snapshot, invite `bozdogalex`, and open lab PRs inside that repository. Submit each PR link and commit SHA through the university LMS. Do not submit exercise solutions or roster entries as PRs to this public repository.

Start with the **[private repository and submission guide](docs/git-workflow.md)**, then [set up your environment](docs/onboarding.md). Open Codespaces on your private repository when working on assessed exercises; the badge above opens the public teaching repository for browsing/demos.

## Labs (index)

- 01 — Databases & GitHub: [labs/01_intro&databases](labs/01_intro&databases)
- 02 — Sequence Alignment: [labs/02_alignment](labs/02_alignment)
- 03 — NGS: [labs/03_formats&NGS](labs/03_formats&NGS)
- 04 — Phylogenetics: [labs/04_phylogenetics](labs/04_phylogenetics)
- 05 — Clustering: [labs/05_clustering](labs/05_clustering)
- 06 — WGCNA + Diseasome: [labs/06_wgcna](labs/06_wgcna)
- 07 — Network Viz & GNN: [labs/07_network_viz](labs/07_network_viz)
- 08 — Machine Learning: [labs/08_ML_flower](labs/08_ML_flower)
- 09 — Drug Repurposing: [labs/09_repurposing](labs/09_repurposing)
- 10 — Integrative Genomics: [labs/10_integrative](labs/10_integrative)
- 11 — Multi‑omics + Quantum (planned; materials not yet published): [labs/11_multiomics](labs/11_multiomics)
- 12 — Generative AI (planned; materials not yet published): [labs/12_genAI](labs/12_genAI)
- Assignment Presentations

---

> Full onboarding (screenshots, tips): **[docs/onboarding.md](docs/onboarding.md)**

---

## Repo map

- `labs/` — all weekly lab content
- `docs/` — onboarding, ANIS pack (before/after, one‑pagers, screenshots)
- `mlops/` — MLflow helpers
- `.devcontainer/` — Codespaces/Devcontainer (pulls GHCR image)
- `.github/workflows/` — CI + image publish
- `Dockerfile`, `requirements.txt` — env definition
- `dev.ps1`, `Makefile` — local helpers

---
# Docs

Supporting material and submission pack live under `docs/`:

- [Onboarding](docs/onboarding.md) — Codespaces & Docker setup, smoke test, troubleshooting
- [One-pagers](docs/lab_onepagers/) — summaries of labs
- [Changelog](docs/changelog.md) — changes across versions
- [Policies](docs/policies.md) — third-party license references, repository policies
- [Resources](docs/resources.md) — recommended readings/tutorials
- [GA4GH](docs/GA4GH_primer.md) - ethical & technical standards for sharing biomedical data
- [GDPR and Data policy](docs/GDPR_and_DataPolicy.md) 
---

### MLflow 
- See **[docs/mlflow.md](docs/mlflow.md)** for quick start, Codespaces UI, and troubleshooting.
---

## Contributing / Policies / Citation

- [Contributing](CONTRIBUTING.md) — contribution rules & PR tips  
- [Citation](citation.cff)  — how to cite this work  
- [Changelog](docs/changelog.md) — changelog (linked from releases)


