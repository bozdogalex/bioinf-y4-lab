# BIOINF-Y4 — Bioinformatica și Genomică Funcțională (Bachelor, Year 4) 

[![Open in Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/numipasaaa/bioinf-y4-lab?quickstart=1)
![CI](https://github.com/numipasaaa/bioinf-y4-lab/actions/workflows/ci.yml/badge.svg?branch=main)


> Laboratoare la nivel de licență (anul IV), combinând bioinformatica clasică cu metode moderne de învățare automată, rețele și GenAI.  
> Mediul este CPU-only și identic între Codespaces și Docker prin imaginea preconstruită `ghcr.io/numipasaaa/bioinf-y4-lab:base`.


## Labs (index)

- 01 — Databases & GitHub: [labs/01_intro&databases](labs/01_intro&databases)
- 02 — Sequence Alignment: [labs/02_alignment](labs/02_alignment)
- 03 — NGS: [labs/03_formats&NGS](labs/03_formats&NGS)
- 04 — Phylogenetics: [labs/04_phylogenetics](labs/04_phylogenetics)
- 05 — Clustering: [labs/05_clustering](labs/05_clustering)
- 06 — WGCNA + Diseasome: [labs/06_wgcna](labs/06_wgcna)
- 07 — Network Viz & GNN: [labs/07_network_viz](labs/07_network_viz)
- 08 — Federated Learning: [labs/08_ML_flower](labs/08_ML_flower)
- 09 — Drug Repurposing: [labs/09_repurposing](labs/09_repurposing)
- 10 — Integrative + Digital Twin: [labs/10_integrative](labs/10_integrative)
- 11 — Multi‑omics + Quantum : [labs/11_multiomics](labs/11_multiomics)
- 12 — Generative AI : [labs/12_genAI](labs/12_genAI)
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
- [One-pagers](docs/lab_onepagers/) — PDF summaries of labs
- [Screenshots](docs/screens/) — environment/UI captures (MLflow, Codespaces, Argo)
- [Changelog](docs/changelog.md) — changes across versions
- [Policies](docs/policies.md) — third-party license references, repository policies
- [Resources](docs/resources.md) — recommended readings/tutorials
- [GA4GH](docs/GA4GH_primer) - ethical & technical standards for sharing biomedical data
- [GDPR and Data policy](docs/GDPR_and_DataPolicy.md) 
---

### MLflow 
- See **[docs/mlflow.md](docs/mlflow.md)** for quick start, Codespaces UI, and troubleshooting.
---

## Contributing / Policies / Citation

- [Contributing](CONTRIBUTING.md) — contribution rules & PR tips  
- [Citation](CITATION.cff)  — how to cite this work  
- [Changelog](docs/changelog.md) — changelog (linked from releases)

