# MLflow — Minimal Standard (BIOINF-Y4)

This lab uses **local, file-based MLflow tracking** so everything works offline and inside the container.

---

## Quick Start

**Log a demo run**
```bash
python labs/00_smoke/mlflow_smoke.py --experiment "BIOINF-Y4 Demo"
```

Expected output:
```bash
OK: metrics/params/artifact logged to ./mlruns (file store)
```

Open the UI
```bash
mlflow ui --backend-store-uri file://$PWD/mlruns --host 0.0.0.0 --port 5000
```

Open port 5000 → browse to the UI → experiment BIOINF-Y4 Demo → run mlflow_smoke.
You’ll see:
- Metrics: auc over steps 0..3
- Params: seed=42, model=toy
- Artifacts: artifact.txt

--- 

### Codespaces

Open the Codespace for this repo.

Terminal:

- Make target
```bash
make mlflow-smoke
```
- or the direct script
```bash
python labs/00_smoke/mlflow_smoke.py --experiment "BIOINF-Y4 Demo"
```
Start the UI:
```bash
mlflow ui --backend-store-uri file://$PWD/mlruns --host 0.0.0.0 --port 5000
```

In the Ports panel, open forwarded port 5000 in the browser and keep its visibility Private.

Local (Docker) mirror of CI
```bash
docker run --rm -v "$PWD":/work -w /work ghcr.io/bozdogalex/bioinf-y4-lab:base \
  python labs/00_smoke/mlflow_smoke.py --experiment "BIOINF-Y4 Demo"

mlflow ui --backend-store-uri file://$PWD/mlruns --host 127.0.0.1 --port 5000
```
### How CI checks MLflow

The CI job runs:

python labs/00_smoke/mlflow_smoke.py --experiment "BIOINF-Y4 Demo"

It does not start a server; it just verifies import + logging to ./mlruns.

### Troubleshooting

1. UI opens but shows no runs
- Ensure you’re pointing the UI to the same store:
```
--backend-store-uri file://$PWD/mlruns
```
- Rerun the smoke to create a run.
2. Port 5000 not reachable (Codespaces)
    - Ensure port 5000 is forwarded, then open its forwarded URL while signed in to your GitHub account.
