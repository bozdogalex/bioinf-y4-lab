from pathlib import Path
import pandas as pd

HANDLE = "olezhka0809"

RAW_CSV = Path(f"/workspaces/bioinf-y4-lab/data/work/olezhka0809/lab06/GSE115469_Data.csv")
OUT_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/olezhka0809/expression_matrix.csv")

# cât să păstrăm (modifici dacă vrei)
MAX_SAMPLES = 1000     # max 1000 celule
MAX_GENES   = 5000     # max 5000 gene

if __name__ == "__main__":
    print(f"[INFO] Citim fișierul brut: {RAW_CSV}")
    df = pd.read_csv(RAW_CSV, index_col=0)

    print("[DEBUG] formă inițială:", df.shape)  # (20007, 8444)

    # Subsample pe probe (coloane)
    if df.shape[1] > MAX_SAMPLES:
        df = df.iloc[:, :MAX_SAMPLES]
        print("[INFO] Am redus numărul de probe la:", df.shape[1])

    # Subsample naiv pe gene (primele MAX_GENES)
    if df.shape[0] > MAX_GENES:
        df = df.iloc[:MAX_GENES, :]
        print("[INFO] Am redus numărul de gene la:", df.shape[0])

    print("[INFO] Formă finală pentru expression_matrix:", df.shape)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV)
    print(f"[OK] Am salvat matricea de expresie în: {OUT_CSV}")