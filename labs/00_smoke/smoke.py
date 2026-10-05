# labs/00_smoke/smoke.py
"""Import smoke test to verify core libs inside the container."""

import duckdb
import matplotlib
import mlflow
import networkx as nx
import numpy as np
import pandas as pd
import polars as pl
import sklearn
import Bio


def main() -> None:
    # minimal no-op uses to satisfy flake8 and verify imports work
    _ = np.array([1, 2, 3]).sum()
    _ = pd.DataFrame({"x": [1, 2]}).shape
    _ = pl.DataFrame({"x": [1, 2]}).height
    _ = duckdb.connect(":memory:").execute("select 1").fetchone()
    _ = sklearn.__version__
    _ = Bio.__version__
    g = nx.Graph()
    g.add_edge("a", "b")
    _ = len(g)
    _ = mlflow.__version__
    _ = matplotlib.__version__
    print("ok")


if __name__ == "__main__":
    main()
