"""Recreate Figures 1-8 of the 2024 report from the original scripts.

The two original scripts (wave_dirichlet.py and wave_neumann.py) are run
unchanged. Their full solution history w[n, j, k] is then plotted at
t = 1, 2, 5 and 10 in the same way as in the report (with plt.imshow, so the
picture has the same orientation as the report's figures).

Run from this folder:  python3 reproduce_report_figures.py
Takes a few minutes, because the original code uses plain Python loops.
"""
import os
import runpy
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                    # no windows; plt.show() does nothing
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TIMES = [1, 2, 5, 10]
DT = 0.01                                # time step used in the original scripts

results = {}
for name in ["dirichlet", "neumann"]:
    print(f"Running wave_{name}.py (original 2024 code)...")
    with tempfile.TemporaryDirectory() as tmp:
        cwd = os.getcwd()
        os.chdir(tmp)                    # keep the scripts' own PNG files out of the repo
        try:
            results[name] = runpy.run_path(str(HERE / f"wave_{name}.py"))["w"]
        finally:
            os.chdir(cwd)
        plt.close("all")

fig, axes = plt.subplots(2, 4, figsize=(16, 7.5))
number = 1
for row, name in enumerate(["dirichlet", "neumann"]):
    w = results[name]
    for col, t in enumerate(TIMES):
        ax = axes[row, col]
        im = ax.imshow(w[round(t / DT)], extent=[0, 10, 0, 10])
        fig.colorbar(im, ax=ax, shrink=0.8)
        ax.set_title(f"Figure {number}: {name.capitalize()}, t = {t}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        number += 1

fig.tight_layout()
out = HERE / "report-figures-1-8.png"
fig.savefig(out, dpi=120)
print(f"Saved {out}")
