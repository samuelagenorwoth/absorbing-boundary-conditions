"""Compare the four boundary conditions with a reflection-free reference.

The reference uses the same scheme on a much larger domain, so no reflection
reaches the window [0, 10]^2 before t = 10. Inside the window, any difference
from the reference is caused by the boundary condition alone.
"""
from pathlib import Path
import numpy as np
from wave import solve, BOUNDARY_CONDITIONS, NAMES

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    print("Reference solution (large domain)...")
    ref = solve("dirichlet", margin=8.0)
    scale = np.linalg.norm(ref["frames"][0])        # size of the initial pulse

    data = {"t": ref["t"], "x": ref["x"], "ref": ref["frames"]}
    for bc in BOUNDARY_CONDITIONS:
        sol = solve(bc)
        err = np.linalg.norm(sol["frames"] - ref["frames"], axis=(1, 2)) / scale
        data[bc] = sol["frames"]
        data[f"err_{bc}"] = err
        print(f"{NAMES[bc]:28s} largest reflection error: {err.max():.4f}")

    RESULTS.mkdir(exist_ok=True)
    np.savez_compressed(RESULTS / "compare.npz", **data)
    print(f"Saved {RESULTS / 'compare.npz'}")


if __name__ == "__main__":
    main()
