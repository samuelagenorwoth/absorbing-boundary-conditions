"""Leapfrog finite-difference solver for the 2D wave equation u_tt = u_xx + u_yy
with reflecting or absorbing boundary conditions.

Boundary conditions (the same one is used on all four edges):
    'dirichlet'  u = 0                                   (reflects, sign flips)
    'neumann'    du/dn = 0                               (reflects, same sign)
    'em1'        Engquist-Majda first order,  u_t - u_x = 0 at x = 0
    'em2'        Engquist-Majda second order, u_xt - u_tt + u_yy/2 = 0 at x = 0
"""
import numpy as np

BOUNDARY_CONDITIONS = ["dirichlet", "neumann", "em1", "em2"]
NAMES = {
    "dirichlet": "Dirichlet",
    "neumann": "Neumann",
    "em1": "Engquist-Majda, 1st order",
    "em2": "Engquist-Majda, 2nd order",
}


def laplacian(u):
    """Five-point Laplacian times h^2, zero on the edges."""
    d = np.zeros_like(u)
    d[1:-1, 1:-1] = (u[2:, 1:-1] + u[:-2, 1:-1] + u[1:-1, 2:] + u[1:-1, :-2]
                     - 4 * u[1:-1, 1:-1])
    return d


# Each edge is turned into the LEFT edge (row 0), updated, and turned back.
# Arrays are stored as u[i, j] = u(x_i, y_j).
_TO_LEFT = {
    "left":   (lambda a: a,             lambda a: a),
    "right":  (lambda a: a[::-1, :],    lambda a: a[::-1, :]),
    "bottom": (lambda a: a.T,           lambda a: a.T),
    "top":    (lambda a: a.T[::-1, :],  lambda a: a[::-1, :].T),
}


def apply_edge(u_new, u, u_old, side, bc, r):
    """Set the boundary values of u_new (time n+1) on one edge, in place.

    The formulas are written for the left edge x = 0, where outgoing waves
    travel in the -x direction. r = dt/h is the Courant number.
    """
    to_left, _ = _TO_LEFT[side]
    # These are views into the original arrays, so writing to wn writes to u_new.
    wn, w, wo = to_left(u_new), to_left(u), to_left(u_old)
    c1 = (r - 1) / (r + 1)

    if bc == "dirichlet":
        wn[0, :] = 0
    elif bc == "neumann":
        wn[0, :] = wn[1, :]
    elif bc in ("em1", "em2"):
        # First order (Mur 1981): u_t - u_x = 0 centred between the first two
        # grid lines and between the two time levels.
        wn[0, :] = w[1, :] + c1 * (wn[1, :] - w[0, :])
        if bc == "em2":
            # Second order on the edge; the corners keep the first-order value.
            k = slice(1, -1)

            def d2(v):  # second difference along the edge
                return v[2:] - 2 * v[1:-1] + v[:-2]

            wn[0, k] = (-wo[1, k]
                        + c1 * (wn[1, k] + wo[0, k])
                        + 2 / (r + 1) * (w[0, k] + w[1, k])
                        + r**2 / (2 * (r + 1)) * (d2(w[0, :]) + d2(w[1, :])))
    else:
        raise ValueError(f"Unknown boundary condition {bc!r}")


def solve(bc, L=10.0, h=0.05, r=0.5, T=10.0, dt_save=0.1, margin=0.0):
    """Solve on the window [0, L]^2 with boundary condition bc on all edges.

    margin adds extra space around the window. With a large margin no
    reflection reaches the window before time T, which gives a reference
    solution. Returns a dict with x, t and frames[s, i, j] on the window.
    """
    dt = r * h
    m = round(margin / h)                    # extra grid points on each side
    n_win = round(L / h)                     # the window has n_win + 1 points
    x = h * np.arange(-m, n_win + m + 1)     # grid aligned with the window
    win = slice(m, m + n_win + 1)
    X, Y = np.meshgrid(x, x, indexing="ij")

    # Initial condition from the 2024 report: a Gaussian pulse at (5, 5), at rest
    u0 = np.exp(-((X - 5) ** 2 + (Y - 5) ** 2))
    u_old = u0
    u = u0 + 0.5 * r**2 * laplacian(u0)      # second-order accurate first step

    n_steps = round(T / dt)
    every = round(dt_save / dt)
    frames, times = [u0[win, win].astype(np.float32)], [0.0]

    for n in range(1, n_steps):
        u_new = 2 * u - u_old + r**2 * laplacian(u)      # interior points
        for side in ("left", "right", "bottom", "top"):
            apply_edge(u_new, u, u_old, side, bc, r)
        u_old, u = u, u_new                              # u is now at (n+1)*dt
        if (n + 1) % every == 0:
            frames.append(u[win, win].astype(np.float32))
            times.append((n + 1) * dt)

    return {"bc": bc, "x": x[win], "t": np.array(times), "frames": np.array(frames)}
