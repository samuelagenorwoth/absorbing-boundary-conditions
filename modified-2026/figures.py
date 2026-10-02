"""Make the figures and the animation for the README (saved in ../figures)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from wave import BOUNDARY_CONDITIONS, NAMES

ROOT = Path(__file__).resolve().parent.parent
FIGURES = ROOT / "figures"
INK, TEAL, ORANGE, GREY = "#15222E", "#0B6E78", "#D9822B", "#8C8C8C"
CMAP = LinearSegmentedColormap.from_list("wave", [TEAL, "white", ORANGE])
LIM = 0.15                                   # common colour scale

plt.rcParams.update({"font.size": 11, "text.color": INK,
                     "axes.labelcolor": INK, "axes.titlecolor": INK})


def show(ax, frame, x):
    return ax.imshow(frame.T, origin="lower", extent=[x[0], x[-1], x[0], x[-1]],
                     cmap=CMAP, vmin=-LIM, vmax=LIM)


def snapshots(d):
    times = [4, 6, 8]
    idx = [int(np.argmin(abs(d["t"] - tt))) for tt in times]
    rows = [(NAMES[bc], d[bc]) for bc in BOUNDARY_CONDITIONS]
    rows.append(("Reference\n(no boundary)", d["ref"]))

    fig, axes = plt.subplots(len(rows), len(times), figsize=(8, 12.5),
                             constrained_layout=True)
    for b, (name, frames) in enumerate(rows):
        for c, s in enumerate(idx):
            ax = axes[b, c]
            im = show(ax, frames[s], d["x"])
            ax.set_xticks([]); ax.set_yticks([])
            if b == 0:
                ax.set_title(f"t = {times[c]}")
            if c == 0:
                ax.set_ylabel(name.replace(", ", ",\n"), fontweight="bold")
    fig.colorbar(im, ax=axes, shrink=0.4, label="u")
    fig.savefig(FIGURES / "snapshots.png", dpi=150)
    plt.close(fig)


def reflection_error(d):
    styles = {"dirichlet": (INK, "-"), "neumann": (GREY, "--"),
              "em1": (ORANGE, "-"), "em2": (TEAL, "-")}
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    for bc in BOUNDARY_CONDITIONS:
        color, ls = styles[bc]
        ax.semilogy(d["t"], np.maximum(d[f"err_{bc}"], 1e-6), ls, color=color,
                    lw=2, label=NAMES[bc])
    ax.set_ylim(1e-4, 2)
    ax.set_xlabel("Time t")
    ax.set_ylabel("Reflection error\n(relative to the initial pulse)")
    ax.set_title("How much of the wave is reflected back into the domain")
    ax.grid(True, which="major", alpha=0.4)
    ax.legend(loc="lower right")
    fig.savefig(FIGURES / "reflection_error.png", dpi=200)
    plt.close(fig)


def animation(d):
    panels = [("Reflecting boundary", d["dirichlet"]),
              ("Absorbing boundary", d["em2"])]
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.3), dpi=70, constrained_layout=True)
    ims = []
    for ax, (title, frames) in zip(axes, panels):
        ims.append(show(ax, frames[0], d["x"]))
        ax.set_title(title, fontsize=14)
        ax.axis("off")
    label = fig.suptitle("t = 0.0")

    steps = range(0, len(d["t"]), 3)         # every 0.3 time units

    def update(s):
        for im, (_, frames) in zip(ims, panels):
            im.set_data(frames[s].T)
        label.set_text(f"t = {d['t'][s]:.1f}")
        return ims

    # Draw each frame and store it with a small colour palette (no dithering),
    # which keeps the GIF light enough for a web page.
    images = []
    for s in steps:
        update(s)
        fig.canvas.draw()
        rgb = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3])
        images.append(rgb.quantize(colors=48, method=Image.Quantize.MEDIANCUT,
                                   dither=Image.Dither.NONE))
    images[0].save(FIGURES / "reflect-vs-absorb.gif", save_all=True,
                   append_images=images[1:], duration=120, loop=0, optimize=True)
    plt.close(fig)


def main():
    d = dict(np.load(ROOT / "results" / "compare.npz"))
    FIGURES.mkdir(exist_ok=True)
    snapshots(d)
    reflection_error(d)
    animation(d)
    print(f"Figures saved in {FIGURES}")


if __name__ == "__main__":
    main()
