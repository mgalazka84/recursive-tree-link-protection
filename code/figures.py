"""Rebuild the two vector figures from the reported model and CSV files."""
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import FixedLocator, ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.65, "pdf.fonttype": 42,
    "axes.labelcolor": "0.12", "text.color": "0.12",
})


def records(name):
    with (ROOT / "results" / name).open(newline="") as source:
        return [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(source)]


def tree_example():
    pos = {1: (0, 1.8), 2: (-0.92, 0.9), 3: (0.6, 0.9),
           7: (-0.92, 0.0), 4: (-0.14, 0.0), 5: (0.6, 0.0), 6: (1.34, 0.0)}
    edges = [(1, 2), (1, 3), (2, 7), (3, 4), (3, 5), (3, 6)]
    fig, axes = plt.subplots(1, 2, figsize=(6.25, 1.94), layout="constrained")
    for ax, protected, title in zip(axes, [(1, 2), (1, 3)],
                                    ["(a) Installation order", "(b) Observed tree"]):
        for u, v in edges:
            a, b = pos[u], pos[v]
            selected = (u, v) == protected
            ax.plot([a[0], b[0]], [a[1], b[1]],
                    lw=3.7 if selected else 1.15,
                    color="0.12" if selected else "0.58", solid_capstyle="round")
        for label, (x, y) in pos.items():
            ax.add_patch(Circle((x, y), 0.16, facecolor="0.93" if label==1 else "white",
                                edgecolor="0.18", lw=0.95, zorder=3))
            ax.text(x, y, str(label), ha="center", va="center", zorder=4, fontsize=8.2)
        ax.set(xlim=(-1.43, 1.65), ylim=(-0.24, 2.07), aspect="equal")
        ax.set_title(title, fontsize=9, pad=3)
        ax.axis("off")
    fig.savefig(OUT / "tree_example.pdf", bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


def budget_figure():
    finite = records("finite_budgets.csv")
    limits = records("limiting_constants.csv")
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.0, 2.63), layout="constrained")
    for ax in (left, right):
        ax.set_xscale("log", base=2)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.grid(axis="y", color="0.89", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlabel("Protected links, $k$")
        ax.tick_params(direction="out", length=3, color="0.4")
    for n, marker, style, gray in [(128, "s", "--", "0.47"), (1024, "^", "-.", "0.3")]:
        rows = [r for r in finite if r["n"] == n and r["c"] == 1]
        left.errorbar([r["k"] for r in rows], [100*(r["ratio"]-1) for r in rows],
                      yerr=[196*r["ratio_se"] for r in rows], color=gray,
                      fmt=marker, linestyle=style, lw=1, ms=3.6, mfc="white",
                      capsize=1.8, label=f"$n={n}$")
    rows = [r for r in limits if r["k"] <= 64]
    left.errorbar([r["k"] for r in rows], [100*(r["ratio"]-1) for r in rows],
                  yerr=[196*r["ratio_se"] for r in rows], color="0.05", fmt="o-",
                  lw=1.15, ms=3.2, capsize=1.8, label="Fixed-budget limit")
    left.set_ylabel("Extra gain from the map (%)")
    left.set_title("(a) Relative benefit", fontsize=9.2)
    left.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8, 16, 32, 64]))
    left.set_ylim(0, 28)
    left.legend(fontsize=7, loc="lower left", frameon=False, handlelength=2.2)
    ks = [r["k"] for r in limits]
    right.errorbar(ks, [r["gap"] for r in limits], yerr=[1.96*r["map_se"] for r in limits],
                   color="0.08", fmt="o-", lw=1.1, ms=3.4, capsize=1.8, label="Estimated gap")
    right.plot(ks, [math.log1p(r["k"])-r["age_constant"] for r in limits],
               color="0.45", ls=":", lw=1.3, label="Upper bound")
    right.axhline(1-0.5772156649015329, color="0.25", ls="--", lw=1.0, label="$1-\\gamma$")
    right.set_title("(b) Additive benefit in the limit", fontsize=9.2)
    right.set_ylabel(r"Gain gap, $v_k-a_k$")
    right.xaxis.set_major_locator(FixedLocator([1, 4, 16, 64, 256, 1024]))
    right.set_ylim(0.09, 0.465)
    right.legend(fontsize=7, loc="lower right", frameon=False, handlelength=2.2)
    fig.savefig(OUT / "budget_comparison.pdf", bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


if __name__ == "__main__":
    tree_example()
    budget_figure()
