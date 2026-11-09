"""Rebuild the three vector figures from the model, its worked examples and the CSV files."""
import csv
import math
from tempfile import NamedTemporaryFile
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


def save_plot(fig, filename, padding):
    with NamedTemporaryFile(dir=OUT, suffix=".pdf", delete=False) as output:
        temp = Path(output.name)
    try:
        fig.savefig(temp, bbox_inches="tight", pad_inches=padding,
                    metadata={"CreationDate": None, "ModDate": None})
        temp.replace(OUT / filename)
    finally:
        temp.unlink(missing_ok=True)
        plt.close(fig)


def draw_tree(ax, positions, parents, highlighted, radius):
    for child, parent in parents.items():
        a, b = positions[parent], positions[child]
        chosen = child in highlighted
        ax.plot([a[0], b[0]], [a[1], b[1]], color="0.13" if chosen else "0.55",
                lw=3.3 if chosen else 1.15, solid_capstyle="round", zorder=1)
    for site, (x, y) in positions.items():
        ax.add_patch(Circle((x, y), radius, facecolor="0.92" if site == 1 else "white",
                            edgecolor="0.2", lw=0.9, zorder=2))
        ax.text(x, y, str(site), ha="center", va="center", fontsize=8, zorder=3)
    ax.set_aspect("equal")
    ax.axis("off")


def root_link_gain(parents, child, reliability):
    # Protect the root edge ending at child. A site gains only if its path uses that edge.
    gain = 0.0
    for site in parents:
        node, depth = site, 0
        while parents[node] != 1:
            node = parents[node]
            depth += 1
        if node == child:
            gain += (1 - reliability) * reliability ** depth
    return gain


def tree_example():
    fig, axes = plt.subplots(1, 2, figsize=(6.65, 2.0), layout="constrained")
    for ax, larger, panel in zip(axes, (2, 3), ("a", "b")):
        smaller = 5 - larger
        parents = {2: 1, 3: 1, 4: larger, 5: larger, 6: larger, 7: smaller}
        positions = {1: (0, 1.65), 2: (-0.87, 0.77), 3: (0.87, 0.77)}
        for label, offset in zip((4, 5, 6), (-0.56, 0, 0.56)):
            positions[label] = (positions[larger][0] + offset, -0.03)
        positions[7] = (positions[smaller][0], -0.03)
        gains = {j: root_link_gain(parents, j, 0.8) for j in (2, 3)}
        best = max(gains, key=gains.get)
        draw_tree(ax, positions, parents, {best}, radius=0.14)
        for child in (2, 3):
            ax.text(positions[child][0], -0.51, f"Link {child}   +{gains[child]:.2f}",
                    ha="center", va="center", fontsize=8,
                    weight="bold" if child == best else "normal")
        ax.set_title(f"({panel}) Relays 4–6 attach to {larger}", fontsize=9, pad=3)
        ax.set(xlim=(-1.7, 1.7), ylim=(-0.70, 1.87))
    save_plot(fig, "tree_example.pdf", 0.03)


def subtree_thresholds():
    # Child labels increase along each path. Accumulate finite subtree sizes from the leaves.
    parents = {2: 1, 3: 1, 4: 2, 5: 2, 6: 4, 7: 4}
    n = len(parents) + 1
    sizes = {j: 1 for j in range(1, n + 1)}
    for child in sorted(parents, reverse=True):
        sizes[parents[child]] += sizes[child]
    ranked = sorted(parents, key=sizes.get, reverse=True)
    top = ranked[:2]
    fig, (left, right) = plt.subplots(1, 2, figsize=(6.9, 2.35),
                                     gridspec_kw={"width_ratios": [0.94, 1.06]},
                                     layout="constrained")
    positions = {1: (0, 2.56), 2: (-0.55, 1.75), 3: (0.92, 1.75),
                 4: (-0.91, 0.94), 5: (0.18, 0.94), 6: (-1.30, 0.10), 7: (-0.54, 0.10)}
    draw_tree(left, positions, parents, set(top), radius=0.145)
    for child in top:
        x, y = positions[child]
        left.text(x - 0.31, y + 0.04, f"{sizes[child]}/{n}", ha="right", va="center", fontsize=8)
    left.set(xlim=(-1.82, 1.37), ylim=(-0.12, 2.78))
    left.set_title("(a) Two nested subtrees", fontsize=9, pad=3)

    fractions = [sizes[j] / n for j in parents]
    breaks = [0, *sorted(set(fractions)), 1]
    midpoints = [(a + b) / 2 for a, b in zip(breaks, breaks[1:])]
    counts = [sum(s >= x for s in fractions) for x in midpoints]
    capped = [min(len(top), count) for count in counts]
    right.stairs(counts, breaks, baseline=None, color="0.43", linestyle="--",
                 linewidth=1.4, label="All qualifying subtrees", zorder=4)
    right.stairs(capped, breaks, fill=True, color="0.84", label="Area for the two largest")
    right.stairs(capped, breaks, baseline=None, color="0.12", linewidth=1.7)
    # Included endpoints correspond to the weak threshold in the paper.
    for x in sorted(set(fractions)):
        h = sum(s >= x for s in fractions)
        right.plot(x, h, "o", color="0.4", markersize=3.1, zorder=5)
        right.plot(x, min(len(top), h), "o", color="0.12", markersize=3.1, zorder=6)
    right.set(xlim=(0, 1), ylim=(-0.12, 6.7), xlabel="Threshold (fraction of all sites)",
              ylabel="Number of subtrees", yticks=[0, 1, 2, 4, 6])
    right.set_xticks(breaks, ["0", "1/7", "3/7", "5/7", "1"])
    right.tick_params(direction="out", length=3, color="0.4")
    right.legend(loc="upper right", fontsize=7.2, frameon=False, handlelength=2)
    right.set_title("(b) Count above each threshold", fontsize=9, pad=3)
    save_plot(fig, "subtree_thresholds.pdf", 0.035)


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
    save_plot(fig, "budget_comparison.pdf", 0.025)


if __name__ == "__main__":
    tree_example()
    subtree_thresholds()
    budget_figure()
