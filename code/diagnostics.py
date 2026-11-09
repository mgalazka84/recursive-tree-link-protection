"""Compute comparisons to the critical-window constants from the recorded draws."""
import csv
import math
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[1] / "results"
GOLOMB_DICKMAN = 0.6243299885435508


def read_rows(filename):
    with (RESULTS / filename).open(newline="") as stream:
        return list(csv.DictReader(stream))


def comparisons():
    limits = {int(row["k"]): row for row in read_rows("limiting_constants.csv")}
    observations = [(row, row["k"], float(row["c"]), source)
                    for source in ("finite_budgets.csv", "single_link.csv")
                    for row in read_rows(source)]
    for row in read_rows("one_link_sensitivity.csv"):
        if row["regime"] == "critical" and row["growth"] == "uniform":
            observations.append((row, "1", 1.0, "one_link_sensitivity.csv"))

    for row, budget, c, source in observations:
        n, k = int(row["n"]), int(budget)
        q = c / math.log(n)
        baseline = n * math.exp(math.fsum(math.log1p(-q / j) for j in range(2, n + 1)))
        age, full, se = (float(row[key]) for key in ("age_exact", "map_mean", "map_se"))
        scale = n * q * math.exp(-c)
        limit = limits[k]
        gap_limit = GOLOMB_DICKMAN - 0.5 if k == 1 else float(limit["gap"])
        limit_se = 0.0 if k == 1 else float(limit["map_se"])
        yield (source, n, k, c, row["trees"], row["seed"], baseline,
               age / (n*q), math.exp(-c) * float(limit["age_constant"]),
               (full-age) / scale, se / scale, (full-age) / (baseline*q),
               se / (baseline*q), gap_limit, limit_se)


def main():
    fields = ("source n k c trees seed baseline age_over_nq age_limit gap_over_nq_exp "
              "gap_over_nq_exp_se gap_over_Bq gap_over_Bq_se fixed_k_gap_limit limit_se").split()
    with (RESULTS / "normalizations.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        for record in comparisons():
            writer.writerow([format(x, ".13g") if isinstance(x, float) else x for x in record])


if __name__ == "__main__":
    main()
