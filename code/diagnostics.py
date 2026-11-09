"""Compute finite and limiting comparisons from the recorded uniform-tree draws."""
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


def baseline_for_reliability(n, p):
    """Expected unprotected number served in the uniform tree."""
    return math.exp(math.fsum(math.log1p(p / t) for t in range(1, n)))


def fixed_risk_bounds():
    """Approximate pointwise intervals for the limiting one-link ratio."""
    for row in read_rows("one_link_sensitivity.csv"):
        if row["growth"] != "uniform" or row["regime"] != "fixed":
            continue
        n = int(row["n"])
        q = float(row["q"])
        p = 1 - q
        base = baseline_for_reliability(n, p)
        squared_base = baseline_for_reliability(n, p * p)
        tail_bound = (1 + p) / base * math.sqrt(
            squared_base / (2 * p - p * p) + (1 + p) / (2 * p)
        )
        ratio = float(row["ratio"])
        se = float(row["ratio_se"])
        yield (n, row["trees"], row["seed"], p, q, ratio, se, tail_bound,
               ratio - 1.96 * se, ratio + 1.96 * se + tail_bound)


def write_rows(filename, fields, records):
    with (RESULTS / filename).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields.split())
        for record in records:
            writer.writerow([format(x, ".13g") if isinstance(x, float) else x for x in record])


def main():
    write_rows("normalizations.csv",
               "source n k c trees seed baseline age_over_nq age_limit gap_over_nq_exp "
               "gap_over_nq_exp_se gap_over_Bq gap_over_Bq_se fixed_k_gap_limit limit_se",
               comparisons())
    write_rows("fixed_risk_limits.csv",
               "n trees seed p q finite_ratio ratio_se finite_size_bound "
               "limit_lower_approx95 limit_upper_approx95", fixed_risk_bounds())


if __name__ == "__main__":
    main()
