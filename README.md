# Link protection in random recursive trees

Computational materials for *The value of topology information for link protection in random recursive trees* by **Marek Gałązka** and **Hanna Wdowicka**.

The experiments compare two ways to protect a fixed number of links before an outage. One uses installation order. The other observes the tree and maximizes the expected number of vertices that retain a path to the root. Each reported gain is relative to service with ordinary, unprotected links.

## Reproduce the calculations

From the repository root, use a C++17 compiler, GNU Make and Python 3 with matplotlib. The recorded results were generated with GCC 13.3.0, Python 3.12.14 and matplotlib 3.10.8.

~~~sh
make verify
make reproduce
~~~

The first command checks the optimizer against exhaustive subset enumeration on every recursive tree of size 2 through 8, and checks the service-threshold ordering through size 10. An independent Python enumeration also checks every candidate one-link label in the preferential model through size 8. The second command reruns the checks and all experiments, updates the files in `results/` and regenerates all three vector figures in `figures/`. Use `make figures` to plot the supplied results directly.

## Files and experiments

| File | Contents |
|---|---|
| `code/experiments.cpp` | C++ simulation, exact conditional optimizer, exhaustive checks and limiting mass sampler |
| `code/figures.py` | Python script for the two possible relay maps, the subtree-threshold example and the budget comparison |
| `results/verification.txt` | Enumeration counts and maximum numerical discrepancies |
| `results/single_link.csv` | 40,000 independent trees of size 10,000, one protected link, failure probability 1/log(10,000) |
| `results/finite_budgets.csv` | 35 combinations of network size, failure probability and budget |
| `results/limiting_constants.csv` | Direct estimates for 11 budgets, using 20,000 independent limiting trees |
| `results/one_link_sensitivity.csv` | Two growth rules, three constant failure risks and a critical risk, at four sizes from 128 to 65,536 |
| `code/diagnostics.py`, `results/normalizations.csv` | Compares finite and limiting gains at each budget using both the exact finite baseline and its large-network approximation |
| `code/check_preferential.py`, `results/preferential_verification.txt` | Exhaustive, probability-weighted check of the preferential one-link optimum and its exact denominator |

The sensitivity runs use independent trees for each size and growth rule. Within a draw we reuse the same tree across risks and compare the fully informed one-link choice with protecting the largest branch at the root by vertex count (earliest label breaks ties). The `paired_gap_se` column is the standard error of the within-tree difference between those choices. The exact oldest-link mean supplies the ratio denominator under either growth rule. In `normalizations.csv` a zero limiting standard error at budget one denotes the analytic Golomb–Dickman constant. For larger budgets the limiting comparison is a separate estimate based on 20,000 draws.

The executable also accepts `verify`, `one`, `finite`, `limit` and `sensitivity` as individual modes. Each writes its output to standard output. The seeds are fixed in the source and recorded in the CSV files. The map gain is averaged over independent trees. The age-based denominator is computed analytically. No external data are needed.
