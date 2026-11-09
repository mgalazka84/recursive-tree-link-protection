# Link protection in random recursive trees

Computational materials for *The value of topology information for link protection in random recursive trees* by **Marek Gałązka** and **Hanna Wdowicka**.

The experiments compare two ways to protect a fixed number of links before an outage. One uses installation order. The other observes the tree and maximizes the expected number of vertices that retain a path to the root. Each reported gain is relative to service with ordinary, unprotected links.

## Reproduce the calculations

From the repository root, use a C++17 compiler, GNU Make and Python 3 with matplotlib. The recorded results were generated with GCC 13.3.0, Python 3.12.14 and matplotlib 3.10.8.

~~~sh
make verify
make reproduce
~~~

The first command checks the optimizer against exhaustive subset enumeration on every recursive tree of size 2 through 8, and checks the service-threshold ordering through size 10. The second command reruns the checks and all experiments, updates the files in `results/` and regenerates both vector figures in `figures/`. Use `make figures` to plot the supplied results directly.

## Files and experiments

| File | Contents |
|---|---|
| `code/experiments.cpp` | C++ simulation, exact conditional optimizer, exhaustive checks and limiting mass sampler |
| `code/figures.py` | Python script for the network example and budget comparison |
| `results/verification.txt` | Enumeration counts and maximum numerical discrepancies |
| `results/single_link.csv` | 40,000 independent trees of size 10,000, one protected link, failure probability 1/log(10,000) |
| `results/finite_budgets.csv` | 35 combinations of network size, failure probability and budget |
| `results/limiting_constants.csv` | Direct estimates for 11 budgets, using 20,000 independent limiting trees |

The executable also accepts `verify`, `one`, `finite` and `limit` as individual modes. Each writes its output to standard output. The seeds are fixed in the source and recorded in the CSV files. The map gain is averaged over independent trees. The age-based denominator is computed analytically. No external data are needed.
