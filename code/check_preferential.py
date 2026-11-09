"""Independently enumerate the preferential model and check its one-link denominator."""
import itertools
import math


def verify():
    maximum_error = 0.0
    trees = 0
    for n in range(2, 9):
        for p in (0.25, 0.5, 0.8):
            gains = [0.0] * n
            total_probability = 0.0
            for parents in itertools.product(*(range(j) for j in range(1, n))):
                degrees = [0] * n
                depth = [0] * n
                probability = 1.0
                for child, parent in enumerate(parents, 1):
                    probability *= (1 + degrees[parent]) / (2 * child - 1)
                    degrees[parent] += 1
                    depth[child] = depth[parent] + 1
                if p == 0.25:
                    trees += 1
                total_probability += probability
                # Every non-root edge on a vertex's path has this same marginal gain.
                for vertex in range(1, n):
                    contribution = probability * (1-p) * p**(depth[vertex]-1)
                    ancestor = vertex
                    while ancestor:
                        gains[ancestor] += contribution
                        ancestor = parents[ancestor-1]
            expected_weight = 1.0
            for t in range(2, n):
                expected_weight += ((1+p)*expected_weight - 1) / (2*t-1)
            expected = (1-p) * expected_weight
            maximum_error = max(maximum_error, abs(total_probability - 1),
                                abs(expected - gains[1]), abs(expected - max(gains)))
    if not math.isfinite(maximum_error) or maximum_error > 1e-10:
        raise AssertionError(f"Preferential one-link check failed: {maximum_error}")
    print(f"preferential_trees={trees}")
    print(f"max_preferential_denominator_or_optimum_error={maximum_error:.12g}")


if __name__ == "__main__":
    verify()
