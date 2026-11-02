"""Print the six-gene slide example using all 48 synthetic samples.

The submatrix is an excerpt from the full correlation network. Restricting
this printout does not change how correlations are estimated across samples.
"""
import argparse

from ex01_gce_networks import (
    DEFAULT_INPUT, read_expression_matrix, correlation_matrix,
    adjacency_from_correlation,
)

GENES = ["G001", "G002", "G032", "G033", "G061", "G062"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threshold", type=float, default=0.7)
    parser.add_argument("--network", choices=["positive", "unsigned"], default="positive")
    args = parser.parse_args()
    expression = read_expression_matrix(DEFAULT_INPUT)
    corr = correlation_matrix(expression, method="spearman").loc[GENES, GENES]
    adjacency = adjacency_from_correlation(corr, args.threshold, weighted=True,
                                           use_abs=args.network == "unsigned")
    print("Corelații Spearman între șase gene; fiecare profil are 48 de probe:")
    print(corr.round(3).to_string())
    print(f"\nAdiacență ponderată, {args.network}, prag {args.threshold}:")
    print(adjacency.round(3).to_string())
    print("\nExemplu sintetic: o muchie exprimă asocierea profilurilor, nu reglare directă.")


if __name__ == "__main__":
    main()
