import argparse
from os import path

import matplotlib.pyplot as plt
import numpy as np


def file_type(string):
    if path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)


parser = argparse.ArgumentParser(description="Plot interpolation results.")
parser.add_argument(
    "points_file",
    type=file_type,
    help="path to CSV file containing the input points ('knots') on which the interpolation is based.",
)

parser.add_argument(
    "results_file",
    type=file_type,
    help="path to CSV file containing the interpolation results alongside the actual function values.",
)


def main(args: argparse.Namespace) -> None:
    data = np.loadtxt(args.results_file, delimiter=",")
    points = np.loadtxt(args.points_file, delimiter=",")

    fig = plt.figure()
    ax = fig.add_subplot(111)

    n_func = points.shape[1] - 1

    cmap = plt.get_cmap("tab10")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Interpolation Results")
    for i in range(n_func):
        # Reference function
        ax.plot(data[:, 0], data[:, 1 + i], color=cmap(i), label=f"$f_{i}(x)$")
        # Interpolation result
        ax.plot(
            data[:, 0],
            data[:, 1 + n_func + i],
            marker="x",
            color=cmap(i),
            label=f"$\\tilde{{f_{i}}}(x)$",
        )
        # Knots
        ax.scatter(points[:, 0], points[:, 1 + i], marker="^", color="r", zorder=100)

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main(parser.parse_args())
