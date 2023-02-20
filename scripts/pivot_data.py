import argparse
import re
from io import StringIO
from os import path

import pandas as pd


def check_file(string: str) -> str:
    if path.isfile(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"'{string}' is not a file.")


def check_directory(string: str) -> str:
    if path.exists(path.dirname(string)):
        return string
    else:
        raise argparse.ArgumentTypeError(
            f"Directory '{path.dirname(string)}' doesn't exist."
        )


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Utility for pivoting the data csv file."
    )
    parser.add_argument(
        "input_file",
        type=check_file,
        help="path to CSV file containing the input data",
    )

    parser.add_argument(
        "output_file",
        type=check_directory,
        help="path to CSV file where the pivoted csv data should be stored",
    )
    return parser


def main(args: argparse.Namespace) -> None:
    if path.exists(args.output_file):
        print(
            f"The provided output file {args.output_file} already exists.\nShould I overwrite it? (y/n)"
        )
        response = input()
        if response != "y":
            return

    with open(args.input_file, "r") as input_file:
        # Extract the original file header
        comments = "".join(
            [l for l in input_file.readlines() if l.startswith("#")][-3:]
        )
        # And update the input dimension
        comments = comments.replace("InputDimension=3", "InputDimension=2")
        input_file.seek(0)

        # Load the csv file
        df = pd.read_csv(input_file, comment="#", header=None)
        df.rename(
            columns={0: "taperAngle", 1: "stickingProbability", 2: "timestep"},
            inplace=True,
        )

        # Pivot the table
        df.pivot(index=["taperAngle", "stickingProbability"], columns=["timestep"])

        # Store the pivoted table in a string buffer
        buffer = StringIO()
        df.to_csv(buffer, header=False, index=True)

        # Write the modified header as well as the pivoted table to the output file
        with open(args.output_file, "w") as output_file:
            output_file.write(comments)
            output_file.write(buffer.getvalue())


if __name__ == "__main__":
    main(create_parser().parse_args())
