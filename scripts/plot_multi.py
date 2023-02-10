import argparse
import logging
from dataclasses import dataclass
from os import path
from typing import Optional

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

logging.basicConfig(level=logging.INFO)


@dataclass(slots=True, frozen=True)
class InputDescriptor:
    slices: list[int]
    x0: int
    x1: Optional[int] = None

    def __len__(self) -> int:
        return len(self.slices)

    def is_x(self, axis: int) -> bool:
        if axis == self.x0:
            return True
        if self.x1 is not None:
            if axis == self.x1:
                return True
        return False

    @property
    def plot_dim(self) -> int:
        if self.x1 is not None:
            return 2
        return 1


def create_parser() -> argparse.ArgumentParser:
    def file_type(string: str) -> str:
        if path.isfile(string):
            return string
        else:
            raise FileNotFoundError(string)

    def descriptor_type(string: str) -> InputDescriptor:
        tokens = string.strip().replace(",", " ").replace(";", " ").split()
        x0 = None
        x1 = None
        slices = []
        for i, tok in enumerate(tokens):
            slices.append(0)
            if tok.lower() == "x0" and x0 is None:
                x0 = i
            elif tok.lower() == "x1" and x1 is None:
                x1 = i
            else:
                slices[-1] = int(tok)
        # Handle the case where x1 was set, but not x0
        if x0 is None and x1 is not None:
            x0 = x1
            x1 = None
        if x0 is None:
            raise ValueError("No 'x0' was provided.")

        return InputDescriptor(slices=slices, x0=x0, x1=x1)

    parser = argparse.ArgumentParser(description="Plot interpolation results.")

    parser.add_argument(
        "data_file",
        type=file_type,
        help="path to CSV file containing the interpolated data.",
    )

    parser.add_argument(
        "-i",
        "--input-descriptor",
        dest="input_descriptor",
        type=descriptor_type,
        help="List of indices used as the input from which the output data was generated.",
        required=True,
    )

    parser.add_argument(
        "-t",
        "--plot-target-dim",
        dest="target_dim",
        type=int,
        default=None,
        help="The index of the dimension used as the output dimension",
        required=True,
    )

    parser.add_argument(
        "-p",
        "--points",
        dest="points_file",
        type=file_type,
        help="path to CSV file containing the input points.",
    )
    return parser


def plot_2d(
    fig: mpl.figure.Figure,  # figure to plot on
    gspec: gridspec.SubplotSpec,  # subplot specification
    data: np.ndarray,  # data to plot, as a NumPy array
    extent: tuple[float, ...],  # extent of the plotted data
    points: Optional[np.ndarray] = None,  # points to plot, as NumPy array
) -> tuple[mpl.axes.Axes, mpl.cm.ScalarMappable]:
    """Plots an image of the given data on a subplot of the given figure with the data
    points overlayed as scatter plot.

    Args:
        fig: The figure to plot on.
        gspec: The subplot specification for the figure.
        points: The points to plot, as a NumPy array.
        data: The data to plot, as a NumPy array.
        extent: The extent of the plotted data, as a tuple of floats in the form (xmin,
            xmax, ymin, ymax).

    Returns:
        A tuple containing the Axes object for the subplot and a scalar mappable object
        for the plotted data.
    """
    # create a normalization object for the data
    norm = mpl.colors.Normalize(vmin=np.min(data), vmax=np.max(data))

    # create a colormap object
    # cmap = mpl.colormaps["jet"]
    cmap = mpl.colormaps["viridis"]

    # add a subplot to the figure and assign it to the variable ax
    ax = fig.add_subplot(gspec)

    # set the title and labels for the subplot
    ax.set_title("Predicted value")
    ax.set_xlabel("$x_0$")
    ax.set_ylabel("$x_1$")

    # plot the data as an image on the subplot
    ax.imshow(
        data.T,
        extent=extent,
        interpolation="none",
        origin="lower",
        aspect=(extent[1] - extent[0]) / (extent[3] - extent[2]),
        norm=norm,
        cmap=cmap,
    )

    if points is not None:
        ax.scatter(
            points[:, 0],
            points[:, 1],
            # 1,
            color="k",
            marker="^",
        )

    # return the subplot and the scalar mappable object
    return ax, mpl.cm.ScalarMappable(norm=norm, cmap=cmap)


def plot_3d(
    fig: mpl.figure.Figure,  # figure to plot on
    gspec: gridspec.SubplotSpec,  # subplot specification
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,  # data to plot, as a NumPy array
    # extent: tuple[float, ...],  # extent of the plotted data
    points: Optional[np.ndarray] = None,  # points to plot, as NumPy array
) -> tuple[mpl.axes.Axes, mpl.cm.ScalarMappable]:
    """Plots an image of the given data on a subplot of the given figure with the data
    points overlayed as scatter plot.

    Args:
        fig: The figure to plot on.
        gspec: The subplot specification for the figure.
        points: The points to plot, as a NumPy array.
        data: The data to plot, as a NumPy array.
        extent: The extent of the plotted data, as a tuple of floats in the form (xmin,
            xmax, ymin, ymax).

    Returns:
        A tuple containing the Axes object for the subplot and a scalar mappable object
        for the plotted data.
    """
    # create a normalization object for the data
    norm = mpl.colors.Normalize(vmin=np.min(z), vmax=np.max(z))

    # create a colormap object
    # cmap = mpl.colormaps["jet"]
    cmap = plt.get_cmap("viridis")

    # add a subplot to the figure and assign it to the variable ax
    ax = fig.add_subplot(gspec, projection="3d")

    # set the title and labels for the subplot
    ax.set_title("Predicted value")
    ax.set_xlabel("$x_0$")
    ax.set_ylabel("$x_1$")
    ax.set_zlabel("$target$")

    # plot the data as an image on the subplot
    ax.plot_surface(
        x,
        y,
        z,
        rstride=1,
        cstride=1,
        cmap=cmap,
        norm=norm,
    )

    if points is not None:
        ax.scatter3D(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            # 1,
            color="k",
            marker="^",
        )

    # return the subplot and the scalar mappable object
    return ax, mpl.cm.ScalarMappable(norm=norm, cmap=cmap)


def main(args: argparse.Namespace) -> None:
    data = np.loadtxt(args.data_file, delimiter=",")

    points = None
    if args.points_file:
        points = np.loadtxt(args.points_file, delimiter=",")

    input_descriptor = args.input_descriptor

    if len(input_descriptor) >= data.shape[1]:
        logging.error(
            "The provided data has less columns than required by the input descriptor"
        )
        return

    if points is not None:
        if points.shape[1] < len(input_descriptor):
            logging.error(
                "The provided points file does not cover enough dimensions for the provided input descriptor."
            )
            return

    # Gather all unique values in the input dimensions
    unique_inputs = [np.unique(data[:, i]) for i in range(len(input_descriptor))]

    # Create a mask
    mask = np.full(data.shape[0], fill_value=True)

    points_mask = None
    if points is not None:
        points_mask = np.full(points.shape[0], fill_value=True)

    for i, (s, u) in enumerate(zip(input_descriptor.slices, unique_inputs)):
        if input_descriptor.is_x(i):
            # Skip the line corresponding to x0 or x1
            continue

        if s >= len(u):
            logging.error(
                f"The provided input slice at position {i} is greater than or equal to number of unique values in that dimension."
            )
            return
        mask = mask & (data[:, i] == u[s])
        if points is not None:
            points_mask = points_mask & (points[:, i] == u[s])

    target_dim = args.target_dim
    if 0 <= target_dim < len(input_descriptor.slices):
        logging.warning(
            "The provided target dimension lies within the range of input dimensions."
        )
    elif target_dim >= data.shape[1]:
        logging.error(
            "The provided target dimension is larger than the number of columns in the provided data file."
        )
        return

    if input_descriptor.plot_dim == 1:
        logging.info("Line Plot")
    elif input_descriptor.plot_dim == 2:
        logging.info(input_descriptor)

        # Infer plot extent and original 2d data shape
        extent_list: list[float] = []
        shape_list: list[int] = []
        for i in [input_descriptor.x0, input_descriptor.x1]:
            shape_list.append(len(unique_inputs[i]))
            extent_list.extend([np.min(data[mask, i]), np.max(data[mask, i])])
        extent = tuple(extent_list)
        shape = tuple(shape_list)

        dat = data[mask, target_dim]
        if input_descriptor.x1 < input_descriptor.x0:
            extent = (*extent[2:], *extent[:2])
            shape = (shape[1], shape[0])
            dat = dat.reshape(shape).T
        else:
            dat = dat.reshape(shape)

        if points is not None:
            points = points[points_mask][
                :, [input_descriptor.x0, input_descriptor.x1, target_dim]
            ]

        fig = plt.figure()

        # Create a GridSpec object for the figure
        gspec = gridspec.GridSpec(1, 8, figure=fig)

        x, y = np.meshgrid(
            unique_inputs[input_descriptor.x0], unique_inputs[input_descriptor.x1]
        )
        _, sm = plot_3d(
            fig=fig,
            gspec=gspec[0, :-1],
            x=x,
            y=y,
            z=dat.T,
            points=points,
        )

        # Plot the value data on the first subplot
        # _, sm = plot_2d(
        #     fig=fig,
        #     gspec=gspec[0, :-1],
        #     data=dat,
        #     extent=extent,
        #     points=points,
        # )

        # Add a colorbar for the value data
        fig.colorbar(
            sm,
            cax=fig.add_subplot(gspec[0, -1]),
            orientation="vertical",
            label="Predicted Value",
        )
        plt.tight_layout()
        plt.show()
    else:
        logging.error("Currently only 1D and 2D plots are supported.")


if __name__ == "__main__":
    main(create_parser().parse_args())
