from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

plt.rcParams.update({"text.usetex": True})


def two_sided_value(
    fun: Callable[[NDArray | float, NDArray | float], NDArray],
    num_obs_per_axis: int = 100,
    lambda_2: float = 1,
):
    a1 = np.linspace(0, 1, num_obs_per_axis)
    a2 = np.linspace(0, 1, num_obs_per_axis)
    a1_grid, a2_grid = np.meshgrid(a1, a2)
    z = fun(a1_grid, a2_grid)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_proj_type("ortho")

    ax.plot_surface(
        a1_grid,
        a2_grid,
        z,
        alpha=0.5,
        rstride=10,
        cstride=10,
        edgecolor="k",
        linewidth=0.5,
    )

    a1_diag = a1
    a2_diag = a1**lambda_2

    ax.plot(a1_diag, a2_diag, np.zeros_like(a1_diag), color="black", linewidth=2)
    ax.plot(
        a1_diag,
        a2_diag,
        fun(a1_diag, a2_diag),
        color="black",
        linewidth=2,
    )

    # Fill the area between the diagonal and the function along the diagonal
    for x, y in zip(a1_diag, a2_diag):
        ax.plot(
            [x, x],
            [y, y],
            [0, fun(x, y)],
            color="blue",
            alpha=0.5,
        )

    ax.set_xlabel("$a_1$")
    ax.set_ylabel("$a_2$")

    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_zlim([0, fun(1, 1)])

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_zticks([0, fun(1, 1)])  # type: ignore
    ax.set_zticklabels([0, "$f(1, 1)$"])  # type: ignore

    ax.view_init(elev=20, azim=-70)

    return fig, ax


if __name__ == "__main__":
    out_path = snakemake.output[0]  # type: ignore  # noqa: F821
    lambda_2 = float(snakemake.wildcards["lambda_2"])  # type: ignore  # noqa: F821

    fig, _ = two_sided_value(
        fun=lambda a1, a2: a1 * a2,  # type: ignore
        num_obs_per_axis=100,
        lambda_2=lambda_2,
    )

    fig.set_size_inches(3, 3)
    fig.set_facecolor("white")
    fig.savefig(out_path, bbox_inches="tight")
