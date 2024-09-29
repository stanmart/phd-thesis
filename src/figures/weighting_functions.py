from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"text.usetex": True})


def plot_weighting_functions(
    param_list: list[int | float], num_obs=500, type: str = "multiple"
):
    x = np.linspace(0, 1, num_obs)

    linestyles = cycle(["solid", "dashed", "dotted", "dashdot"])

    fig, ax = plt.subplots()
    for par, linestyle in zip(param_list, linestyles):
        if type == "multiple_p":
            y = par * (1 - x) ** (par - 1)
            label = f"$m={par}$"
        elif type == "weighted_shapley":
            y = par * x ** (par - 1)
            label = f"$\\lambda_P={par}$"
        else:
            raise ValueError(f"Unknown type: {type}")
        ax.plot(x, y, label=label, color="black", linestyle=linestyle)

    ax.set_xlabel("$h(x)$")
    ax.set_xticks([0, 1])
    if type == "multiple_p":
        ax.set_yticks([0] + param_list)

    ax.set_xlim(0, 1)
    if type == "multiple_p":
        ax.set_ylim(0, ax.get_ylim()[1])
    else:
        ax.set_ylim(0, 3)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend()

    return fig, ax


if __name__ == "__main__":
    param_list = [float(par) for par in snakemake.wildcards["par_list"].split(",")]  # type: ignore  # noqa F821
    weighting_type = snakemake.wildcards["type"]  # type: ignore  # noqa F821
    out_path = snakemake.output[0]  # type: ignore # noqa F821

    fig, _ = plot_weighting_functions(
        param_list=param_list,
        num_obs=500,
        type=weighting_type,
    )

    fig.set_size_inches(3, 3)
    fig.set_facecolor("white")
    fig.tight_layout()

    fig.savefig(out_path, dpi=300)
