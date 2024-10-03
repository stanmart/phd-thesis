import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

plt.rcParams.update({"text.usetex": True})


VARS = {
    "entry-fee": {
        "benchmark": "K_F_opt",
        "bargaining": "K_F_implied",
        "pure_retail": None,
        "y_label": "$K_F$",
    },
    "platform-profit": {
        "benchmark": "pi_P_bench",
        "bargaining": "pi_P",
        "pure_retail": "pi_P_noF",
        "y_label": "$pi_P$",
    },
    "fringe-number": {
        "benchmark": "N_F_bench",
        "bargaining": "N_F",
        "pure_retail": None,
        "y_label": "$N_F$",
    },
    "aggregate": {
        "benchmark": "A_bench",
        "bargaining": "A",
        "pure_retail": "A_noF",
        "y_label": "$A$",
    },
    "consumer-surplus": {
        "benchmark": "CS_bench",
        "bargaining": "CS",
        "pure_retail": "CS_noF",
        "y_label": "$CS$",
    },
}


def plot_equilibrium_outcomes(
    df: pd.DataFrame,
    benchmark_var: str,
    bargaining_var: str | None = None,
    pure_retail_var: str | None = None,
    y_label: str = "",
    hybrid_indicator_benchmark: str | None = None,
    hybrid_indicator_bargaining: str | None = None,
) -> tuple[Figure, Axes]:
    fig, ax = plt.subplots()

    ax.plot(df["N_P"], df[benchmark_var], label="Benchmark")

    if bargaining_var is not None:
        ax.plot(df["N_P"], df[bargaining_var], label="Bargaining")

    if pure_retail_var is not None:
        ax.plot(
            df["N_P"],
            df[pure_retail_var],
            label="Pure retail",
            color="black",
            linestyle=":",
        )

    if hybrid_indicator_benchmark is not None:
        shade_end = df.loc[df[hybrid_indicator_benchmark] != 0, "N_P"].max()
        ax.axvspan(0, shade_end, color="gray", alpha=0.2)

    if hybrid_indicator_bargaining is not None:
        shade_end = df.loc[df[hybrid_indicator_bargaining] != 0, "N_P"].max()
        ax.axvspan(0, shade_end, color="gray", alpha=0.2)

    ax.set_xlabel("$N_P$")
    ax.set_ylabel(y_label)

    ax.set_xlim(0, df["N_P"].max())
    ax.set_xticks([0, df["N_P"].max()])

    ax.spines[["right", "top"]].set_visible(False)
    ax.legend()

    return fig, ax


if __name__ == "__main__":
    input_data = snakemake.input.csv  # type: ignore  # noqa: F821
    output_figure = snakemake.output[0]  # type: ignore  # noqa: F821
    var = snakemake.wildcards.var  # type: ignore  # noqa: F821
    plot_bargaining = snakemake.wildcards.add_bargaining == "with"  # type: ignore  # noqa: F821
    # TODO: handle vars

    df = pd.read_csv(input_data)
    fig, _ = plot_equilibrium_outcomes(
        df=df,
        benchmark_var=VARS[var]["benchmark"],
        bargaining_var=VARS[var]["bargaining"] if plot_bargaining else None,
        pure_retail_var=VARS[var]["pure_retail"] if plot_bargaining else None,
        y_label=VARS[var]["y_label"],
        hybrid_indicator_benchmark="hybrid_bench",
        hybrid_indicator_bargaining="hybrid" if plot_bargaining else None,
    )

    fig.set_size_inches(3, 2.5)
    fig.tight_layout()
    fig.savefig(output_figure, bbox_inches="tight", dpi=300)
