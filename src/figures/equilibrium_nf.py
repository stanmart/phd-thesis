import matplotlib.pyplot as plt
import numpy as np
import typer


def plot_equilibrium_nf(N_P_list, mu, V_P, V_F, I_F, N_F_range, num_obs=500):
    N_F = np.linspace(N_F_range[0], N_F_range[1], num_obs)
    rhs = I_F * N_F
    lhs_list = []
    for N_P in N_P_list:
        lhs = mu * (
            (N_P * V_P + N_F * V_F) / (N_P * V_P + N_F * V_F + 1)
            - 1
            - (np.log((N_P * V_P + 1) / (N_P * V_P + N_F * V_F + 1))) / (N_F * V_F)
        )
        lhs_list.append(lhs)

    fig, ax = plt.subplots()
    ax.plot(N_F, rhs, label=r"$N_F$", color="black")
    for lhs, N_P in zip(lhs_list, N_P_list):
        ax.plot(N_F, lhs, label=r"$\phi_F(N_F) \mid N_P = {}$".format(N_P))

    ax.set_xlabel(r"$\bar N_F$")
    ax.legend()

    return fig, ax


def create_plot(
    out_path: str,
    mu: float = 1,
    V_P: float = 1,
    V_F: float = 1,
    I_F: float = 0.2,
    N_P: list[float] = [],
    N_F_range: tuple[float, float] = (0, 1),
    num_obs: int = 500,
    width: float = 5,
    height: float = 3.5,
    dpi: int = 300,
):
    if not N_P:
        raise typer.BadParameter("At least one value for N_P must be provided")

    fig, ax = plot_equilibrium_nf(
        N_P_list=N_P,
        mu=mu,
        V_P=V_P,
        V_F=V_F,
        I_F=I_F,
        N_F_range=N_F_range,
        num_obs=num_obs,
    )

    fig.set_size_inches(width, height)
    fig.set_facecolor("white")

    fig.savefig(out_path, dpi=dpi)


if __name__ == "__main__":
    typer.run(create_plot)
