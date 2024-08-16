from collections.abc import Iterable

import numpy as np
import pandas as pd
import typer
from scipy.integrate import quad
from scipy.optimize import fsolve
from sympy import Symbol, lambdify, log  # type: ignore


def create_plot_data(
    out_path: str,
    mu: float = 1,
    V_P: float = 1,
    V_F: float = 1,
    I_F: float = 0.2,
    N_C: float = 1,
    N_P_range: tuple[float, float] = (0, 1),
    num_obs: int = 500,
    lambda_P: float = 1,
    value_function: str = "profit",
    bargaining: str = "onesided",
):
    N_P_sp, N_F_sp = Symbol("N_P"), Symbol("N_F")
    s = Symbol("s")

    f_profit = (
        N_C
        * mu
        * (s * N_F_sp * V_F + N_P_sp * V_P)  # type: ignore
        / (s * N_F_sp * V_F + N_P_sp * V_P + 1)  # type: ignore
    )
    f_surplus = N_C * mu * log(s * N_F_sp * V_F + N_P_sp * V_P + 1)  # type: ignore

    if value_function == "profit":
        f = f_profit
    elif value_function == "welfare":
        f = f_profit + f_surplus
    else:
        raise typer.BadParameter(f"Unknown value function: {value_function}")

    f_num = lambdify((N_P_sp, N_F_sp, s), f)
    f_diff_num = lambdify((N_P_sp, N_F_sp, s), f.diff(s))

    if bargaining == "onesided":

        def pi_P_scalar(N_P, N_F):
            return quad(
                lambda s: f_num(N_P, N_F, s) * lambda_P * s ** (lambda_P - 1), 0, 1
            )[0]

        def pi_F_scalar(N_P, N_F):
            return quad(lambda s: f_diff_num(N_P, N_F, s) * s**lambda_P, 0, 1)[0]
    elif bargaining == "twosided":

        def pi_P_scalar(N_P, N_F):
            return quad(lambda s: f_num(N_P, N_F, s) * lambda_P * s**lambda_P, 0, 1)[0]

        def pi_F_scalar(N_P, N_F):
            return quad(lambda s: f_diff_num(N_P, N_F, s) * s ** (lambda_P + 1), 0, 1)[
                0
            ]
    else:
        raise typer.BadParameter(f"Unknown bargaining mode: {bargaining}")

    pi_P = np.vectorize(pi_P_scalar)
    pi_F = np.vectorize(pi_F_scalar)

    def pi_F_t(N_P, N_F):
        return pi_F(N_P, N_F) - I_F * N_F

    def N_F_opt(N_P):
        if isinstance(N_P, Iterable):
            return np.array([N_F_opt(n) for n in N_P])
        else:
            return fsolve(lambda N_F: pi_F_t(N_P, N_F), 3)[0]

    def K_F_opt():
        return np.sqrt(mu * I_F * V_F) - I_F

    def N_F_opt_bench(N_P):
        K_F = K_F_opt()
        N_F_candidate = mu / (K_F + I_F) - N_P * V_P / V_F - 1 / V_F
        return np.maximum(N_F_candidate, 0)

    N_P_vec = np.linspace(N_P_range[0], N_P_range[1], num_obs)
    N_F_vec = N_F_opt(N_P_vec)

    pi_F_vec = pi_F(N_P_vec, N_F_vec)
    pi_P_vec = pi_P(N_P_vec, N_F_vec)
    A_vec = N_P_vec * V_P + N_F_vec * V_F + 1
    CS_vec = np.log(A_vec)

    N_F_cf = np.maximum(N_F_vec[0] * V_F - N_P_vec * V_P, 0)
    A_cf = N_P_vec * V_P + N_F_cf * V_F + 1
    CS_cf = np.log(A_cf)

    A_noF = N_P_vec * V_P + 1
    CS_noF = np.log(A_noF)
    pi_P_noF = mu * N_P_vec * V_P / A_noF

    pi_P_var_vec = mu * N_P_vec * V_P / A_vec

    if bargaining == "onesided" and value_function == "profit":
        N_F_bench = N_F_opt_bench(N_P_vec)
        A_bench = N_P_vec * V_P + N_F_bench * V_F + 1
        CS_bench = np.log(A_bench)
        pi_P_bench = K_F_opt() * N_F_bench + mu * N_P_vec * V_P / A_bench
        K_F_opt_vec = np.where(N_F_bench > 1e-5, K_F_opt(), np.nan)
        hybrid_mode_bench = np.where(N_F_bench > 1e-5, 10, 0)
    else:
        N_F_bench = np.ones_like(N_P_vec, dtype=float) * np.nan
        A_bench = np.ones_like(N_P_vec, dtype=float) * np.nan
        CS_bench = np.ones_like(N_P_vec, dtype=float) * np.nan
        pi_P_bench = np.ones_like(N_P_vec, dtype=float) * np.nan
        K_F_opt_vec = np.ones_like(N_P_vec, dtype=float) * np.nan
        hybrid_mode_bench = np.ones_like(N_P_vec, dtype=float) * np.nan

    K_F_implied = np.where(N_F_vec > 1e-5, (pi_P_vec - pi_P_var_vec) / N_F_vec, np.nan)

    hybrid_mode = np.where(N_F_vec > 1e-5, 10, 0)

    A_noF = N_P_vec * V_P + 1
    CS_noF = np.log(A_noF)
    pi_P_noF = mu * N_P_vec * V_P / A_noF

    data = pd.DataFrame(
        {
            "N_P": N_P_vec,
            "N_F": N_F_vec,
            "pi_F": pi_F_vec,
            "pi_P": pi_P_vec,
            "A": A_vec,
            "CS": CS_vec,
            "N_F_cf": N_F_cf,
            "A_cf": A_cf,
            "CS_cf": CS_cf,
            "A_noF": A_noF,
            "CS_noF": CS_noF,
            "pi_P_noF": pi_P_noF,
            "K_F_implied": K_F_implied,
            "K_F_opt": K_F_opt_vec,
            "hybrid": hybrid_mode,
            "N_F_bench": N_F_bench,
            "A_bench": A_bench,
            "CS_bench": CS_bench,
            "pi_P_bench": pi_P_bench,
            "hybrid_bench": hybrid_mode_bench,
        }
    )

    data.to_csv(out_path, index=False)


if __name__ == "__main__":
    typer.run(create_plot_data)
