from manim import (
    BLACK,
    BLUE_D,
    DOWN,
    RED_D,
    UP,
    WHITE,
    Axes,
    Brace,
    BraceBetweenPoints,
    Create,
    Line,
    MathTex,
    ParametricFunction,
    Scene,
    Tex,
    Text,
    Transform,
    Write,
)
from numpy import log
from scipy.optimize import fsolve


def pi_F(N_P, N_F):
    numerator = N_F * (N_F + N_P) - (N_F + N_P + 1) * (
        N_F + log(N_P + 1) - log(N_F + N_P + 1)
    )
    denominator = N_F * (N_F + N_P + 1)
    return numerator / denominator


class Baseline(Scene):
    def construct(self):
        self.next_section("draw_graph")

        self.camera.background_color = WHITE  # type: ignore
        Text.set_default(color=BLACK)
        Line.set_default(color=BLACK)
        Tex.set_default(color=BLACK)
        MathTex.set_default(color=BLACK)
        Brace.set_default(color=BLACK)
        ParametricFunction.set_default(color=BLACK)

        # Create plane
        ax = Axes(
            x_range=[0.005, 30, 0.1],
            y_range=[0, 0.25, 0.01],
            x_length=12,
            y_length=8,
            x_axis_config={"include_ticks": False},
            y_axis_config={"include_ticks": False},
        )
        x_label = ax.get_x_axis_label(r"N_F")
        y_label = ax.get_y_axis_label(r"\pi_F")

        N_P_0 = 0
        N_P_1 = 0.5
        I_F = 0.005

        # Intersections
        N_F_opt_0 = fsolve(lambda N_F: pi_F(N_P_0, N_F) - I_F * N_F, 20)[0]
        N_F_opt_0_val = pi_F(N_P_0, N_F_opt_0)
        N_F_opt_point_0 = ax.c2p(N_F_opt_0, N_F_opt_0_val)
        N_F_opt_1 = fsolve(lambda N_F: pi_F(N_P_1, N_F) - I_F * N_F, 20)[0]
        N_F_opt_1_val = pi_F(N_P_1, N_F_opt_1)
        N_F_opt_point_1 = ax.c2p(N_F_opt_1, N_F_opt_1_val)

        # Destinations
        pi_F_orig = ax.plot(lambda x: pi_F(N_P_0, x), color=BLUE_D)
        pi_F_alt = ax.plot(lambda x: pi_F(N_P_1, x), color=RED_D)
        pi_F_orig_label = ax.get_graph_label(
            pi_F_orig,
            r"\pi_F(N_P, N_F)",
            direction=UP,  # type: ignore
        )
        pi_F_alt_label = ax.get_graph_label(
            pi_F_alt,
            r"\pi_F(N_P', N_F)",
            direction=DOWN,  # type: ignore
        )

        N_F_opt_0_bar = ax.get_vertical_line(N_F_opt_point_0, color=BLACK)
        N_F_opt_1_bar = ax.get_vertical_line(N_F_opt_point_1, color=BLACK)
        # N_F_opt_0_label = MathTex(r"N_F^*(N_P)", color=BLUE_D)
        # N_F_opt_0_label.next_to(N_F_opt_0_bar, DOWN)
        # N_F_opt_1_label = MathTex(r"N_F^*(N_P')", color=RED_D)
        # N_F_opt_1_label.next_to(N_F_opt_1_bar, DOWN)

        investment_cost = ax.plot(lambda x: I_F * x)
        investment_cost_label = ax.get_graph_label(
            investment_cost,
            r"I_F N_F",
            direction=UP,  # type: ignore
        )
        brace_loss = BraceBetweenPoints(
            N_F_opt_1_bar.get_bottom(),  # type: ignore
            N_F_opt_0_bar.get_bottom(),  # type: ignore
        )
        brace_label = MathTex(r"> N_P' - N_P", color=BLACK)
        brace_label.next_to(brace_loss, DOWN)

        # Moving objects
        pi_F_plot = pi_F_orig.copy()
        pi_F_label = pi_F_orig_label.copy()
        N_F_opt_bar = N_F_opt_0_bar.copy()
        # N_F_opt_label = N_F_opt_0_label.copy()

        # First phase: Draw pi_F and investment cost
        self.play(Create(ax), Write(x_label), Write(y_label))
        self.play(Create(pi_F_plot), Create(pi_F_label))
        self.play(Create(investment_cost), Create(investment_cost_label))

        # Second phase: mark equilibrium
        self.next_section("mark_equilibrium")
        self.play(Create(N_F_opt_bar))

        # Third phase: move to alternative equilibrium
        self.next_section("alternate_equilibrium")
        self.add(pi_F_orig, pi_F_orig_label, N_F_opt_0_bar)
        self.play(Transform(pi_F_plot, pi_F_alt), Transform(pi_F_label, pi_F_alt_label))
        self.wait(0.5)
        self.play(
            Transform(N_F_opt_bar, N_F_opt_1_bar),
            # Transform(N_F_opt_label, N_F_opt_1_label)
        )

        # Fourth phase: show loss
        self.next_section("show_loss")
        self.play(Create(brace_loss), Write(brace_label))
        self.wait(0.1)
