from manim import (
    BLACK,
    BLUE,
    DOWN,
    RED,
    RIGHT,
    WHITE,
    Axes,
    Brace,
    BraceBetweenPoints,
    Create,
    DrawBorderThenFill,
    FadeIn,
    FadeOut,
    Line,
    ParametricFunction,
    ReplacementTransform,
    Scene,
    Tex,
    Text,
    Write,
)
from numpy import sqrt


class Baseline(Scene):
    def construct(self):
        self.next_section("draw_bars")

        self.camera.background_color = WHITE  # type: ignore
        Text.set_default(color=BLACK)
        Line.set_default(color=BLACK)
        Tex.set_default(color=BLACK)
        Brace.set_default(color=BLACK)
        ParametricFunction.set_default(color=BLACK)

        # Create plane
        ax = Axes(x_range=[0, 1, 0.05], y_range=[0, 1, 0.05], x_length=10, y_length=10)
        ax_label = ax.get_axis_labels(x_label=r"t")

        self.play(Create(ax), Write(ax_label))

        f = ax.plot(lambda x: sqrt(x))
        area = ax.get_area(f, (0, 1), opacity=0.5, color=BLUE)  # type: ignore

        # Discrete number of firms
        rects = ax.get_riemann_rectangles(
            f, (0.2, 1), dx=0.2, input_sample_type="right"
        )
        self.play(DrawBorderThenFill(rects))

        self.next_section("add_brace")

        # Brace
        brace = BraceBetweenPoints(ax.c2p(0.2, 0), ax.c2p(0.4, 0))
        label = Tex(r"$\Pr[P \text{ is second}] = \frac{1}{5}$")
        label.next_to(brace, DOWN)
        self.play(Create(brace), Write(label))

        self.next_section("more_firms")

        # More firms
        rects2 = ax.get_riemann_rectangles(
            f, (0.05, 1), dx=0.05, input_sample_type="right"
        )

        # Brace again
        brace2 = BraceBetweenPoints(ax.c2p(0.05, 0), ax.c2p(0.1, 0))
        label2 = Tex(r"$\Pr[P \text{ is second}] = \frac{1}{20}$")
        label2.next_to(brace2, DOWN)

        self.play(
            ReplacementTransform(rects, rects2),
            ReplacementTransform(brace, brace2),
            ReplacementTransform(label, label2),
        )

        self.next_section("continuous_approximation")

        # Continuous approximation

        self.play(FadeOut(brace), FadeOut(label), FadeOut(brace2), FadeOut(label2))
        f_label = ax.get_graph_label(f, "f(t)", direction=RIGHT)  # type: ignore
        self.play(Create(f), Create(f_label))
        self.play(FadeOut(rects2), FadeIn(area))

        self.next_section("surplus_division")

        # Surplus division

        top = ax.plot(lambda x: 1)
        self.play(Create(top))

        area_fringe = ax.get_area(f, (0, 1), bounded_graph=top, opacity=0.5, color=RED)  # type: ignore
        self.play(FadeIn(area_fringe))
