from manim import (
    BLACK,
    BLUE,
    DOWN,
    RED,
    RIGHT,
    WHITE,
    Axes,
    Brace,
    Create,
    FadeIn,
    Line,
    MathTex,
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
        self.next_section("draw_graph")

        self.camera.background_color = WHITE  # type: ignore
        Text.set_default(color=BLACK)
        Line.set_default(color=BLACK)
        Tex.set_default(color=BLACK)
        MathTex.set_default(color=BLACK)
        Brace.set_default(color=BLACK)
        ParametricFunction.set_default(color=BLACK)

        # Create plane
        ax = Axes(x_range=[0, 1, 0.05], y_range=[0, 1, 0.05], x_length=10, y_length=10)
        ax_label = ax.get_x_axis_label(r"t")

        self.play(Create(ax), Write(ax_label))

        f = ax.plot(lambda x: x**2)
        f_label = ax.get_graph_label(f, "f(t)", direction=RIGHT)  # type: ignore
        self.play(Create(f), Create(f_label))
        self.wait(0.5)

        area = ax.get_area(f, (0, 1), opacity=0.5, color=BLUE)  # type: ignore
        top = ax.plot(lambda x: 1)
        area_fringe = ax.get_area(f, (0, 1), bounded_graph=top, opacity=0.5, color=RED)  # type: ignore

        self.play(Create(top))

        self.play(FadeIn(area), FadeIn(area_fringe))

        # Convex function
        self.next_section("complements")

        f2 = ax.plot(lambda x: sqrt(x))
        f2_label = ax.get_graph_label(f2, "f(t)", direction=RIGHT)  # type: ignore
        area2 = ax.get_area(f2, (0, 1), opacity=0.5, color=BLUE)  # type: ignore
        area_fringe2 = ax.get_area(
            f2,
            (0, 1),
            bounded_graph=top,
            opacity=0.5,
            color=RED,  # type: ignore
        )

        self.play(
            ReplacementTransform(f, f2),
            ReplacementTransform(f_label, f2_label),
            ReplacementTransform(area, area2),
            ReplacementTransform(area_fringe, area_fringe2),
        )

        # Platform has intrinsic value
        self.next_section("intrinsic_value")

        f3 = ax.plot(lambda x: sqrt(0.1 + 0.9 * x))
        area3 = ax.get_area(f3, (0, 1), opacity=0.5, color=BLUE)  # type: ignore
        area_fringe3 = ax.get_area(
            f3,
            (0, 1),
            bounded_graph=top,
            opacity=0.5,
            color=RED,  # type: ignore
        )

        self.play(
            ReplacementTransform(f2, f3),
            ReplacementTransform(area2, area3),
            ReplacementTransform(area_fringe2, area_fringe3),
        )

        # Outside option for firms
        self.next_section("outside_option")

        f4 = ax.plot(lambda x: sqrt(0.1 + 0.9 * x))
        fo4 = ax.plot(lambda x: sqrt(x) / 2)
        fo4_label = ax.get_graph_label(fo4, "f_0(t)", direction=DOWN)  # type: ignore
        area4 = ax.get_area(f4, (0, 1), bounded_graph=fo4, opacity=0.5, color=BLUE)  # type: ignore
        areao4 = ax.get_area(fo4, (0, 1), opacity=0.5, color=RED)  # type: ignore
        area_fringe4 = ax.get_area(
            f4,
            (0, 1),
            bounded_graph=top,
            opacity=0.5,
            color=RED,  # type: ignore
        )

        self.play(ReplacementTransform(f3, f4), Create(fo4), Create(fo4_label))

        self.play(
            ReplacementTransform(area3, area4),
            ReplacementTransform(area_fringe3, area_fringe4),
            FadeIn(areao4),
        )
        self.wait(1)
