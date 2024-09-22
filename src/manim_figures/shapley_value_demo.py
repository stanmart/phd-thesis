from manim import (
    BLACK,
    BLUE_D,
    DOWN,
    LEFT,
    RED,
    RIGHT,
    UP,
    WHITE,
    Brace,
    BraceBetweenPoints,
    FadeIn,
    Line,
    MoveToTarget,
    Rectangle,
    ReplacementTransform,
    Scene,
    Tex,
    Text,
    Transform,
    VGroup,
    Write,
)


class ShapleyDemo(Scene):
    def construct(self):
        box_size = 2
        label_size = 150
        horiz_space = 7
        up_shift = 3

        self.camera.background_color = WHITE  # type: ignore

        Text.set_default(color=BLACK)
        Line.set_default(color=BLACK)
        Tex.set_default(color=BLACK)
        Brace.set_default(color=BLACK)

        self.next_section("draw_objects")

        P_text = Tex("$P$", color=BLACK, font_size=label_size)
        P_text.bg = Rectangle(  # type: ignore
            height=box_size,
            width=box_size,
            color=BLUE_D,  # type: ignore
            fill_opacity=0.25,
        )
        P = VGroup(P_text, P_text.bg)

        A1_text = Tex("$A_1$", color=BLACK, font_size=label_size)
        A1_text.bg = Rectangle(  # type: ignore
            height=box_size,
            width=box_size,
            color=RED,  # type: ignore
            fill_opacity=0.25,
        )
        A1 = VGroup(A1_text, A1_text.bg)

        A2_text = Tex("$A_2$", color=BLACK, font_size=label_size)
        A2_text.bg = Rectangle(  # type: ignore
            height=box_size,
            width=box_size,
            color=RED,  # type: ignore
            fill_opacity=0.25,
        )
        A2 = VGroup(A2_text, A2_text.bg)

        P.next_to(A1, horiz_space * LEFT)
        A2.next_to(A1, horiz_space * RIGHT)

        P.shift(up_shift * UP)
        A1.shift(up_shift * UP)
        A2.shift(up_shift * UP)

        eq = Tex(
            "$\\varphi_P$ = ",
            "$\\frac{1}{6}$",
            "$\\bigg($",
            "$0$",
            "$+0$",
            "$+2$",
            "$+2$",
            "$+4$",
            "$+4$",
            "$\\bigg)$" "$=2$",
        )
        eq.next_to(A1, up_shift * 6 * DOWN)
        eq.font_size = 80

        self.play(FadeIn(P, A1, A2))

        left = P.get_center()
        middle = A1.get_center()
        right = A2.get_center()

        pl0 = P.get_corner(DOWN + LEFT)
        pl1 = P.get_corner(DOWN + RIGHT)
        pm1 = A1.get_corner(DOWN + RIGHT)
        pr1 = A2.get_corner(DOWN + LEFT)
        pr1 = A2.get_corner(DOWN + RIGHT)

        self.next_section("value_P12")

        brace_1 = BraceBetweenPoints(pl0, pl1)  # type: ignore
        brace_1_text = brace_1.get_text("0")
        brace_1_text.font_size = 80

        self.play(FadeIn(brace_1), FadeIn(brace_1_text))
        self.play(Write(eq[0]), Write(eq[3]))

        self.next_section("value_P21")

        A1.generate_target()
        A1.target.shift(right - middle)  # type: ignore
        A2.generate_target()
        A2.target.shift(middle - right)  # type: ignore

        self.play(MoveToTarget(A1), MoveToTarget(A2))
        self.play(Write(eq[4]))

        self.next_section("value_2P1")

        P.generate_target()
        P.target.shift(middle - left)  # type: ignore
        A2.generate_target()
        A2.target.shift(left - middle)  # type: ignore

        brace_2 = BraceBetweenPoints(pl0, pm1)  # type: ignore
        brace_2.shift(DOWN)
        brace_2_text = brace_2.get_text("3")
        brace_2_text.font_size = 80

        brace_1_text_new = brace_1.get_text("1")
        brace_1_text_new.font_size = 80

        self.play(MoveToTarget(P), MoveToTarget(A2))
        self.play(
            ReplacementTransform(brace_1_text, brace_1_text_new),
            FadeIn(brace_2),
            FadeIn(brace_2_text),
        )
        self.play(Write(eq[5]))

        self.next_section("value_1P2")

        A1.generate_target()
        A1.target.shift(left - right)  # type: ignore
        A2.generate_target()
        A2.target.shift(right - left)  # type: ignore

        self.play(MoveToTarget(A1), MoveToTarget(A2))
        self.play(Write(eq[6]))

        self.next_section("value_12P")

        P.generate_target()
        P.target.shift(right - middle)  # type: ignore
        A2.generate_target()
        A2.target.shift(middle - right)  # type: ignore

        brace_2_up = BraceBetweenPoints(pl0, pm1)  # type: ignore
        brace_2_up_text = brace_2_up.get_text("2")
        brace_2_up_text.font_size = 80

        brace_3 = BraceBetweenPoints(pl0, pr1)  # type: ignore
        brace_3.shift(DOWN)
        brace_3_text = brace_3.get_text("6")
        brace_3_text.font_size = 80

        self.play(
            MoveToTarget(P),
            MoveToTarget(A2),
            Transform(brace_1, brace_2_up),
            Transform(brace_2, brace_3),
            Transform(brace_1_text_new, brace_2_up_text),
            Transform(brace_2_text, brace_3_text),
        )
        self.play(Write(eq[7]))

        self.next_section("value_21P")

        A1.generate_target()
        A1.target.shift(middle - left)  # type: ignore
        A2.generate_target()
        A2.target.shift(left - middle)  # type: ignore

        self.play(MoveToTarget(A1), MoveToTarget(A2))
        self.play(Write(eq[8]))

        self.next_section("final")

        self.play(Write(eq[1]), Write(eq[2]), Write(eq[9]))
