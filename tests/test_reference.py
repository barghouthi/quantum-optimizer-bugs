import importlib.util
import io
import random
import sys
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "bugs", Path(__file__).parents[1] / "bugs.py"
)
bugs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = bugs
spec.loader.exec_module(bugs)


class ReferenceTests(unittest.TestCase):
    def test_echo_phase_and_even_x_requirement(self):
        d = 3e-6
        echo = [("p", d), ("p", d), ("x", 0), ("p", -d), ("p", -d), ("x", 0)]
        self.assertEqual(bugs.phase(echo * 100), bugs.Phase(400 * Fraction(d)))
        with self.assertRaises(ValueError):
            bugs.phase([("x", 0)])

    def test_rational_sum_retains_large_angle_residual(self):
        self.assertEqual(
            bugs.phase([("rz", x) for x in (1e16, 1, 1, 1, -1e16)]),
            bugs.Phase(Fraction(3)),
        )
        self.assertEqual(
            bugs.phase([("rzpi", x) for x in (1e16, 1, -1e16)]),
            bugs.Phase(pi=Fraction(1)),
        )

    def test_phase_against_independent_dense_matrix_products(self):
        ctx = bugs.mpmath.mp.clone()
        ctx.dps = 100
        rng = random.Random(42)
        for _ in range(30):
            gates = [
                (
                    rng.choice(("p", "rz", "rzpi", "x")),
                    Fraction(rng.randrange(-10, 11), 7),
                )
                for _ in range(20)
            ]
            if sum(name == "x" for name, _ in gates) % 2:
                gates.append(("x", 0))
            matrix = ctx.eye(2)
            for name, angle in gates:
                a = (
                    ctx.mpf(angle.numerator) / angle.denominator
                    if isinstance(angle, Fraction)
                    else ctx.mpf(angle)
                )
                if name == "x":
                    gate = ctx.matrix([[0, 1], [1, 0]])
                elif name == "p":
                    gate = ctx.diag([1, ctx.exp(1j * a)])
                else:
                    if name == "rzpi":
                        a *= ctx.pi
                    gate = ctx.diag([ctx.exp(-1j * a / 2), ctx.exp(1j * a / 2)])
                matrix = gate * matrix
            self.assertLess(abs(matrix[0, 1]) + abs(matrix[1, 0]), ctx.mpf("1e-90"))
            self.assertLess(
                abs(
                    matrix[1, 1] / matrix[0, 0]
                    - ctx.exp(1j * bugs.phase(gates).value(ctx))
                ),
                ctx.mpf("1e-90"),
            )

    def test_distance_global_phase_periodicity_and_known_values(self):
        ctx = bugs.mpmath.mp.clone()
        ctx.dps = 100
        identity = bugs.Phase()
        self.assertLess(
            bugs.distance(ctx, identity, bugs.Phase(pi=Fraction(2))), ctx.mpf("1e-90")
        )
        self.assertLess(
            abs(bugs.distance(ctx, identity, bugs.Phase(pi=Fraction(1))) - ctx.sqrt(2)),
            ctx.mpf("1e-90"),
        )
        tiny = Fraction(1, 10**50)
        error = bugs.distance(ctx, identity, bugs.Phase(tiny))
        self.assertGreater(error, 0)
        self.assertLess(abs(error - ctx.mpf("5e-51")), ctx.mpf("1e-90"))

    def test_unsupported_and_nonfinite_angles_fail(self):
        for gates in ([("h", 0)], [("rz", float("inf"))], [("rz", float("nan"))]):
            with self.assertRaises((ValueError, OverflowError)):
                bugs.phase(gates)

    def test_distance_matches_optimal_matrix_phase_and_beats_sampled_phases(self):
        ctx = bugs.mpmath.mp.clone()
        ctx.dps = 100
        for angle in (Fraction(1, 3), Fraction(-7, 3), Fraction(10), Fraction(10**16)):
            theta = bugs.Phase(angle).value(ctx)
            delta = ctx.atan2(ctx.sin(theta), ctx.cos(theta))
            error = bugs.distance(ctx, bugs.Phase(angle), bugs.Phase())

            # Independently compute the largest singular value of the diagonal
            # matrix U - exp(i*phi) I; for diagonal matrices it is max(abs(diag)).
            def matrix_error(phi, theta=theta):
                return max(
                    abs(1 - ctx.exp(1j * phi)),
                    abs(ctx.exp(1j * theta) - ctx.exp(1j * phi)),
                )

            self.assertLess(abs(matrix_error(delta / 2) - error), ctx.mpf("1e-80"))
            for k in range(64):
                self.assertGreaterEqual(
                    matrix_error(k * ctx.pi / 32) + ctx.mpf("1e-80"), error
                )

    def test_report_preserves_qasm_count_when_import_erases_the_gate(self):
        stdout = io.StringIO()
        with (
            patch.object(
                bugs, "run_native", return_value=([], [], "pass", "loader", "note")
            ),
            patch.object(bugs.importlib.metadata, "version", return_value="test"),
            patch("sys.stdout", stdout),
        ):
            self.assertFalse(bugs.main(["--tools", "quizx"]))
        self.assertIn("1 / 0 / 0 gates", stdout.getvalue())

    def test_measurement_witness_and_stage_attribution(self):
        ideal = bugs.Phase(pi=Fraction(1))
        for imported, expected_import, expected_optimization in (
            (ideal, "0.0", False),
            (bugs.Phase(), None, True),
        ):
            result = bugs.metrics(ideal, imported, bugs.Phase())
            self.assertEqual(result["measure"], "X")
            self.assertEqual(float(result["probability_change"]), 1)
            if expected_import:
                self.assertEqual(result["import_error"], expected_import)
            if expected_optimization:
                self.assertEqual(result["optimization_error"], "0.0")

    def test_literal_roundtrips_and_fixture_scope(self):
        from decimal import Decimal

        for value in (3e-6, 1e-4):
            self.assertEqual(Fraction(Decimal(bugs.literal(value))), Fraction(value))
        self.assertEqual(len(bugs.TOOLS), 6)
        for tool in bugs.TOOLS:
            source, _, _, _ = bugs.fixture(tool)
            self.assertTrue(source.startswith("OPENQASM 2.0;"))
            self.assertNotIn("1e308", source)
            self.assertNotIn("inf", source)


if __name__ == "__main__":
    unittest.main()
