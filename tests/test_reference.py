import importlib.util
import io
import math
import unittest
from decimal import Decimal
from fractions import Fraction
from unittest.mock import patch

import bugs


class MeasurementTests(unittest.TestCase):
    def test_identity_z_and_s_measurements(self):
        matrices = (
            ([[1, 0], [0, 1]], 1, 0.5),
            ([[1, 0], [0, -1]], 0, 0.5),
            ([[1, 0], [0, 1j]], 0.5, 1),
            ([[1, 0], [0, -1j]], 0.5, 0),
        )
        for entries, x, y in matrices:
            matrix = bugs.mp.matrix(entries)
            self.assertAlmostEqual(bugs.output_probability(matrix, "X"), x)
            self.assertAlmostEqual(bugs.output_probability(matrix, "Y"), y)
            # The physical probabilities ignore the circuit's global phase.
            matrix *= bugs.mp.exp(0.73j)
            self.assertAlmostEqual(bugs.output_probability(matrix, "X"), x)
            self.assertAlmostEqual(bugs.output_probability(matrix, "Y"), y)

    def test_ideal_pi_and_three_radian_residual(self):
        self.assertEqual(bugs.ideal_probability("pi", "X"), 0)
        expected = (1 + bugs.mp.cos(3)) / 2
        self.assertEqual(bugs.ideal_probability(Fraction(3), "X"), expected)
        self.assertAlmostEqual(float(100 * (1 - expected)), 99.499624830, places=8)

    def test_fixture_precision_stability(self):
        for _, angle, axis, _, _ in bugs.examples().values():
            with bugs.mp.workdps(100):
                before = bugs.ideal_probability(angle, axis)
            with bugs.mp.workdps(160):
                after = bugs.ideal_probability(angle, axis)
                self.assertLess(abs(after - before), bugs.mp.mpf("1e-70"))

    def test_quizx_reference_against_exact_sum_and_matrix_product(self):
        phases = (Fraction(2147483647, 4294967295), Fraction(2147483648, 4294967297))
        total = sum(phases)
        self.assertEqual(total, 1 - Fraction(2**32, 2**64 - 1))
        _, angle, axis, _, _ = bugs.examples()["quizx"]
        ctx = bugs.mp.mp.clone()
        ctx.dps = 160
        matrix = ctx.eye(2)
        for phase in phases:
            theta = ctx.pi * phase.numerator / phase.denominator
            matrix = (
                ctx.diag([ctx.exp(-1j * theta / 2), ctx.exp(1j * theta / 2)]) * matrix
            )
        expected = abs(sum(matrix[i, j] for i in range(2) for j in range(2)) / 2) ** 2
        reference = bugs.ideal_probability(angle, axis)
        self.assertGreater(reference, 0)
        self.assertLess(abs(reference / expected - 1), ctx.mpf("1e-70"))

    @unittest.skipUnless(importlib.util.find_spec("quizx"), "requires quizx")
    def test_quizx_import_is_exact_but_simplification_erases_rotation(self):
        import quizx

        source, angle, axis, _, _ = bugs.examples()["quizx"]
        graph = quizx.qasm(source)
        self.assertEqual(
            [graph.phase(v) for v in graph.vertices() if graph.type(v) == 1],
            [Fraction(2147483647, 4294967295), Fraction(2147483648, 4294967297)],
        )
        before = quizx.extract_circuit(graph)
        self.assertLess(
            abs(
                bugs.output_probability(before.to_matrix(), axis)
                - bugs.ideal_probability(angle, axis)
            ),
            bugs.mp.mpf("1e-24"),
        )
        quizx.full_simp(graph)
        self.assertEqual(len(list(graph.vertices())), 2)
        self.assertTrue(all(graph.type(v) == 0 for v in graph.vertices()))
        after = quizx.extract_circuit(graph)
        self.assertEqual(len(after.gates), 0)
        self.assertEqual(bugs.output_probability(after.to_matrix(), axis), 1)

    def test_echo_reference_against_independent_matrix_product(self):
        ctx = bugs.mp.mp.clone()
        ctx.dps = 100
        d = Fraction(3e-6)
        a = ctx.mpf(d.numerator) / d.denominator
        positive = ctx.diag([1, ctx.exp(1j * a)])
        negative = ctx.diag([1, ctx.exp(-1j * a)])
        x = ctx.matrix([[0, 1], [1, 0]])
        matrix = (x * negative * negative * x * positive * positive) ** 100_000
        self.assertAlmostEqual(
            bugs.output_probability(matrix, "Y"),
            float(bugs.ideal_probability(d * 400_000, "Y")),
            places=14,
        )

    def test_exact_literals_and_six_qasm_inputs(self):
        for value in (3e-6, 1e-4):
            self.assertEqual(Fraction(Decimal(bugs.literal(value))), Fraction(value))
        self.assertEqual(tuple(bugs.examples()), bugs.TOOLS)
        for source, _, _, _, _ in bugs.examples().values():
            self.assertTrue(source.startswith("OPENQASM 2.0;"))
            self.assertNotIn("1e308", source)

    def test_nonfinite_output_fails(self):
        with self.assertRaises(ValueError):
            bugs.output_probability(bugs.mp.matrix([[math.inf, 0], [0, 1]]), "X")

    def test_report_uses_actual_output_matrix_and_only_measurement_metric(self):
        stdout = io.StringIO()
        # Returning the correct Z output should give zero discrepancy for tket.
        with (
            patch.object(
                bugs, "optimize", return_value=(bugs.mp.matrix([[1, 0], [0, -1]]), 1)
            ),
            patch("sys.stdout", stdout),
        ):
            bugs.main(["--tools", "tket"])
        report = stdout.getvalue()
        self.assertIn("output: 1 gates", report)
        self.assertIn("Change:   0.000000 percentage points", report)
        self.assertNotIn("UNITARY", report)
        self.assertNotIn("Import error", report)


if __name__ == "__main__":
    unittest.main()
