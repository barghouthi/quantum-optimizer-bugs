#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "qiskit==2.5.1", "pennylane==0.45.1", "cirq-core==1.7.0",
#   "bqskit==1.2.1", "bqskitrs==0.4.1", "quizx==0.3.0",
#   "pyzx==0.9.0", "pytket==2.18.5", "mpmath==1.3.0", "ply==3.11",
# ]
# ///
"""Six finite QASM reproducers. Run: uv run bugs.py [--tools cirq tket ...].

All circuits are diagonal up to global phase after their X gates cancel.
The reference sums rational angles exactly, keeping a separate coefficient of
mathematical pi. Only final trigonometric evaluations use mpmath (100/160 digits).
This specialization avoids general circuit simulation and is not an interval proof.
"""

import argparse
import importlib.metadata
import json
import time
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import mpmath

TOOLS = ("qiskit", "pennylane", "cirq", "bqskit", "quizx", "tket")
HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
REFERENCE = "https://arxiv.org/abs/2302.14066 (Definition 1.4, Eq. 3)"


def literal(value):
    """An exact decimal spelling of a fixture's stored binary64 value."""
    return format(Decimal.from_float(value), "f")


@dataclass(frozen=True)
class Phase:
    radians: Fraction = Fraction(0)
    pi: Fraction = Fraction(0)

    def __sub__(self, other):
        return Phase(self.radians - other.radians, self.pi - other.pi)

    def value(self, ctx):
        def rational(q):
            return ctx.mpf(q.numerator) / q.denominator

        return rational(self.radians) + rational(self.pi) * ctx.pi


def phase(gates):
    """Relative phase of P/RZ/pi-power/X gates, ignoring one global phase.

    Maintain U = global * X^b P(theta). P(a) adds (-1)^b * a; X toggles b.
    Count signed raw angles before rationalizing, so even
    the 600,000-gate input needs only a few Fraction conversions.
    """
    counts, sign = Counter(), 1
    for name, angle in gates:
        name = name.lower()
        if name == "x":
            sign = -sign
        elif name in ("p", "u1", "rz", "rzpi", "zpow"):
            counts[(name in ("rzpi", "zpow"), angle)] += sign
        else:
            raise ValueError(f"Unsupported reference gate: {name}")
    if sign != 1:
        raise ValueError("This reference requires an even number of X gates")
    totals = [Fraction(0), Fraction(0)]
    for (pi_units, angle), count in counts.items():
        totals[pi_units] += Fraction(angle) * count
    return Phase(*totals)


def distance(ctx, before, after):
    """min_phi ||U - exp(i phi) V||_op for these diagonal qubit unitaries."""
    delta = (before - after).value(ctx)
    delta = ctx.atan2(ctx.sin(delta), ctx.cos(delta))  # Principal [-pi, pi].
    return 2 * ctx.sin(abs(delta) / 4)


def metrics(ideal, imported, output):
    evaluations = []
    for digits in (100, 160):
        ctx = mpmath.mp.clone()
        ctx.dps = digits
        theta, eta = ideal.value(ctx), output.value(ctx)
        probabilities = {
            "X": ((1 + ctx.cos(theta)) / 2, (1 + ctx.cos(eta)) / 2),
            "Y": ((1 + ctx.sin(theta)) / 2, (1 + ctx.sin(eta)) / 2),
        }
        values = {
            "unitary_error": distance(ctx, ideal, output),
            "import_error": distance(ctx, ideal, imported),
            "optimization_error": distance(ctx, imported, output),
        }
        for axis, (before, after) in probabilities.items():
            values[f"expected_{axis}"] = before
            values[f"output_{axis}"] = after
        evaluations.append((ctx, values))
    ctx, values = evaluations[1]
    previous_ctx, previous = evaluations[0]
    if any(
        abs(value - ctx.mpf(previous_ctx.nstr(previous[key], 100))) >= ctx.mpf("1e-70")
        for key, value in values.items()
    ):
        raise ArithmeticError("Reference failed the 100/160-digit stability check")
    axis = max("XY", key=lambda a: abs(values[f"expected_{a}"] - values[f"output_{a}"]))
    values["probability_change"] = abs(
        values[f"expected_{axis}"] - values[f"output_{axis}"]
    )
    return {key: ctx.nstr(value, 80) for key, value in values.items()} | {
        "measure": axis
    }


def fixture(tool):
    """QASM and independent analytic ideal, specialized to six fixed examples."""
    if tool == "qiskit":
        d = literal(3e-6)
        source = (
            HEADER
            + "gate p(a) q { u1(a) q; }\n"
            + f"gate echo q {{ p({d}) q; p({d}) q; x q; p(-{d}) q; p(-{d}) q; x q; }}\n"
            + "qreg q[1];\n"
            + "echo q[0];\n" * 100_000
        )
        return (
            source,
            Phase(4 * 100_000 * Fraction(3e-6)),
            "[P(d); P(d); X; P(-d); P(-d); X] x 100,000, d=3e-6",
            "RZ(1.2 approximately)",
        )
    if tool == "pennylane":
        source = (
            HEADER
            + "qreg q[1];\n"
            + "".join(f"rz({a}) q[0];\n" for a in (10**16, 1, 1, 1, -(10**16)))
        )
        return (
            source,
            Phase(Fraction(3)),
            "RZ(1e16); RZ(1); RZ(1); RZ(1); RZ(-1e16)",
            "RZ(3)",
        )
    if tool in ("cirq", "tket"):
        source = (
            HEADER
            + "qreg q[1];\nrz(10000000000000000*pi) q[0];\nrz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n"
        )
        return (
            source,
            Phase(pi=Fraction(1)),
            "RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)",
            "Z, up to global phase",
        )
    if tool == "bqskit":
        source = (
            HEADER
            + "gate audit_block(a) q { rz(a) q; }\nqreg q[1];\n"
            + f"audit_block({literal(1e-4)}) q[0];\n" * 10_000
        )
        return (
            source,
            Phase(10_000 * Fraction(1e-4)),
            "[RZ(1e-4)] x 10,000, one gate per block",
            "RZ(1.0 approximately)",
        )
    source = HEADER + "qreg q[1];\nrz(10000000000000114) q[0];\n"
    return (
        source,
        Phase(Fraction(10000000000000114)),
        "RZ(10000000000000114 radians)",
        "RZ(-3.1330955872497763...), up to global phase",
    )


def qiskit_load(source):
    from qiskit import qasm2
    from qiskit.circuit.library import PhaseGate

    instructions = [i for i in qasm2.LEGACY_CUSTOM_INSTRUCTIONS if i.name != "p"]
    instructions.append(qasm2.CustomInstruction("p", 1, 1, PhaseGate))
    return qasm2.loads(source, custom_instructions=instructions).decompose(
        gates_to_decompose=["echo"]
    )


def qiskit_gates(circuit):
    return [
        (i.operation.name, i.operation.params[0] if i.operation.params else 0)
        for i in circuit.data
    ]


def run_native(tool, source):
    """Every native input originates in a QASM loader, never a gate builder."""
    if tool == "qiskit":
        from qiskit.transpiler import PassManager
        from qiskit.transpiler.passes import (
            CommutativeCancellation,
            InverseCancellation,
        )

        circuit = qiskit_load(source)
        imported = qiskit_gates(circuit)
        output = PassManager(
            [CommutativeCancellation(basis_gates=["p", "x"]), InverseCancellation()]
        ).run(circuit)
        return (
            imported,
            qiskit_gates(output),
            "CommutativeCancellation + InverseCancellation",
            "qiskit.qasm2.loads",
            "Small local cutoffs accumulate into a large rotation loss.",
        )
    if tool == "pennylane":
        import pennylane as qml

        parsed = qiskit_load(source)
        imported = qiskit_gates(parsed)
        if any(name != "rz" for name, _ in imported):
            raise ValueError("PennyLane adapter accepts only RZ for this fixture")
        tape = qml.tape.QuantumScript([qml.RZ(angle, 0) for _, angle in imported])
        output = qml.transforms.merge_rotations(tape, atol=0)[0][0]
        gates = [(op.name.lower(), op.parameters[0]) for op in output.operations]
        return (
            imported,
            gates,
            "merge_rotations(atol=0)",
            "Qiskit QASM parser -> explicit PennyLane RZ adapter",
            "Floating-point addition loses the three-radian residual, even at atol=0.",
        )
    if tool == "cirq":
        import cirq
        from cirq.contrib.qasm_import import circuit_from_qasm

        def describe(circuit):
            result = []
            for op in circuit.all_operations():
                if not isinstance(op.gate, cirq.ZPowGate):
                    raise TypeError(f"Unsupported Cirq gate: {op.gate}")
                result.append(("zpow", op.gate.exponent))
            return result

        circuit = circuit_from_qasm(source)
        imported = describe(circuit)
        return (
            imported,
            describe(cirq.eject_z(circuit, atol=0)),
            "eject_z(atol=0)",
            "cirq.contrib.qasm_import.circuit_from_qasm",
            "Large-exponent addition erases the middle Z.",
        )
    if tool == "tket":
        from pytket.passes import RemoveRedundancies
        from pytket.qasm import circuit_from_qasm_str

        def describe(circuit):
            result = []
            for command in circuit.get_commands():
                if command.op.type.name != "Rz":
                    raise ValueError(f"Unsupported tket gate: {command.op.type.name}")
                result.append(("rzpi", command.op.params[0]))
            return result

        circuit = circuit_from_qasm_str(source)
        imported = describe(circuit)
        RemoveRedundancies().apply(circuit)
        return (
            imported,
            describe(circuit),
            "RemoveRedundancies",
            "pytket.qasm.circuit_from_qasm_str",
            "Import preserves Z; the optimization pass erases it.",
        )
    if tool == "quizx":
        import quizx

        def describe(graph):
            circuit = quizx.extract_circuit(graph)
            if circuit.qubits != 1:
                raise ValueError("QuiZX extraction lost the one-qubit interface")
            gates = []
            for gate in circuit.gates:
                if gate.name != "ZPhase":
                    raise ValueError(f"Unsupported QuiZX gate: {gate.name}")
                angle = Fraction(gate.phase) * (
                    -1 if getattr(gate, "adjoint", False) else 1
                )
                gates.append(("rzpi", angle))
            return gates

        imported = describe(quizx.qasm(source))
        graph = quizx.qasm(
            source
        )  # Fresh parse: graph.copy can lose interface metadata.
        quizx.full_simp(graph)
        return (
            imported,
            describe(graph),
            "qasm -> full_simp -> extract_circuit",
            "quizx.qasm",
            "The error occurs during import; full_simp adds no error.",
        )
    from bqskit.compiler import Compiler
    from bqskit.ir.lang.qasm2 import OPENQASM2Language
    from bqskit.passes import ForEachBlockPass, ScanningGateRemovalPass, UnfoldPass

    def describe(circuit):
        flat = circuit.copy()
        flat.unfold_all()
        return [
            (op.gate.name.removesuffix("Gate").lower(), op.params[0]) for op in flat
        ]

    circuit = OPENQASM2Language().decode(source)
    imported = describe(circuit)
    with Compiler(num_workers=1) as compiler:
        output, data = compiler.compile(
            circuit,
            [
                ForEachBlockPass(ScanningGateRemovalPass(), calculate_error_bound=True),
                UnfoldPass(),
            ],
            request_data=True,
        )
    note = f"Configured local approximation workflow; recorded bound {data.error:.6g} is not enforced as a global limit."
    return (
        imported,
        describe(output),
        "ForEachBlockPass(scan) + UnfoldPass",
        "bqskit.OPENQASM2Language.decode",
        note,
    )


def short(value):
    ctx = mpmath.mp.clone()
    ctx.dps = 100
    return ctx.nstr(ctx.mpf(value), 12)


def print_result(row):
    print(
        f"\n{'=' * 78}\n{row['tool'].upper()} {row['version']} | {row['seconds']:.2f} seconds\n{'=' * 78}"
    )

    def line(label, value):
        print(f"{label:<23} {value}")

    for label, key in (
        ("Circuit", "circuit"),
        ("Expected", "expected"),
        ("Native pass", "optimizer"),
        ("QASM loader", "loader"),
        ("Output", "output"),
    ):
        line(label, row[key])
    line(
        "QASM / parsed / output",
        f"{row['input_gates']:,} / {row['imported_gates']:,} / {row['output_gates']:,} gates",
    )
    for label, key in (
        ("UNITARY ERROR E", "unitary_error"),
        ("Import error E", "import_error"),
        ("Optimization error E", "optimization_error"),
    ):
        line(label, short(row[key]))
    axis = row["measure"]
    line("Measurement", f"Prepare |+X>; measure {axis}; probability of outcome +1")
    line(
        "Expected -> output",
        f"{100 * float(row[f'expected_{axis}']):.9f}% -> {100 * float(row[f'output_{axis}']):.9f}%",
    )
    line(
        "PROBABILITY CHANGE",
        f"{100 * float(row['probability_change']):.9f} percentage points",
    )
    line("Reference check", "100 / 160 digits agree within absolute 1e-70")
    line("Interpretation", row["note"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tools", nargs="+", choices=TOOLS, default=list(TOOLS))
    parser.add_argument(
        "--qasm-dir", type=Path, help="Export the exact standalone QASM inputs"
    )
    parser.add_argument(
        "--json", type=Path, help="Save full precision results and QASM inputs"
    )
    args = parser.parse_args(argv)
    print("Six finite quantum optimizer reproducers (all inputs are OpenQASM 2.0)")
    print("E = min over real phi of ||U_expected - exp(i*phi) U_output||_op.")
    print(
        "phi is ONE global phase for the entire matrix. E is unitless, not a percentage."
    )
    print("Reference:", REFERENCE)
    print("Measurement changes are separate physical witnesses, not worst-case bounds.")
    rows, failures = [], 0
    for tool in dict.fromkeys(args.tools):
        print(f"\nRunning {tool}...", flush=True)
        started = time.perf_counter()
        source, ideal, circuit, expected = fixture(tool)
        if args.qasm_dir:
            args.qasm_dir.mkdir(parents=True, exist_ok=True)
            (args.qasm_dir / f"{tool}.qasm").write_text(source)
        try:
            imported, output, optimizer, loader, note = run_native(tool, source)
            package = {"cirq": "cirq-core", "tket": "pytket"}.get(tool, tool)
            row = {
                "tool": tool,
                "version": importlib.metadata.version(package),
                "circuit": circuit,
                "expected": expected,
                "input_qasm": source,
                "optimizer": optimizer,
                "loader": loader,
                "note": note,
                "input_gates": {
                    "qiskit": 600_000,
                    "pennylane": 5,
                    "cirq": 3,
                    "bqskit": 10_000,
                    "quizx": 1,
                    "tket": 3,
                }[tool],
                "imported_gates": len(imported),
                "output_gates": len(output),
                "output": "identity" if not output else repr(output),
                **metrics(ideal, phase(imported), phase(output)),
                "seconds": time.perf_counter() - started,
            }
            print_result(row)
        except Exception as exc:  # noqa: BLE001 -- Report native failures and continue the other tools.
            failures += 1
            row = {
                "tool": tool,
                "input_qasm": source,
                "error": f"{type(exc).__name__}: {exc}",
            }
            print(f"FAILED {tool}: {row['error']}")
        rows.append(row)
    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
    print(f"\nCompleted {len(rows)} examples; {failures} execution/reference failures.")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
