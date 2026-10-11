#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "qiskit==2.5.1", "pennylane==0.45.1", "cirq-core==1.7.0",
#   "bqskit==1.2.1", "bqskitrs==0.4.1", "quizx==0.3.0",
#   "pyzx==0.9.0", "pytket==2.18.5", "mpmath==1.3.0", "ply==3.11",
# ]
# ///
"""Six QASM reproducers. Compare P(+1) after preparing |+X> and measuring X or Y."""

import argparse
import math
from decimal import Decimal
from fractions import Fraction

import mpmath as mp

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
TOOLS = ("qiskit", "pennylane", "cirq", "bqskit", "quizx", "tket")
mp.mp.dps = 100


def literal(value):
    return format(Decimal.from_float(value), "f")


def examples():
    d = literal(3e-6)
    echo = f"gate echo q {{ p({d}) q; p({d}) q; x q; p(-{d}) q; p(-{d}) q; x q; }}\n"
    rotations = "rz(10000000000000000*pi) q[0];\nrz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n"
    return {
        # tool: (QASM, exact ideal angle, measurement, readable input, ideal output)
        "qiskit": (
            HEADER
            + "gate p(a) q { u1(a) q; }\n"
            + echo
            + "qreg q[1];\n"
            + "echo q[0];\n" * 100_000,
            Fraction(3e-6) * 400_000,
            "Y",
            "[P(d); P(d); X; P(-d); P(-d); X] x 100,000; d=3e-6",
            "RZ(~1.2)",
        ),
        "pennylane": (
            HEADER
            + "qreg q[1];\n"
            + "".join(f"rz({a}) q[0];\n" for a in (10**16, 1, 1, 1, -(10**16))),
            Fraction(3),
            "X",
            "RZ(1e16); RZ(1); RZ(1); RZ(1); RZ(-1e16)",
            "RZ(3)",
        ),
        "cirq": (
            HEADER + "qreg q[1];\n" + rotations,
            "pi",
            "X",
            "RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)",
            "Z",
        ),
        "bqskit": (
            HEADER
            + "gate audit_block(a) q { rz(a) q; }\nqreg q[1];\n"
            + f"audit_block({literal(1e-4)}) q[0];\n" * 10_000,
            Fraction(1e-4) * 10_000,
            "Y",
            "RZ(1e-4) x 10,000; one gate per block",
            "RZ(~1)",
        ),
        "quizx": (
            HEADER + "qreg q[1];\nrz(10000000000000114) q[0];\n",
            Fraction(10000000000000114),
            "X",
            "RZ(10000000000000114)",
            "RZ(~-3.13309558725)",
        ),
        "tket": (
            HEADER + "qreg q[1];\n" + rotations,
            "pi",
            "X",
            "RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)",
            "Z",
        ),
    }


def ideal_probability(angle, axis):
    theta = mp.pi if angle == "pi" else mp.mpf(angle.numerator) / angle.denominator
    return (1 + (mp.cos(theta) if axis == "X" else mp.sin(theta))) / 2


def output_probability(unitary, axis):
    # <+X|U|+X> or <+Y|U|+X>; square the amplitude for outcome +1.
    a = (unitary[0, 0] + unitary[0, 1]) / 2
    b = (unitary[1, 0] + unitary[1, 1]) / 2
    probability = float(abs(a + (b if axis == "X" else -1j * b)) ** 2)
    if not math.isfinite(probability):
        raise ValueError("Nonfinite output probability")
    return probability


def qiskit_load(source):
    from qiskit import qasm2
    from qiskit.circuit.library import PhaseGate

    instructions = [i for i in qasm2.LEGACY_CUSTOM_INSTRUCTIONS if i.name != "p"]
    instructions.append(qasm2.CustomInstruction("p", 1, 1, PhaseGate))
    return qasm2.loads(source, custom_instructions=instructions).decompose(
        gates_to_decompose=["echo"]
    )


def optimize(tool, source):
    """Return the output's native matrix and gate count. All inputs are parsed QASM."""
    if tool == "qiskit":
        from qiskit.quantum_info import Operator
        from qiskit.transpiler import PassManager
        from qiskit.transpiler.passes import (
            CommutativeCancellation,
            InverseCancellation,
        )

        output = PassManager(
            [CommutativeCancellation(basis_gates=["p", "x"]), InverseCancellation()]
        ).run(qiskit_load(source))
        return Operator(output).data, len(output.data)
    if tool == "pennylane":
        import pennylane as qml

        # Explicit QASM-parser adapter: this fixture contains only RZ gates.
        parsed = qiskit_load(source)
        assert all(i.operation.name == "rz" for i in parsed.data)
        tape = qml.tape.QuantumScript(
            [qml.RZ(i.operation.params[0], 0) for i in parsed.data]
        )
        output = qml.transforms.merge_rotations(tape, atol=0)[0][0]
        return qml.matrix(output, wire_order=[0]), len(output.operations)
    if tool == "cirq":
        import cirq
        from cirq.contrib.qasm_import import circuit_from_qasm

        circuit = circuit_from_qasm(source)
        output = cirq.eject_z(circuit, atol=0)
        return output.unitary(qubit_order=sorted(circuit.all_qubits())), len(
            list(output.all_operations())
        )
    if tool == "bqskit":
        from bqskit.compiler import Compiler
        from bqskit.ir.lang.qasm2 import OPENQASM2Language
        from bqskit.passes import ForEachBlockPass, ScanningGateRemovalPass, UnfoldPass

        circuit = OPENQASM2Language().decode(source)
        with Compiler(num_workers=1) as compiler:
            output = compiler.compile(
                circuit,
                [
                    ForEachBlockPass(ScanningGateRemovalPass()),
                    UnfoldPass(),
                ],
            )
        return output.get_unitary(), output.num_operations
    if tool == "quizx":
        import quizx

        graph = quizx.qasm(source)
        quizx.full_simp(graph)
        output = quizx.extract_circuit(graph)
        assert output.qubits == 1
        return output.to_matrix(), len(output.gates)
    from pytket.passes import RemoveRedundancies
    from pytket.qasm import circuit_from_qasm_str

    circuit = circuit_from_qasm_str(source)
    RemoveRedundancies().apply(circuit)
    return circuit.get_unitary(), circuit.n_gates


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tools", nargs="+", choices=TOOLS, default=TOOLS)
    args = parser.parse_args(argv)
    cases = examples()
    print("Prepare |+X>; measure X or Y. Compare probability of outcome +1.")
    for tool in dict.fromkeys(args.tools):
        source, angle, axis, description, expected = cases[tool]
        print(f"\n{tool}\n  Input:    {description}", flush=True)
        unitary, count = optimize(tool, source)
        ideal = ideal_probability(angle, axis)
        actual = output_probability(unitary, axis)
        print(
            f"  Expected: {expected}; output: {'identity' if count == 0 else f'{count} gates'}"
        )
        print(
            f"  Measure {axis}: expected {float(100 * ideal):.6f}% -> output {100 * actual:.6f}%"
        )
        print(f"  Change:   {float(100 * abs(ideal - actual)):.6f} percentage points")


if __name__ == "__main__":
    main()
