#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "qiskit==2.5.1", "pennylane==0.45.1", "cirq-core==1.7.0",
#   "bqskit==1.2.1", "bqskitrs==0.4.1", "quizx==0.3.0",
#   "pyzx==0.9.0", "pytket==2.18.5", "mpmath==1.3.0", "ply==3.11",
# ]
# ///
"""Six QASM reproducers. Compare P(+1) after preparing |+> and measuring X."""

import math
from fractions import Fraction

import mpmath as mp

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
mp.mp.dps = 100


def examples():
    rotations = "rz(10000000000000000*pi) q[0];\nrz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n"
    return {
        # tool: (QASM, exact ideal angle, readable input, ideal output)
        "qiskit": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "pi",
            "RZ(1e17); RZ(pi); RZ(-1e17)",
            "Z",
        ),
        "pennylane": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "pi",
            "RZ(1e17); RZ(pi); RZ(-1e17)",
            "Z",
        ),
        "cirq": (
            HEADER + "qreg q[1];\n" + rotations,
            "pi",
            "RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)",
            "Z",
        ),
        "bqskit": (
            HEADER
            + "gate audit_block(a) q { rz(a) q; }\nqreg q[1];\n"
            + "audit_block(pi/12000) q[0];\n" * 12_000,
            "pi",
            "RZ(pi/12000) x 12,000; one gate per block",
            "Z",
        ),
        "quizx": (
            HEADER
            + "qreg q[1];\nrz(2147483647*pi/4294967295) q[0];\nrz(2147483648*pi/4294967297) q[0];\n",
            lambda: mp.pi * (1 - mp.mpf(2**32) / (2**64 - 1)),
            "RZ(2147483647*pi/4294967295); RZ(2147483648*pi/4294967297)",
            "RZ(~pi)",
        ),
        "tket": (
            HEADER + "qreg q[1];\n" + rotations,
            "pi",
            "RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)",
            "Z",
        ),
    }


def ideal_probability(angle):
    theta = angle() if callable(angle) else mp.pi
    return (1 + mp.cos(theta)) / 2


def output_probability(unitary):
    # <+|U|+>; square the amplitude for outcome +1.
    amplitude = sum(unitary[i, j] for i in range(2) for j in range(2)) / 2
    probability = float(abs(amplitude) ** 2)
    if not math.isfinite(probability):
        raise ValueError("Nonfinite output probability")
    return probability


def qiskit_load(source):
    from qiskit import qasm2

    return qasm2.loads(source)


def optimize(tool, source):
    """Return the output's native matrix and gate count. All inputs are parsed QASM."""
    if tool == "qiskit":
        from qiskit.quantum_info import Operator
        from qiskit.transpiler import PassManager
        from qiskit.transpiler.passes import CommutativeCancellation

        output = PassManager([CommutativeCancellation(basis_gates=["rz"])]).run(
            qiskit_load(source)
        )
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
        # Check that import preserves both phases before testing simplification.
        assert [graph.phase(v) for v in graph.vertices() if graph.type(v) == 1] == [
            Fraction(2147483647, 4294967295),
            Fraction(2147483648, 4294967297),
        ]
        quizx.full_simp(graph)
        output = quizx.extract_circuit(graph)
        assert output.qubits == 1
        return output.to_matrix(), len(output.gates)
    from pytket.passes import RemoveRedundancies
    from pytket.qasm import circuit_from_qasm_str

    circuit = circuit_from_qasm_str(source)
    RemoveRedundancies().apply(circuit)
    return circuit.get_unitary(), circuit.n_gates


def main():
    print("Prepare |+>; measure X. Compare probability of outcome +1.")
    for tool, (source, angle, description, expected) in examples().items():
        print(f"\n{tool}\n  Input:    {description}", flush=True)
        unitary, count = optimize(tool, source)
        ideal = ideal_probability(angle)
        actual = output_probability(unitary)
        print(
            f"  Expected: {expected}; output: {'identity' if count == 0 else f'{count} gates'}"
        )
        print(
            f"  Measure X: expected {float(100 * ideal):.6f}% -> output {100 * actual:.6f}%"
        )
        print(f"  Change:   {float(100 * abs(ideal - actual)):.6f} percentage points")


if __name__ == "__main__":
    main()
