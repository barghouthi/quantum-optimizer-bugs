#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "qiskit==2.5.1", "pennylane==0.45.1", "cirq-core==1.7.0",
#   "bqskit==1.2.1", "bqskitrs==0.4.1", "quizx==0.3.0",
#   "pyzx==0.9.0", "pytket==2.18.5", "mpmath==1.3.0", "ply==3.11",
#   "sympy==1.14.0", "mqt.qmap==3.10.0", "mqt.core==3.10.1",
# ]
# ///
"""QASM reproducers. Prepare every qubit in |+>; measure X on q[0]."""

import math
import re
import sys
from fractions import Fraction

import mpmath as mp

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
mp.mp.dps = 100
# MQT QMAP 3.10.0 requires MQT Core 3.10, so MQT Core 4 runs in its own environment.
# Pin every package in that environment, including indirect dependencies.
MQT_CORE4 = [
    "mqt.core==4.0.0", "qiskit==2.5.1", "mpmath==1.3.0", "numpy==2.5.4",
    "scipy==1.18.1", "rustworkx==0.18.1", "nanobind-backend==1.0.0",
    "dill==0.4.1", "stevedore==5.9.1", "typing-extensions==4.16.0",
]


def examples():
    rotations = "rz(10000000000000000*pi) q[0];\nrz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n"
    return {
        # tool: (QASM, readable input, ideal output)
        "qiskit": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "rz(1e17); rz(pi); rz(-1e17)",
            "Z",
        ),
        "pennylane": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "rz(1e17); rz(pi); rz(-1e17)",
            "Z",
        ),
        "cirq": (
            HEADER + "qreg q[1];\n" + rotations,
            "rz(1e16*pi); rz(pi); rz(-1e16*pi)",
            "Z",
        ),
        "quizx": (
            HEADER
            + "qreg q[1];\nrz(2147483647*pi/4294967295) q[0];\nrz(2147483648*pi/4294967297) q[0];\n",
            "rz(2147483647*pi/4294967295); rz(2147483648*pi/4294967297)",
            "rz(~pi)",
        ),
        "tket": (
            HEADER + "qreg q[1];\n" + rotations,
            "rz(1e16*pi); rz(pi); rz(-1e16*pi)",
            "Z",
        ),
        "bqskit": (
            HEADER
            + "gate g(a) x,y { rz(a) x; id y; }\nqreg q[2];\n"
            + "g(0) q[1],q[0];\ng(pi) q[0],q[1];\n",
            "g(0) q[1],q[0]; g(pi) q[0],q[1]; g(a) = rz(a) on its first qubit",
            "Z on q[0]",
        ),
        "mqt": (
            HEADER + "qreg q[2];\ncrz(2*pi) q[0],q[1];\n",
            "crz(2*pi) q[0],q[1]",
            "Z on q[0]",
        ),
        "mqt core 4": (
            'OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[1] q;\n'
            + "ctrl @ pow(80599999999716496) @ gphase(1) q[0];\n",
            "ctrl @ pow(80599999999716496) @ gphase(1) q[0] (= p(80599999999716496))",
            "p(80599999999716496), almost identity",
        ),
    }


def ideal_probability(source):
    """Exact diagonal phases from fixture QASM, then 100-digit probability.

    This small reference supports rz, crz, id, ctrl @ pow(n) @ gphase(a),
    and the fixture's g(a) macro.
    Unsupported syntax raises an error instead of supplying a guessed reference.
    """
    import sympy as sp  # not installed in the MQT Core 4 environment

    definition = re.search(r"gate g\(a\) x,y \{(.*?)\}", source, re.DOTALL)
    if definition:
        body = definition[1]
        source = source.replace(definition[0], "")

        def expand(match):
            bindings = dict(zip(("a", "x", "y"), match.groups()))
            bindings["a"] = f"({bindings['a']})"
            return re.sub(r"\b(a|x|y)\b", lambda m: bindings[m[0]], body)

        source = re.sub(r"g\(([^;]*)\) (q\[\d+\]),(q\[\d+\]);", expand, source)
    qubits = int(re.search(r"(?:qreg q\[|qubit\[)(\d+)\]", source)[1])
    phases = [sp.S.Zero] * (1 << qubits)
    for statement in source.split(";"):
        statement = statement.strip()
        if not statement or re.fullmatch(
            r'OPENQASM [23]\.0|include "(qelib1|stdgates)\.inc"|qreg q\[\d+\]|qubit\[\d+\] q|id q\[\d+\]',
            statement,
        ):
            continue
        # ctrl @ pow(n) @ gphase(a) applies phase n*a when the control is 1.
        controlled = re.fullmatch(r"ctrl @ pow\((\d+)\) @ gphase\((.*)\) q\[(\d+)\]", statement)
        if controlled:
            angle = int(controlled[1]) * sp.sympify(controlled[2], locals={"pi": sp.pi}, rational=True)
            for basis in range(len(phases)):
                if (basis >> (qubits - 1 - int(controlled[3]))) & 1:
                    phases[basis] += angle
            continue
        gate = re.fullmatch(r"(rz|crz)\((.*)\) (q\[\d+\](?:,q\[\d+\])?)", statement)
        if not gate:
            raise ValueError(f"Unsupported reference gate: {statement}")
        angle = sp.sympify(gate[2], locals={"pi": sp.pi}, rational=True)
        wires = [int(q) for q in re.findall(r"q\[(\d+)\]", gate[3])]
        if len(wires) != (1 if gate[1] == "rz" else 2):
            raise ValueError(f"Wrong gate arity: {statement}")
        for basis in range(len(phases)):
            bits = [(basis >> (qubits - 1 - q)) & 1 for q in wires]
            if gate[1] == "rz" or bits[0]:
                phases[basis] += angle * (bits[-1] - sp.Rational(1, 2))
    half = len(phases) // 2
    probability = (
        1 + sum(sp.cos(phases[i] - phases[i + half]) for i in range(half)) / half
    ) / 2
    return mp.mpf(str(probability.evalf(100)))


def output_probability(unitary):
    # Prepare every qubit in |+>; measure X on the most significant qubit.
    dimension = len(unitary)
    state = [
        sum(unitary[i, j] for j in range(dimension)) / math.sqrt(dimension)
        for i in range(dimension)
    ]
    half = dimension // 2
    probability = float(
        sum(abs(state[i] + state[i + half]) ** 2 for i in range(half)) / 2
    )
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

        # Explicit QASM-parser adapter: this fixture contains only rz gates.
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
        from bqskit.passes import BlockConversionPass
        from bqskit.passes.processing.extract_diagonal import ExtractDiagonalPass

        circuit = OPENQASM2Language().decode(source)
        with Compiler(num_workers=1) as compiler:
            output = compiler.compile(
                circuit,
                [
                    BlockConversionPass("variable"),
                    ExtractDiagonalPass(
                        success_threshold=1e-12,
                        instantiate_options={"seed": 75, "multistarts": 1},
                    ),
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
    if tool == "tket":
        from pytket.passes import RemoveRedundancies
        from pytket.qasm import circuit_from_qasm_str

        circuit = circuit_from_qasm_str(source)
        RemoveRedundancies().apply(circuit)
        return circuit.get_unitary(), circuit.n_gates
    if tool == "mqt":
        from mqt.core import load
        from mqt.core.dd import build_unitary
        from mqt.qmap import sc

        circuit = load(source)
        config = sc.Configuration()
        config.add_measurements_to_mapped_circuit = False
        config.add_barriers_between_layers = False
        output, _ = sc.map_(circuit, sc.Architecture(2, {(0, 1), (1, 0)}), config)
        assert output.initial_layout == output.output_permutation == {0: 0, 1: 1}
        # MQT matrices place q[0] in the least significant bit. Reorder to MSB.
        matrix = build_unitary(output)
        order = [0, 2, 1, 3]
        return matrix[order][:, order], len(output)
    if tool == "mqt core 4":
        import json
        import os
        import subprocess

        import numpy as np

        command = [os.environ.get("UV", "uv"), "run", "--quiet", "--isolated", "--no-project"]
        command += ["--python", "3.12", *(f"--with={p}" for p in MQT_CORE4)]
        command += ["python", __file__, "--mqt-core4"]
        result = subprocess.run(command, input=source, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stderr)
        matrix, count = json.loads(result.stdout)
        return np.array([[complex(*z) for z in row] for row in matrix]), count
    raise ValueError(f"Unknown tool: {tool}")


def mqt_core4(source):
    """Run in the MQT Core 4 environment; print the output matrix and gate count."""
    import json

    from mqt.core import mlir
    from qiskit.quantum_info import Operator

    imported = mlir.compile_program(source, output=mlir.OutputFormat.QC_IMPORT).ir
    assert "0x4371E5945E0E12C9" in imported  # import preserves K's binary64 bits
    output = mlir.compile_program(
        source, output=mlir.OutputFormat.QCO_OPTIMIZED
    ).to_qiskit()
    matrix = [[[z.real, z.imag] for z in row] for row in Operator(output).data]
    print(json.dumps([matrix, len(output.data)]))


def main():
    print("Prepare every qubit in |+>; measure X on q[0]. Compare P(+1).")
    reproduced_all = True
    for tool, (source, description, expected) in examples().items():
        print(f"\n{tool}\n  Input:    {description}", flush=True)
        unitary, count = optimize(tool, source)
        ideal = ideal_probability(source)
        actual = output_probability(unitary)
        print(f"  Expected: {expected}; output: {count} gates")
        print(
            f"  Measure X: expected {float(100 * ideal):.6f}% -> output {100 * actual:.6f}%"
        )
        error = float(100 * abs(ideal - actual))
        reproduced = error > 99
        reproduced_all &= reproduced
        print(f"  Error:    {error:.6f} percentage points")
        print(f"  Reproduced (>99 pp): {'yes' if reproduced else 'no'}")
    if not reproduced_all:
        raise SystemExit("At least one reported failure did not reproduce.")


if __name__ == "__main__":
    if sys.argv[1:] == ["--mqt-core4"]:
        mqt_core4(sys.stdin.read())
    else:
        main()
