#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "qiskit==2.5.1", "pennylane==0.45.1", "cirq-core==1.7.0",
#   "bqskit==1.2.1", "bqskitrs==0.4.1", "quizx==0.3.0",
#   "pyzx==0.9.0", "pytket==2.18.5", "mpmath==1.3.0", "ply==3.11",
#   "sympy==1.13.0", "mqt.qmap==3.10.0", "mqt.core==3.10.1",
#   "qrisp==0.9.9", "cuda-quantum-cu13==0.16.0",
#   "qbraid==0.12.1", "pyqasm==1.0.4",
# ]
# ///
"""Small-angle QASM reproducers. Prepare every qubit in |+>; measure X on q[0]."""

import math
import re
import sys

import mpmath as mp

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
mp.mp.dps = 100


def examples():
    tket_block = (
        "rz(6.134e-10) q[0];\nrx(1.064e-10) q[0];\nry(6.317e-11) q[0];\nz q[0];\n"
        "rx(8.644e-11) q[0];\ny q[0];\ncx q[0],q[1];\nx q[0];\ncx q[0],q[1];\n"
    )
    return {
        # tool: (QASM, readable input, ideal output)
        "qiskit": (
            HEADER
            + "qreg q[3];\nh q[1];\nh q[2];\n"
            + "".join(f"crz(2*pi/44430) q[0],q[{1 + k % 2}];\n" for k in range(44430)),
            "h q[1]; h q[2]; 44430 x crz(2*pi/44430), alternating targets q[1] and q[2]",
            "Z on q[0]",
        ),
        "pennylane": (
            HEADER + "qreg q[2];\nh q[1];\nx q[1];\n" + "cu1(pi/351) q[0],q[1];\n" * 351,
            "h q[1]; x q[1]; 351 x cu1(pi/351) q[0],q[1]",
            "Z on q[0]",
        ),
        "cirq": (
            HEADER
            + "qreg q[2];\nh q[1];\ncx q[0],q[1];\n"
            + "rz(0.00001) q[1];\n" * 300000
            + "cx q[0],q[1];\n",
            "h q[1]; cx q[0],q[1]; 300000 x rz(1e-5) q[1]; cx q[0],q[1]",
            "rz(3) on q[0]",
        ),
        "tket": (
            HEADER + "qreg q[2];\ns q[0];\n" + tket_block * 5 + "s q[0];\n",
            "s q[0]; 5 x (four rotations up to 6.134e-10, z, y, cx, x, cx); s q[0]",
            "Z on q[0]",
        ),
        "qrisp": (
            HEADER + "qreg q[1];\n" + "u3(0.00007,0,0) q[0];\n" * 44880,
            "44880 x u3(7e-5,0,0) q[0], a Y rotation by 7e-5",
            "ry(3.1416), almost Y",
        ),
    }


def matrix(name, angles):
    """Exact 100-digit matrix of a qelib1.inc gate. The first operand is the MSB."""
    i = mp.mpc(0, 1)
    one, zero = mp.mpf(1), mp.mpf(0)

    def rotation(axis, theta):
        c, s = mp.cos(theta / 2), mp.sin(theta / 2)
        return {
            "x": [[c, -i * s], [-i * s, c]],
            "y": [[c, -s], [s, c]],
            "z": [[mp.expj(-theta / 2), zero], [zero, mp.expj(theta / 2)]],
        }[axis]

    def u3(theta, phi, lam):
        c, s = mp.cos(theta / 2), mp.sin(theta / 2)
        return [
            [c, -mp.expj(lam) * s],
            [mp.expj(phi) * s, mp.expj(phi + lam) * c],
        ]

    def controlled(m):
        return [
            [one, zero, zero, zero],
            [zero, one, zero, zero],
            [zero, zero, m[0][0], m[0][1]],
            [zero, zero, m[1][0], m[1][1]],
        ]

    r2 = 1 / mp.sqrt(2)
    single = {
        "id": lambda: [[one, zero], [zero, one]],
        "x": lambda: [[zero, one], [one, zero]],
        "y": lambda: [[zero, -i], [i, zero]],
        "z": lambda: [[one, zero], [zero, -one]],
        "h": lambda: [[r2, r2], [r2, -r2]],
        "s": lambda: [[one, zero], [zero, i]],
        "sdg": lambda: [[one, zero], [zero, -i]],
        "t": lambda: [[one, zero], [zero, mp.expj(mp.pi / 4)]],
        "tdg": lambda: [[one, zero], [zero, mp.expj(-mp.pi / 4)]],
        "rx": lambda a: rotation("x", a),
        "ry": lambda a: rotation("y", a),
        "rz": lambda a: rotation("z", a),
        "u1": lambda a: [[one, zero], [zero, mp.expj(a)]],
        "u2": lambda phi, lam: u3(mp.pi / 2, phi, lam),
        "u3": u3,
    }
    if name in single:
        return 1, single[name](*angles)
    # cu3 is left out: Qiskit's qelib1.inc adds a control phase that the
    # official OpenQASM 2 qelib1.inc does not have.
    if name in ("cx", "cy", "cz", "crz", "cu1"):
        return 2, controlled(single[name[1:]](*angles))
    raise ValueError(f"Unsupported reference gate: {name}")


def ideal_probability(source):
    """Exact 100-digit state-vector simulation of simple OpenQASM 2.

    The reference supports only the gates in `matrix` on a register named q.
    It rejects any other statement instead of guessing.
    """
    import sympy as sp  # not installed in the MQT Core 4 environment

    qubits = int(re.search(r"qreg q\[(\d+)\];", source)[1])
    state = [mp.mpc(1) / mp.sqrt(2**qubits)] * (1 << qubits)  # |+...+>
    cache = {}
    for statement in source.split(";"):
        statement = " ".join(statement.split())
        if not statement or re.fullmatch(
            r'OPENQASM 2\.0|include "qelib1\.inc"|qreg q\[\d+\]', statement
        ):
            continue
        gate = re.fullmatch(r"([a-z0-9]+)(?:\((.*)\))? (q\[\d+\](?: ?, ?q\[\d+\])*)", statement)
        if not gate:
            raise ValueError(f"Unsupported reference statement: {statement}")
        key = (gate[1], gate[2])
        if key not in cache:
            angles = [
                mp.mpf(str(sp.sympify(a, locals={"pi": sp.pi}, rational=True).evalf(110)))
                for a in (gate[2].split(",") if gate[2] else [])
            ]
            cache[key] = matrix(gate[1], angles)
        arity, m = cache[key]
        wires = [int(q) for q in re.findall(r"q\[(\d+)\]", gate[3])]
        if len(wires) != arity or len(set(wires)) != arity:
            raise ValueError(f"Wrong gate arity: {statement}")
        shifts = [qubits - 1 - w for w in wires]
        new = list(state)
        for basis in range(1 << qubits):
            if any((basis >> s) & 1 for s in shifts):
                continue  # visit each group once, from its all-zero member
            group = [
                basis | sum(((k >> (arity - 1 - j)) & 1) << shifts[j] for j in range(arity))
                for k in range(1 << arity)
            ]
            for r, out in enumerate(group):
                new[out] = sum(m[r][c] * state[g] for c, g in enumerate(group))
        state = new
    half = len(state) // 2
    return sum(abs(state[k] + state[k + half]) ** 2 for k in range(half)) / 2


def output_probability(output):
    # Prepare every qubit in |+>; measure X on the most significant qubit.
    # The output is a unitary, or the output state for the |+> input.
    dimension = len(output)
    if getattr(output, "ndim", 2) == 1:
        state = output
    else:
        state = [
            sum(output[i, j] for j in range(dimension)) / math.sqrt(dimension)
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
    """Return the output's native matrix (or state) and gate count. All inputs are parsed QASM."""
    if tool == "qiskit":
        from qiskit import transpile
        from qiskit.quantum_info import Operator

        circuit = qiskit_load(source)
        # Check that import preserves every angle before transpiling.
        angles = [i.operation.params[0] for i in circuit.data if i.operation.name == "crz"]
        assert len(angles) == 44430 and set(angles) == {2 * math.pi / 44430}
        output = transpile(circuit, basis_gates=["cx", "u"])
        # Qiskit matrices place q[0] in the least significant bit. Reorder to MSB.
        return Operator(output).reverse_qargs().data, len(output.data)
    if tool == "pennylane":
        import pennylane as qml
        from qiskit.quantum_info import Operator

        # Explicit QASM-parser adapter. Like PennyLane's own QASM importer
        # (pennylane-qiskit), it converts cu1 to a QubitUnitary.
        parsed = qiskit_load(source)
        ops = []
        for i in parsed.data:
            wires = [parsed.find_bit(q).index for q in i.qubits]
            if i.operation.name == "h":
                ops.append(qml.Hadamard(wires))
            elif i.operation.name == "x":
                ops.append(qml.PauliX(wires))
            else:
                assert i.operation.name == "cu1"
                assert i.operation.params == [math.pi / 351]
                matrix = Operator(i.operation).reverse_qargs().data
                ops.append(qml.QubitUnitary(matrix, wires))
        output = qml.transforms.unitary_to_rot(qml.tape.QuantumScript(ops))[0][0]
        return qml.matrix(output, wire_order=[0, 1]), len(output.operations)
    if tool == "cirq":
        import cirq
        from cirq.contrib.qasm_import import circuit_from_qasm

        circuit = circuit_from_qasm(source)
        # Check that import preserves every angle before reordering.
        rotations = [op.gate for op in circuit.all_operations() if isinstance(op.gate, cirq.Rz)]
        assert len(rotations) == 300000 and {g._rads for g in rotations} == {1e-5}
        output = cirq.transformers.insertion_sort_transformer(circuit)
        return output.unitary(qubit_order=sorted(circuit.all_qubits())), len(
            list(output.all_operations())
        )
    if tool == "tket":
        from pytket.passes import FullPeepholeOptimise
        from pytket.qasm import circuit_from_qasm_str

        circuit = circuit_from_qasm_str(source)
        # tket stores angles in half-turns. Check that import preserves them.
        literals = [float(a) for a in re.findall(r"r[xyz]\(([^)]*)\)", source)]
        params = [float(c.op.params[0]) * math.pi for c in circuit if c.op.params]
        assert len(params) == len(literals) == 20
        assert all(abs(p - a) <= 1e-15 * a for p, a in zip(params, literals))
        FullPeepholeOptimise().apply(circuit)
        return circuit.get_unitary(), circuit.n_gates
    if tool == "qrisp":
        import numpy as np
        from qrisp import QuantumCircuit, QuantumVariable

        circuit = QuantumCircuit.from_qasm_str(source)
        # Check that import preserves all gates and angles before compiling.
        assert len(circuit.data) == 44880
        assert all(i.op.name == "u3" and list(i.op.params) == [7e-5, 0, 0] for i in circuit.data)
        variable = QuantumVariable(1)
        variable.qs.append(circuit.to_gate(), [variable[0]])
        output = variable.qs.compile()
        gates = [i for i in output.data if i.op.name != "qb_alloc"]
        # The compiled unitary is complex64; widen it before computing probabilities.
        return np.asarray(output.get_unitary(), dtype=np.complex128), len(gates)
    raise ValueError(f"Unknown tool: {tool}")


def main():
    print("Prepare every qubit in |+>; measure X on q[0]. Compare P(+1).")
    reproduced_all = True
    selected = sys.argv[1:]
    for tool, (source, description, expected) in examples().items():
        if selected and tool not in selected:
            continue
        print(f"\n{tool}\n  Input:    {description}", flush=True)
        output, count = optimize(tool, source)
        ideal = ideal_probability(source)
        actual = output_probability(output)
        gates = "state only" if count is None else f"{count} gates"
        print(f"  Expected: {expected}; output: {gates}")
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
    main()
