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
        "pyzx": (
            HEADER
            + "qreg q[1];\ns q[0];\nx q[0];\n"
            + "rz(0.01) q[0];\n" * 155
            + "x q[0];\nrz(0.01) q[0];\ns q[0];\nx q[0];\nsdg q[0];\nrz(0.01) q[0];\n",
            "s; x; 155 x rz(0.01); x; rz(0.01); s; x; sdg; rz(0.01)",
            "+1 probability 0.010812%",
        ),
        "quizx": (
            HEADER
            + "qreg q[1];\n"
            + "rz(0.01) q[0];\n" * 314
            + "rz(0.001) q[0];\nrz(0.002) q[0];\nrz(0.003) q[0];\n",
            "314 x rz(0.01); rz(0.001); rz(0.002); rz(0.003)",
            "rz(3.146), almost Z",
        ),
        "tket": (
            HEADER + "qreg q[2];\ns q[0];\n" + tket_block * 5 + "s q[0];\n",
            "s q[0]; 5 x (four rotations up to 6.134e-10, z, y, cx, x, cx); s q[0]",
            "Z on q[0]",
        ),
        "bqskit": (
            HEADER
            + "qreg q[4];\n"
            + "rz(0.00028) q[0];\ncx q[0],q[1];\ncx q[0],q[2];\ncx q[0],q[3];\n" * 11220,
            "11220 x (rz(0.00028) q[0]; cx q[0],q[1]; cx q[0],q[2]; cx q[0],q[3])",
            "rz(3.1416) on q[0], almost Z",
        ),
        "qrisp": (
            HEADER + "qreg q[1];\n" + "u3(0.00007,0,0) q[0];\n" * 44880,
            "44880 x u3(7e-5,0,0) q[0], a Y rotation by 7e-5",
            "ry(3.1416), almost Y",
        ),
        "qbraid": (
            HEADER
            + "qreg q[1];\n"
            + "rz(1/408*pi) q[0];\n" * 454
            + "rz(-pi/408) q[0];\n" * 46,
            "454 x rz(1/408*pi); 46 x rz(-pi/408)",
            "rz(pi), equivalent to Z",
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


def exact(value):
    """Exact 100-digit value of a binary64 float, an integer, or a Fraction."""
    from fractions import Fraction

    value = Fraction(value)
    return mp.mpf(value.numerator) / value.denominator


def simulate(qubits, gates):
    """Exact 100-digit state-vector simulation. Return P(+1) of X on q[0].

    Every qubit starts in |+>. Each gate is (name, angles, wires) for a gate in
    `matrix`, or ("unitary", rows, wires) for a raw matrix. The first wire is the MSB.
    """
    state = [mp.mpc(1) / mp.sqrt(2**qubits)] * (1 << qubits)  # |+...+>
    cache = {}
    for name, angles, wires in gates:
        if name == "unitary":
            arity, m = len(wires), angles
        else:
            key = (name, tuple(angles))
            if key not in cache:
                cache[key] = matrix(name, list(angles))
            arity, m = cache[key]
        if len(wires) != arity or len(set(wires)) != arity or max(wires) >= qubits:
            raise ValueError(f"Wrong gate wires: {name} {wires}")
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


def qubit_count(source):
    return int(re.search(r"qreg q\[(\d+)\];", source)[1])


def ideal_probability(source):
    """Exact 100-digit reference for simple OpenQASM 2.

    The reference supports only the gates in `matrix` on a register named q.
    It rejects any other statement instead of guessing.
    """
    import sympy as sp

    gates = []
    angles = {}
    for statement in source.split(";"):
        statement = " ".join(statement.split())
        if not statement or re.fullmatch(
            r'OPENQASM 2\.0|include "qelib1\.inc"|qreg q\[\d+\]', statement
        ):
            continue
        gate = re.fullmatch(r"([a-z0-9]+)(?:\((.*)\))? (q\[\d+\](?: ?, ?q\[\d+\])*)", statement)
        if not gate:
            raise ValueError(f"Unsupported reference statement: {statement}")
        if gate[2] not in angles:
            angles[gate[2]] = [
                mp.mpf(str(sp.sympify(a, locals={"pi": sp.pi}, rational=True).evalf(110)))
                for a in (gate[2].split(",") if gate[2] else [])
            ]
        wires = [int(q) for q in re.findall(r"q\[(\d+)\]", gate[3])]
        gates.append((gate[1], angles[gate[2]], wires))
    return simulate(qubit_count(source), gates)


def float_inputs(source):
    """Replace each angle with the exact value of the binary64 float a parser computes.

    Parsers evaluate an angle such as 2*pi/44430 in binary64 arithmetic.
    The reference on this source shows that the rounding does not explain an error.
    """
    import ast
    import operator
    from decimal import Decimal

    operators = {
        ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
        ast.Div: operator.truediv, ast.USub: operator.neg, ast.UAdd: operator.pos,
    }

    def evaluate(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name) and node.id == "pi":
            return math.pi
        if isinstance(node, ast.BinOp):
            return operators[type(node.op)](evaluate(node.left), evaluate(node.right))
        if isinstance(node, ast.UnaryOp):
            return operators[type(node.op)](evaluate(node.operand))
        raise ValueError(f"Unsupported angle expression: {ast.unparse(node)}")

    def replace(match):
        values = [evaluate(ast.parse(a, mode="eval").body) for a in match[1].split(",")]
        return "(" + ",".join(format(Decimal(v), "f") for v in values) + ") "

    return re.sub(r"\(([^()]*)\) (?=q\[)", replace, source)


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


# Each converter reads a tool's output circuit as a gate list for `simulate`.
# It never calls the tool's own matrix code. `matrix` rejects unknown gate names.


def qiskit_gates(circuit):
    """Qiskit's u and p are u3 and u1 in qelib1.inc. Global phase does not change P(+1)."""
    names = {"u": "u3", "p": "u1", "cp": "cu1"}
    return [
        (
            names.get(i.operation.name, i.operation.name),
            [exact(float(a)) for a in i.operation.params],
            [circuit.find_bit(q).index for q in i.qubits],
        )
        for i in circuit.data
    ]


def pennylane_gates(tape):
    names = {
        "RX": "rx", "RY": "ry", "RZ": "rz", "PhaseShift": "u1", "Hadamard": "h",
        "PauliX": "x", "PauliY": "y", "PauliZ": "z", "CNOT": "cx", "CZ": "cz",
    }
    gates = []
    for op in tape.operations:
        angles = [exact(float(a)) for a in op.parameters]
        wires = [int(w) for w in op.wires]
        if op.name == "Rot":  # Rot(phi, theta, omega) applies RZ(phi), RY(theta), then RZ(omega)
            gates += [("rz", angles[:1], wires), ("ry", angles[1:2], wires), ("rz", angles[2:], wires)]
        elif op.name != "GlobalPhase":
            gates.append((names.get(op.name, op.name), angles, wires))
    return gates


def cirq_gates(circuit):
    import cirq

    names = {cirq.H: "h", cirq.X: "x", cirq.Y: "y", cirq.Z: "z", cirq.CNOT: "cx", cirq.CZ: "cz"}
    gates = []
    for op in circuit.all_operations():
        gate = op.gate
        wires = [int(q.name.removeprefix("q_")) for q in op.qubits]  # QASM import names q_0, q_1, ...
        if isinstance(gate, (cirq.Rx, cirq.Ry, cirq.Rz)):
            gates.append((type(gate).__name__.lower(), [exact(gate._rads)], wires))
        elif gate in names:
            gates.append((names[gate], [], wires))
        elif isinstance(gate, cirq.ZPowGate) and len(wires) == 1:  # u1 up to global phase
            gates.append(("u1", [exact(gate.exponent) * mp.pi], wires))
        else:
            raise ValueError(f"Unsupported Cirq output gate: {gate!r}")
    return gates


def zx_gates(circuit):
    """PyZX and QuiZX extract PyZX circuits. Phases are fractions of pi."""
    names = {"HAD": "h", "NOT": "x", "Z": "z", "S": "s", "T": "t", "CNOT": "cx", "CZ": "cz"}
    gates = []
    for gate in circuit.gates:
        name = type(gate).__name__
        if name in ("ZPhase", "XPhase"):  # rz and rx up to global phase
            rotation = {"ZPhase": "rz", "XPhase": "rx"}[name]
            gates.append((rotation, [exact(gate.phase) * mp.pi], [gate.target]))
        elif name in ("CNOT", "CZ"):
            gates.append((names[name], [], [gate.control, gate.target]))
        elif name in names:
            suffix = "dg" if getattr(gate, "adjoint", False) else ""
            gates.append((names[name] + suffix, [], [gate.target]))
        else:
            raise ValueError(f"Unsupported ZX output gate: {gate}")
    return gates


def tket_gates(circuit):
    """tket angles are half-turns: multiply by pi."""
    names = {
        "Rx": "rx", "Ry": "ry", "Rz": "rz", "U1": "u1", "U3": "u3", "H": "h", "X": "x",
        "Y": "y", "Z": "z", "S": "s", "Sdg": "sdg", "T": "t", "Tdg": "tdg", "CX": "cx", "CZ": "cz",
    }
    gates = []
    for command in circuit:
        kind = command.op.type.name
        angles = [exact(float(a)) * mp.pi for a in command.op.params]
        wires = [q.index[0] for q in command.qubits]
        if kind == "TK1":  # TK1(a, b, c) = Rz(a) Rx(b) Rz(c), so Rz(c) acts first
            gates += [("rz", angles[2:], wires), ("rx", angles[1:2], wires), ("rz", angles[:1], wires)]
        elif kind != "Phase":
            gates.append((names.get(kind, kind), angles, wires))
    return gates


def bqskit_gates(circuit):
    names = {"RXGate": "rx", "RYGate": "ry", "RZGate": "rz", "U3Gate": "u3", "HGate": "h", "CNOTGate": "cx"}
    gates = []
    for op in circuit:
        wires = list(op.location)
        if op.gate.name.startswith("VariableUnitaryGate"):
            # Real parts, then imaginary parts, row by row. BQSKit's own matrix
            # moves this to the closest unitary. The check uses it unchanged.
            size = 2 ** len(wires)
            values = [exact(a) for a in op.params]
            rows = [
                [mp.mpc(values[r * size + c], values[size * size + r * size + c]) for c in range(size)]
                for r in range(size)
            ]
            gates.append(("unitary", rows, wires))
        else:
            gates.append((names.get(op.gate.name, op.gate.name), [exact(a) for a in op.params], wires))
    return gates


def qrisp_gates(circuit):
    gates = []
    for i in circuit.data:
        if i.op.name in ("qb_alloc", "qb_dealloc", "gphase"):
            continue
        gates.append((
            {"p": "u1"}.get(i.op.name, i.op.name),
            [exact(float(a)) for a in i.op.params],
            [circuit.qubits.index(q) for q in i.qubits],
        ))
    return gates


def ionq_gates(program):
    """IonQ JSON gates with rotations in radians."""
    return [(g["gate"], [exact(g["rotation"])], [g["target"]]) for g in program["circuit"]]


def optimize(tool, source):
    """Return the output's native matrix, gate count, and gate list. All inputs are parsed QASM."""
    if tool == "qiskit":
        from qiskit import transpile
        from qiskit.quantum_info import Operator

        circuit = qiskit_load(source)
        # Check that import preserves every angle before transpiling.
        angles = [i.operation.params[0] for i in circuit.data if i.operation.name == "crz"]
        assert len(angles) == 44430 and set(angles) == {2 * math.pi / 44430}
        output = transpile(circuit, basis_gates=["cx", "u"])
        # Qiskit matrices place q[0] in the least significant bit. Reorder to MSB.
        return Operator(output).reverse_qargs().data, len(output.data), qiskit_gates(output)
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
        return qml.matrix(output, wire_order=[0, 1]), len(output.operations), pennylane_gates(output)
    if tool == "cirq":
        import cirq
        from cirq.contrib.qasm_import import circuit_from_qasm

        circuit = circuit_from_qasm(source)
        # Check that import preserves every angle before reordering.
        rotations = [op.gate for op in circuit.all_operations() if isinstance(op.gate, cirq.Rz)]
        assert len(rotations) == 300000 and {g._rads for g in rotations} == {1e-5}
        output = cirq.transformers.insertion_sort_transformer(circuit)
        return (
            output.unitary(qubit_order=sorted(circuit.all_qubits())),
            len(list(output.all_operations())),
            cirq_gates(output),
        )
    if tool == "pyzx":
        from fractions import Fraction

        import pyzx as zx

        circuit = zx.Circuit.from_qasm(source)
        # PyZX stores phases as fractions of pi. Check that import preserves them.
        phases = [g.phase for g in circuit.gates if type(g).__name__ == "ZPhase"]
        assert len(phases) == 157
        assert set(phases) == {Fraction(0.01 / math.pi).limit_denominator(10**8)}
        output = zx.Circuit.from_graph(zx.teleport_reduce(circuit.to_graph()))
        return output.to_matrix(), len(output.gates), zx_gates(output)
    if tool == "quizx":
        from fractions import Fraction

        import quizx

        graph = quizx.qasm(source)
        # QuiZX stores phases as 64-bit rationals of pi. Import approximates each
        # decimal angle with a relative error below 1e-7.
        phases = [graph.phase(v) for v in graph.vertices() if graph.type(v) == 1]
        assert len(phases) == 317 and phases[0] == Fraction(383, 120323)
        assert all(abs(float(p) * math.pi / a - 1) < 1e-7 for p, a in zip(
            phases, [0.01] * 314 + [0.001, 0.002, 0.003]
        ))
        quizx.full_simp(graph)
        output = quizx.extract_circuit(graph)
        assert output.qubits == 1
        return output.to_matrix(), len(output.gates), zx_gates(output)
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
        return circuit.get_unitary(), circuit.n_gates, tket_gates(circuit)
    if tool == "bqskit":
        from bqskit.compiler import Compiler
        from bqskit.ir.lang.qasm2 import OPENQASM2Language
        from bqskit.passes import (
            ForEachBlockPass,
            QuickPartitioner,
            ScanningGateRemovalPass,
            UnfoldPass,
        )

        circuit = OPENQASM2Language().decode(source)
        # Check that import preserves every angle before optimization.
        angles = [op.params for op in circuit if op.gate.name == "RZGate"]
        assert angles == [[0.00028]] * 11220
        with Compiler(num_workers=1) as compiler:
            output = compiler.compile(
                circuit,
                [
                    QuickPartitioner(3),
                    ForEachBlockPass(ScanningGateRemovalPass()),
                    UnfoldPass(),
                ],
            )
        return output.get_unitary(), output.num_operations, bqskit_gates(output)
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
        return np.asarray(output.get_unitary(), dtype=np.complex128), len(gates), qrisp_gates(output)
    if tool == "qbraid":
        import warnings

        from qbraid.runtime.ionq.device import IonQDevice
        from qbraid.runtime.ionq.provider import IonQProvider
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Operator

        # Check the input independently before testing qBraid.
        parsed = qiskit_load(source)
        assert [i.operation.params[0] for i in parsed.data] == (
            [1 / 408 * math.pi] * 454 + [-math.pi / 408] * 46
        )
        # IonQDevice.run applies transform() and then prepare() before it submits
        # the job. Run these two steps offline on the IonQ simulator profile.
        warnings.filterwarnings("ignore", module="qbraid")
        provider = IonQProvider(api_key="offline")
        device = IonQDevice(
            provider._build_profile({"backend": "simulator", "qubits": 29}),
            provider.session,
        )
        output = device.prepare(device.transform(source))
        # Simulate the emitted IonQ rotations directly, without an optimizer.
        # This fixture emits only RZ gates; reject anything else.
        circuit = QuantumCircuit(output["qubits"])
        for gate in output["circuit"]:
            assert set(gate) == {"gate", "target", "rotation"}
            assert gate["gate"] == "rz"
            circuit.rz(gate["rotation"], gate["target"])
        return Operator(circuit).data, len(output["circuit"]), ionq_gates(output)
    raise ValueError(f"Unknown tool: {tool}")


def main():
    print("Prepare every qubit in |+>; measure X on q[0]. Compare P(+1).")
    reproduced_all = True
    selected = sys.argv[1:]
    for tool, (source, description, expected) in examples().items():
        if selected and tool not in selected:
            continue
        print(f"\n{tool}\n  Input:    {description}", flush=True)
        output, count, gates = optimize(tool, source)
        ideal = ideal_probability(source)
        actual = output_probability(output)
        print(f"  Expected: {expected}; output: {count} gate{'' if count == 1 else 's'}")
        print(
            f"  Measure X: expected {float(100 * ideal):.6f}% -> output {100 * actual:.6f}%"
        )
        # Two independent checks. The reference simulates the output gates without
        # the tool's matrix code. The expected value also holds for the exact
        # binary64 angles that the tool gets.
        simulated = simulate(qubit_count(source), gates)
        rounded = ideal_probability(float_inputs(source))
        print(
            f"  Checks:   output gates in the reference {float(100 * simulated):.6f}%; "
            f"expected with float angles {float(100 * rounded):.6f}%"
        )
        error = float(100 * abs(ideal - actual))
        checks = [abs(ideal - simulated), abs(rounded - simulated)]
        reproduced = error > 99 and all(float(100 * e) > 99 for e in checks)
        reproduced_all &= reproduced
        print(f"  Error:    {error:.6f} percentage points")
        print(f"  Reproduced (>99 pp): {'yes' if reproduced else 'no'}")
    if not reproduced_all:
        raise SystemExit("At least one reported failure did not reproduce.")


if __name__ == "__main__":
    main()
