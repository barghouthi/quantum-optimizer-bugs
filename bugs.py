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
"""QASM reproducers. Prepare every qubit in |+>; measure X on q[0]."""

import math
import re
import sys
from fractions import Fraction

import mpmath as mp

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
mp.mp.dps = 100
# quilc 1.26.0, the newest official Docker image, pinned by digest.
QUILC_IMAGE = "rigetti/quilc@sha256:0e65aefeedd467f6dfb4727fc12865f4ab9316d453602a9f7f097ba3b749a424"
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
        "qrisp": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "rz(1e17); rz(pi); rz(-1e17)",
            "Z",
        ),
        "qbraid": (
            HEADER
            + "qreg q[1];\nu3(0,100000000000000000,0) q[0];\n"
            + "rz(-100000000000000000) q[0];\n",
            "u3(0,1e17,0); rz(-1e17)",
            "identity (up to global phase)",
        ),
        "cudaq": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "rz(1e17); rz(pi); rz(-1e17)",
            "Z",
        ),
        "quilc": (
            HEADER
            + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
            "rz(1e17); rz(pi); rz(-1e17)",
            "Z",
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
    return int(re.search(r"(?:qreg q\[|qubit\[)(\d+)\]", source)[1])


def ideal_probability(source):
    """Exact diagonal phases from fixture QASM, then 100-digit probability.

    This small reference supports rz, crz, u3(0,phi,lambda), id, ctrl @ pow(n) @ gphase(a),
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
        gate = re.fullmatch(r"(rz|crz|u3)\((.*)\) (q\[\d+\](?:,q\[\d+\])?)", statement)
        if not gate:
            raise ValueError(f"Unsupported reference gate: {statement}")
        angles = [
            sp.sympify(a, locals={"pi": sp.pi}, rational=True)
            for a in gate[2].split(",")
        ]
        wires = [int(q) for q in re.findall(r"q\[(\d+)\]", gate[3])]
        if len(wires) != (2 if gate[1] == "crz" else 1):
            raise ValueError(f"Wrong gate arity: {statement}")
        if gate[1] == "u3":
            if len(angles) != 3 or angles[0] != 0:
                raise ValueError(f"Unsupported reference gate: {statement}")
        elif len(angles) != 1:
            raise ValueError(f"Wrong parameter count: {statement}")
        for basis in range(len(phases)):
            bits = [(basis >> (qubits - 1 - q)) & 1 for q in wires]
            if gate[1] == "u3":
                phases[basis] += (angles[1] + angles[2]) * bits[0]
            elif gate[1] == "rz" or bits[0]:
                phases[basis] += angles[0] * (bits[-1] - sp.Rational(1, 2))
    half = len(phases) // 2
    probability = (
        1 + sum(sp.cos(phases[i] - phases[i + half]) for i in range(half)) / half
    ) / 2
    return mp.mpf(str(probability.evalf(100)))


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


class Skipped(Exception):
    """A tool that cannot run here, such as quilc without Docker."""


def quil_gates(program):
    """Gates of quilc's native one-qubit output. Other gates raise an error."""
    import sympy as sp

    gates = []
    for line in program.splitlines():
        line = line.split("#")[0].strip()
        if not line or line == "HALT":
            continue
        gate = re.fullmatch(r"(RZ|RX)\((.*)\) 0", line)
        if not gate:
            raise ValueError(f"Unsupported quilc output: {line}")
        angle = sp.sympify(gate[2], locals={"pi": sp.pi}, rational=True).evalf(110)
        gates.append((gate[1].lower(), [mp.mpf(str(angle))], [0]))
    return gates


def quil_unitary(program):
    """Unitary of quilc's native one-qubit output, in binary64."""
    import numpy as np

    unitary = np.eye(2, dtype=complex)
    for name, (angle,), _ in quil_gates(program):
        angle = float(angle)
        c, s = math.cos(angle / 2), math.sin(angle / 2)
        if name == "rz":
            matrix = np.diag([complex(c, -s), complex(c, s)])
        else:
            matrix = np.array([[c, -1j * s], [-1j * s, c]])
        unitary = matrix @ unitary
    return unitary


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


def qir_gates(qir):
    """Gates in CUDA-Q's compiled QIR for one qubit. LLVM may print a double in hex."""
    import struct

    gates = []
    for name, angle in re.findall(r"call void @__quantum__qis__(\w+?)(?:__body)?\((?:double ([^,]+), )?ptr", qir):
        if angle and angle.startswith("0x"):
            angle = struct.unpack(">d", bytes.fromhex(angle[2:]))[0]
        gates.append((name, [exact(float(angle))] if angle else [], [0]))
    return gates


def mqt_gates(circuit):
    names = {"x": "x", "y": "y", "z": "z", "h": "h", "rx": "rx", "ry": "ry", "rz": "rz", "p": "u1"}
    gates = []
    for op in circuit:
        name = op.type_.name
        controls = sorted(c.qubit for c in op.controls)
        if any(c.type_.name != "Pos" for c in op.controls) or len(controls) > 1:
            raise ValueError(f"Unsupported MQT output operation: {op}")
        name = ("c" if controls else "") + names.get(name, name)
        gates.append((name, [exact(a) for a in op.parameter], controls + list(op.targets)))
    return gates


def describe(gates):
    """Readable output gate list, with each angle shown as a float."""
    if not gates:
        return "no gates"
    parts = []
    for name, angles, wires in gates:
        if name == "unitary":
            angles = []
        values = f"({','.join(repr(float(a)) for a in angles)})" if angles else ""
        name = {"u1": "p"}.get(name, name)  # the same gate; p is the newer name
        parts.append(f"{name}{values} " + ",".join(f"q[{w}]" for w in wires))
    return "; ".join(parts)


def optimize(tool, source):
    """Return the output's native matrix, gate count, and gate list. All inputs are parsed QASM.

    CUDA-Q returns its output state and no gate count.
    """
    if tool == "qiskit":
        from qiskit.quantum_info import Operator
        from qiskit.transpiler import PassManager
        from qiskit.transpiler.passes import CommutativeCancellation

        output = PassManager([CommutativeCancellation(basis_gates=["rz"])]).run(
            qiskit_load(source)
        )
        return Operator(output).data, len(output.data), qiskit_gates(output)
    if tool == "pennylane":
        import pennylane as qml

        # Explicit QASM-parser adapter: this fixture contains only rz gates.
        parsed = qiskit_load(source)
        assert all(i.operation.name == "rz" for i in parsed.data)
        tape = qml.tape.QuantumScript(
            [qml.RZ(i.operation.params[0], 0) for i in parsed.data]
        )
        output = qml.transforms.merge_rotations(tape, atol=0)[0][0]
        return qml.matrix(output, wire_order=[0]), len(output.operations), pennylane_gates(output)
    if tool == "cirq":
        import cirq
        from cirq.contrib.qasm_import import circuit_from_qasm

        circuit = circuit_from_qasm(source)
        output = cirq.eject_z(circuit, atol=0)
        return (
            output.unitary(qubit_order=sorted(circuit.all_qubits())),
            len(list(output.all_operations())),
            cirq_gates(output),
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
        return output.get_unitary(), output.num_operations, bqskit_gates(output)
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
        return output.to_matrix(), len(output.gates), zx_gates(output)
    if tool == "tket":
        from pytket.passes import RemoveRedundancies
        from pytket.qasm import circuit_from_qasm_str

        circuit = circuit_from_qasm_str(source)
        RemoveRedundancies().apply(circuit)
        return circuit.get_unitary(), circuit.n_gates, tket_gates(circuit)
    if tool == "qrisp":
        import numpy as np
        from qrisp import QuantumCircuit, QuantumVariable

        circuit = QuantumCircuit.from_qasm_str(source)
        # Check that import preserves all three angles before compiling.
        assert [i.op.params[0] for i in circuit.data] == [1e17, math.pi, -1e17]
        variable = QuantumVariable(1)
        variable.qs.append(circuit.to_gate(), [variable[0]])
        output = variable.qs.compile()
        gates = [i for i in output.data if i.op.name != "qb_alloc"]
        # The compiled unitary is complex64; widen it before computing probabilities.
        return np.asarray(output.get_unitary(), dtype=np.complex128), len(gates), qrisp_gates(output)
    if tool == "qbraid":
        from qbraid.transpiler.conversions.qasm2.qasm2_to_qasm3 import qasm2_to_qasm3
        from qbraid.transpiler.conversions.qasm3.qasm3_to_ionq import qasm3_to_ionq
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Operator

        # Check the input independently before testing qBraid's conversion.
        parsed = qiskit_load(source)
        assert [i.operation.name for i in parsed.data] == ["u3", "rz"]
        assert [list(i.operation.params) for i in parsed.data] == [
            [0, 1e17, 0], [-1e17]
        ]
        assert abs(output_probability(Operator(parsed).data) - 1) < 1e-12
        # Use the SDK's QASM 3 path, which lowers U3 to IonQ's rotation gates.
        # Direct calls avoid discovering unrelated optional backends.
        output = qasm3_to_ionq(qasm2_to_qasm3(source))
        # Simulate the emitted IonQ rotations directly, without an optimizer.
        # This fixture emits only RX and RZ gates; reject anything else.
        circuit = QuantumCircuit(output["qubits"])
        for gate in output["circuit"]:
            assert set(gate) == {"gate", "target", "rotation"}
            assert gate["gate"] in ("rx", "rz")
            getattr(circuit, gate["gate"])(gate["rotation"], gate["target"])
        return Operator(circuit).data, len(output["circuit"]), ionq_gates(output)
    if tool == "cudaq":
        import cudaq
        import numpy as np

        # Explicit QASM-parser adapter: this fixture contains only rz gates.
        parsed = qiskit_load(source)
        assert all(i.operation.name == "rz" for i in parsed.data)
        kernel = cudaq.make_kernel()
        qubit = kernel.qalloc(1)
        kernel.h(qubit[0])  # prepare |+>
        for i in parsed.data:
            kernel.rz(i.operation.params[0], qubit[0])
        # Check that the kernel keeps all three exact angles before optimization.
        assert re.findall(r"arith.constant (\S+) : f64", str(kernel)) == [
            "-1.000000e+17", "3.1415926535897931", "1.000000e+17"
        ]
        # get_state runs the default pipeline. If it removes every gate, the
        # unitary drops the qubit, so return the output state for the |+> input.
        # The compiled QIR gives the output gates. The first gate prepares |+>.
        gates = qir_gates(cudaq.translate(kernel, format="qir"))
        assert gates[0] == ("h", [], [0])
        return np.array(cudaq.get_state(kernel)), None, gates[1:]
    if tool == "quilc":
        import shutil
        import subprocess

        if not shutil.which("docker") or subprocess.run(
            ["docker", "info"], capture_output=True
        ).returncode:
            raise Skipped("Docker is not available")
        # Explicit QASM-parser adapter: this fixture contains only rz gates.
        parsed = qiskit_load(source)
        assert all(i.operation.name == "rz" for i in parsed.data)
        # repr gives each exact 64-bit floating-point angle; quilc parses it back exactly.
        quil = "".join(f"RZ({i.operation.params[0]!r}) 0\n" for i in parsed.data)
        result = subprocess.run(
            ["docker", "run", "--rm", "-i", QUILC_IMAGE],
            input=quil, capture_output=True, text=True,
        )
        if result.returncode:
            raise RuntimeError(result.stderr)
        gates = [
            line for line in result.stdout.splitlines()
            if line.split("#")[0].strip() not in ("", "HALT")
        ]
        return quil_unitary(result.stdout), len(gates), quil_gates(result.stdout)
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
        return matrix[order][:, order], len(output), mqt_gates(output)
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
        matrix, gates = json.loads(result.stdout)
        names = {"p": "u1", "u": "u3"}
        gates = [(names.get(name, name), [exact(a) for a in angles], wires) for name, angles, wires in gates]
        return np.array([[complex(*z) for z in row] for row in matrix]), len(gates), gates
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
    # JSON keeps each float exactly: Python writes the shortest repr that reads back.
    gates = [
        [i.operation.name, [float(a) for a in i.operation.params], [output.find_bit(q).index for q in i.qubits]]
        for i in output.data
    ]
    print(json.dumps([matrix, gates]))


def main():
    print("Prepare every qubit in |+>; measure X on q[0]. Compare P(+1).")
    reproduced_all = True
    for tool, (source, description, expected) in examples().items():
        print(f"\n{tool}\n  Input:    {description}", flush=True)
        try:
            unitary, count, gates = optimize(tool, source)
        except Skipped as reason:
            print(f"  Skipped:  {reason}")
            continue
        ideal = ideal_probability(source)
        actual = output_probability(unitary)
        print(f"  Expected: {expected}")
        print(f"  Output:   {describe(gates)}")
        print(
            f"  Measure X: expected {float(100 * ideal):.6f}% -> output {100 * actual:.6f}%"
        )
        # Two independent checks. The reference simulates the output gates without
        # the tool's matrix code. The expected value also holds for the exact
        # binary64 angles that the tool gets.
        simulated = simulate(qubit_count(source), gates)
        rounded = ideal_probability(float_inputs(source))
        error = float(100 * abs(ideal - actual))
        checks = [abs(ideal - simulated), abs(rounded - simulated)]
        reproduced = error > 99 and all(float(100 * e) > 99 for e in checks)
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
