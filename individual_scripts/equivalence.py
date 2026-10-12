#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""Exact equivalence check for small circuits. It does not use any tool under test.

The check reads every angle exactly: 0.01 is 1/100, and pi is symbolic.
It computes the unitary of each circuit with 100-digit arithmetic and compares
the two unitaries up to global phase.

The scripts in this folder import this module. Run it directly to test the
check on circuit pairs with known answers: uv run equivalence.py
"""
import re

import mpmath as mp

mp.mp.dps = 100
# Largest distance that still counts as equivalent. A missing rotation by theta
# gives a distance of about theta**2 / 8, so this catches rotations above 3e-10.
# Float rounding of printed output angles gives distances near 1e-32.
TOLERANCE = mp.mpf("1e-20")


def matrix(name, angles):
    """100-digit matrix of a gate. The first operand is the most significant qubit."""
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
        return [[c, -mp.expj(lam) * s], [mp.expj(phi) * s, mp.expj(phi + lam) * c]]

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
    name = {"u": "u3", "U": "u3", "p": "u1", "CX": "cx", "cp": "cu1", "cnot": "cx"}.get(name, name)
    if name in single:
        return 1, single[name](*angles)
    if name in ("cx", "cy", "cz", "crz", "cu1"):
        return 2, controlled(single[name[1:]](*angles))
    raise ValueError(f"Unsupported gate: {name}")


def angle(expression, cache={}):
    """Exact 100-digit value of an angle expression such as 2*pi/44430 or 1e-5."""
    if expression not in cache:
        import sympy as sp

        value = sp.sympify(expression, locals={"pi": sp.pi}, rational=True)
        cache[expression] = mp.mpf(str(value.evalf(110)))
    return cache[expression]


def parse(qasm):
    """Read OpenQASM 2, or the OpenQASM 3 subset used here, into (qubits, gates).

    Each gate is (name, angles, wires). Registers get global qubit indices in
    declaration order. Unsupported statements raise an error.
    """
    qasm = re.sub(r"//[^\n]*", "", qasm)
    macros = {}
    for match in re.finditer(r"gate (\w+)(?:\(([^)]*)\))? ([^{]*)\{([^}]*)\}", qasm):
        params = [p.strip() for p in (match[2] or "").split(",") if p.strip()]
        args = [a.strip() for a in match[3].split(",")]
        macros[match[1]] = (params, args, match[4])
    qasm = re.sub(r"gate \w+[^{]*\{[^}]*\}", "", qasm)

    registers, qubits, gates = {}, 0, []

    def wires(text):
        found = re.findall(r"(\w+)\[(\d+)\]", text)
        if not found:
            raise ValueError(f"Unsupported operands: {text}")
        return [registers[name] + int(index) for name, index in found]

    def apply(statement, bindings=None):
        gate = re.fullmatch(r"(\w+)(?:\((.*)\))? (.+)", statement)
        if not gate:
            raise ValueError(f"Unsupported statement: {statement}")
        name, params, operands = gate[1], gate[2], gate[3]
        if bindings:  # inside a macro body: substitute its parameters and operands
            substitute = lambda text: re.sub(r"\b\w+\b", lambda m: bindings.get(m[0], m[0]), text)
            params = substitute(params) if params else params
            operands = substitute(operands)
        values = [p.strip() for p in params.split(",")] if params else []
        if name in macros:
            names, args, body = macros[name]
            inner = dict(zip(names, (f"({v})" for v in values)))
            inner.update(zip(args, (o.strip() for o in operands.split(","))))
            for line in body.split(";"):
                if line.strip():
                    apply(" ".join(line.split()), inner)
            return
        gates.append((name, [angle(v) for v in values], wires(operands)))

    for statement in qasm.split(";"):
        statement = " ".join(statement.split())
        if not statement or re.fullmatch(r'OPENQASM [23](\.0)?|include "[\w.]+"|barrier .*', statement):
            continue
        register = re.fullmatch(r"qreg (\w+)\[(\d+)\]|qubit\[(\d+)\] (\w+)", statement)
        if register:
            name, size = (register[1], register[2]) if register[1] else (register[4], register[3])
            registers[name] = qubits
            qubits += int(size)
            continue
        # ctrl @ pow(n) @ gphase(a) applies the phase n*a when the qubit is 1.
        controlled = re.fullmatch(r"ctrl @ pow\((\d+)\) @ gphase\((.+)\) (.+)", statement)
        if controlled:
            gates.append(("u1", [int(controlled[1]) * angle(controlled[2])], wires(controlled[3])))
            continue
        apply(statement)
    return qubits, gates


def unitary(qubits, gates):
    """100-digit unitary, stored as a list of columns. Qubit 0 is the most significant."""
    dimension = 1 << qubits
    columns = [[mp.mpc(int(r == c)) for r in range(dimension)] for c in range(dimension)]
    cache = {}
    for name, angles, wires in gates:
        if name == "unitary":
            m = angles
        else:
            key = (name, tuple(angles))
            if key not in cache:
                cache[key] = matrix(name, list(angles))[1]
            m = cache[key]
        arity = len(wires)
        if len(m) != 1 << arity or len(set(wires)) != arity or max(wires) >= qubits:
            raise ValueError(f"Wrong gate wires: {name} {wires}")
        entries = [[(c, m[r][c]) for c in range(1 << arity) if m[r][c] != 0] for r in range(1 << arity)]
        shifts = [qubits - 1 - w for w in wires]
        groups = []
        for basis in range(dimension):
            if not any((basis >> s) & 1 for s in shifts):
                groups.append([
                    basis | sum(((k >> (arity - 1 - j)) & 1) << shifts[j] for j in range(arity))
                    for k in range(1 << arity)
                ])
        for state in columns:
            for group in groups:
                old = [state[g] for g in group]
                for r, out in enumerate(group):
                    state[out] = sum((value * old[c] for c, value in entries[r]), mp.mpc(0))
    return columns


def distance(qubits, first, second):
    """1 - |Tr(U^dagger V)| / d. It is 0 exactly when U and V differ only by global phase."""
    u, v = unitary(qubits, first), unitary(qubits, second)
    trace = sum(mp.conj(a) * b for cu, cv in zip(u, v) for a, b in zip(cu, cv))
    return 1 - abs(trace) / (1 << qubits)


def check(input_qasm, output, output_qubits=None):
    """Compare an input circuit with an output circuit and return a one-line verdict.

    `output` is OpenQASM text or a gate list. Missing output qubits count as idle.
    """
    qubits, first = parse(input_qasm)
    if isinstance(output, str):
        output_qubits, second = parse(output)
    else:
        second = output
    if output_qubits is not None and output_qubits > qubits:
        raise ValueError("The output has more qubits than the input")
    value = distance(qubits, first, second)
    verdict = "equivalent" if value <= TOLERANCE else "NOT equivalent"
    return f"// Exact check: {verdict} (distance {mp.nstr(value, 3)}; 0 means equal up to global phase)"


def self_test():
    """Check known pairs. Each pair must give the expected verdict."""
    header = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'
    one = header + "qreg q[1];\n"
    two = header + "qreg q[2];\n"
    large = one + "rz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n"
    tket = one + "rz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n"
    core4 = 'OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[1] q;\nctrl @ pow(80599999999716496) @ gphase(1) q[0];\n'
    macro = header + "gate g(a) x,y { rz(a) x; id y; }\nqreg q[2];\ng(0) q[1],q[0];\ng(pi) q[0],q[1];\n"
    cases = [
        ("1e17 + pi - 1e17 equals rz(pi)", large, one + "rz(pi) q[0];", True),
        ("1e17 + pi - 1e17 is not identity", large, one, False),
        ("rz(pi); rz(-1e16*pi) equals rz(pi)", tket, one + "rz(pi) q[0];", True),
        ("crz(2*pi) equals Z on the control", two + "crz(2*pi) q[0],q[1];", two + "z q[0];", True),
        ("crz(2*pi) is not cx", two + "crz(2*pi) q[0],q[1];", two + "cx q[0], q[1];", False),
        ("ctrl @ pow(K) @ gphase(1) equals p(K)", core4, one + "p(80599999999716496) q[0];", True),
        ("p(K) is not p(-3.1413486180563552)", core4, one + "p(-3.1413486180563552) q[0];", False),
        ("macro g(pi) on q[0] equals Z on q[0]", macro, two + "z q[0];", True),
        ("macro g(pi) on q[0] is not Z on q[1]", macro, two + "z q[1];", False),
        ("other register names", one + "rz(0.3) q[0];", header + "qreg v[1];\nrz(0.3) v[0];", True),
        ("h equals u(pi/2,0,pi)", one + "h q[0];", one + "u(pi/2,0,pi) q[0];", True),
        ("cx direction matters", two + "cx q[0],q[1];", two + "cx q[1],q[0];", False),
        ("cx equals h cz h", two + "cx q[0],q[1];", two + "h q[1];\ncz q[0],q[1];\nh q[1];", True),
        ("rz(1e-5) is not identity", one + "rz(0.00001) q[0];", one, False),
    ]
    failures = 0
    for label, first, second, equivalent in cases:
        verdict = check(first, second)
        ok = ("NOT" not in verdict) == equivalent
        failures += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {label}: {verdict.removeprefix('// Exact check: ')}")
    if failures:
        raise SystemExit(f"{failures} self-test case(s) failed.")


if __name__ == "__main__":
    self_test()
