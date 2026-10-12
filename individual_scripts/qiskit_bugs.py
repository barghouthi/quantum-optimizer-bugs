#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["qiskit==2.5.1"]
# ///
"""Two Qiskit 2.5.1 angle bugs. Transpile each OpenQASM 2 input and print the output QASM.

Usage: uv run qiskit_bugs.py
"""
from qiskit import qasm2, transpile

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def default_transpile(source):
    """transpile() with default settings (optimization level 2)."""
    return qasm2.dumps(transpile(qasm2.loads(source)))


def transpile_to_cx_u(source):
    return qasm2.dumps(transpile(qasm2.loads(source), basis_gates=["cx", "u"]))


EXAMPLES = {
    # Large angles: CommutativeCancellation rounds 1e17 + pi back to 1e17.
    "large angles": (
        HEADER + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
        "rz(pi), equivalent to Z",
        default_transpile,
    ),
    # Small angles: two-qubit synthesis drops each small entangling crz component.
    # The 44,430 crz gates on control q[0] add up to Z on q[0].
    "small angles": (
        HEADER
        + "qreg q[3];\nh q[1];\nh q[2];\n"
        + "".join(f"crz(2*pi/44430) q[0],q[{1 + k % 2}];\n" for k in range(44430)),
        "Z on q[0]",
        transpile_to_cx_u,
    ),
}


def gate_count(source):
    """Count the gate statements after the header (one statement per line)."""
    skip = ("OPENQASM", "include", "qreg", "qubit", "gate ")
    return sum(1 for line in source.splitlines() if line and not line.startswith(skip))


def show(text, limit=40):
    """Print the output. Cut long outputs after `limit` lines."""
    lines = text.strip().splitlines()
    if len(lines) > limit:
        lines = lines[:limit] + [f"// ... {len(lines) - limit} more lines"]
    print("\n".join(lines) + "\n")


for name, (source, expected, optimize) in EXAMPLES.items():
    print(f"// {name}: {gate_count(source)} input gates; expected output {expected}", flush=True)
    show(optimize(source))
