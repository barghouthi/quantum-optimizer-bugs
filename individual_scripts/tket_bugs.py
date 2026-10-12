#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["pytket==2.18.5"]
# ///
"""Two tket (pytket 2.18.5) angle bugs. Optimize each OpenQASM 2 input and print the output QASM.

Usage: uv run tket_bugs.py
"""
from pytket.passes import FullPeepholeOptimise, RemoveRedundancies
from pytket.qasm import circuit_from_qasm_str, circuit_to_qasm_str

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def remove_redundancies(source):
    circuit = circuit_from_qasm_str(source)
    RemoveRedundancies().apply(circuit)
    return circuit_to_qasm_str(circuit)


def full_peephole(source):
    circuit = circuit_from_qasm_str(source)
    FullPeepholeOptimise().apply(circuit)
    return circuit_to_qasm_str(circuit)


BLOCK = (
    "rz(6.134e-10) q[0];\nrx(1.064e-10) q[0];\nry(6.317e-11) q[0];\nz q[0];\n"
    "rx(8.644e-11) q[0];\ny q[0];\ncx q[0],q[1];\nx q[0];\ncx q[0],q[1];\n"
)

EXAMPLES = {
    # Large angles: tket keeps half-turns 1 and -1e16. 1 + (-1e16) rounds to
    # -1e16, a full rotation, so the pass removes it.
    "large angles": (
        HEADER + "qreg q[1];\nrz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n",
        "rz(pi), equivalent to Z",
        remove_redundancies,
    ),
    # Small angles: a shortcut divides two values that are almost zero when it
    # squashes each chain of single-qubit gates.
    "small angles": (
        HEADER + "qreg q[2];\ns q[0];\n" + BLOCK * 5 + "s q[0];\n",
        "Z on q[0]",
        full_peephole,
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
