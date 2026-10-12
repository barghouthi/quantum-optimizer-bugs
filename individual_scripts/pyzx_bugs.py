#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["pyzx==0.9.0", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""A PyZX 0.9.0 small-angle bug. Optimize the OpenQASM 2 input and print the output QASM.

Usage: uv run pyzx_bugs.py
"""
import pyzx as zx

from equivalence import check  # exact 100-digit check, in this folder

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def teleport_reduce(source):
    circuit = zx.Circuit.from_qasm(source)
    return zx.Circuit.from_graph(zx.teleport_reduce(circuit.to_graph())).to_qasm()


EXAMPLES = {
    # Small angles: phase teleportation fuses the 155 rotations into one phase.
    # The x gates make it negative, but the pass loses that sign.
    "small angles": (
        HEADER
        + "qreg q[1];\ns q[0];\nx q[0];\n"
        + "rz(0.01) q[0];\n" * 155
        + "x q[0];\nrz(0.01) q[0];\ns q[0];\nx q[0];\nsdg q[0];\nrz(0.01) q[0];\n",
        "+1 probability 0.010812% for |+> input, X measurement",
        teleport_reduce,
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
    print("\n".join(lines))


for name, (source, expected, optimize) in EXAMPLES.items():
    print(f"// {name}: {gate_count(source)} input gates; expected output {expected}", flush=True)
    output = optimize(source)
    # optimize() returns the output text, or (text, circuit to check) when the
    # tool's output format is not OpenQASM.
    text, circuit = output if isinstance(output, tuple) else (output, output)
    show(text)
    print(check(source, circuit) + "\n", flush=True)
