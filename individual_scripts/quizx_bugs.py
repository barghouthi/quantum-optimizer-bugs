#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["quizx==0.3.0", "pyzx==0.9.0"]
# ///
"""Two QuiZX 0.3.0 angle bugs. Simplify each OpenQASM 2 input and print the output QASM.

Usage: uv run quizx_bugs.py
"""
import quizx

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def full_simp(source):
    graph = quizx.qasm(source)
    quizx.full_simp(graph)
    return quizx.extract_circuit(graph).to_qasm()


EXAMPLES = {
    # Exact rational angles: QuiZX adds phases as 64-bit rationals of pi. The
    # product of the two coprime denominators is 2**64 - 1, so the addition overflows.
    "rational angles": (
        HEADER + "qreg q[1];\nrz(2147483647*pi/4294967295) q[0];\nrz(2147483648*pi/4294967297) q[0];\n",
        "rz(~pi), almost Z",
        full_simp,
    ),
    # Small decimal angles: import gives each decimal angle a large denominator,
    # for example 383/120323 for 0.01. Adding the phases overflows in the same way.
    "small angles": (
        HEADER + "qreg q[1];\n" + "rz(0.01) q[0];\n" * 314 + "rz(0.001) q[0];\nrz(0.002) q[0];\nrz(0.003) q[0];\n",
        "rz(3.146), almost Z",
        full_simp,
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
