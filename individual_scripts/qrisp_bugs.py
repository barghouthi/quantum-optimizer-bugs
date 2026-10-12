#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["qrisp==0.9.9"]
# ///
"""Two Qrisp 0.9.9 compiler bugs. Compile each OpenQASM 2 input and print the output QASM.

Usage: uv run qrisp_bugs.py
"""
from qrisp import QuantumCircuit, QuantumVariable

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def compile_session(source):
    """Compile with Qrisp's default session compiler and return OpenQASM 2."""
    circuit = QuantumCircuit.from_qasm_str(source)
    variable = QuantumVariable(circuit.num_qubits())
    variable.qs.append(circuit.to_gate(), list(variable))
    return variable.qs.compile().qasm()


EXAMPLES = {
    # Large angles: 1e17 + pi rounds back to 1e17, so the pi rotation is lost.
    "large angles": (
        HEADER + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
        "rz(pi), equivalent to Z",
        compile_session,
    ),
    # Small angles: each pair u3(7e-5,0,0) u3(7e-5,0,0) passes an inverse check
    # and cancels. The 44,880 rotations add up to approximately pi.
    "small angles": (
        HEADER + "qreg q[1];\n" + "u3(0.00007,0,0) q[0];\n" * 44880,
        "ry(3.1416), almost Y",
        compile_session,
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
