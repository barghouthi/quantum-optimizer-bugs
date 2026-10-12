#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["cirq-core==1.7.0", "ply==3.11", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""Two Cirq 1.7.0 angle bugs. Transform each OpenQASM 2 input and print the output QASM.

Usage: uv run cirq_bugs.py
"""
import cirq
from cirq.contrib.qasm_import import circuit_from_qasm

from equivalence import check  # exact 100-digit check, in this folder

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def to_qasm(output, circuit):
    # Cirq writes each angle as a multiple of pi, rounded to `precision` decimal
    # places (default 10). 40 places keep all digits of these multiples.
    return output.to_qasm(precision=40, qubit_order=sorted(circuit.all_qubits()))


def eject_z(source):
    circuit = circuit_from_qasm(source)
    return to_qasm(cirq.eject_z(circuit, atol=0), circuit)


def insertion_sort(source):
    circuit = circuit_from_qasm(source)
    return to_qasm(cirq.transformers.insertion_sort_transformer(circuit), circuit)


EXAMPLES = {
    # Large angles: eject_z tracks the phase as exponent/2 and rounds 5e15 + 0.5 to 5e15.
    "large angles": (
        HEADER + "qreg q[1];\nrz(10000000000000000*pi) q[0];\nrz(pi) q[0];\nrz(-10000000000000000*pi) q[0];\n",
        "rz(pi), equivalent to Z",
        eject_z,
    ),
    # Small angles: the commutation check uses a relative tolerance, so rz(1e-5)
    # on the target counts as commuting with cx. The two cx gates become adjacent
    # and cancel, which leaves the 300,000 rotations on q[1] instead of q[0].
    "small angles": (
        HEADER + "qreg q[2];\nh q[1];\ncx q[0],q[1];\n" + "rz(0.00001) q[1];\n" * 300000 + "cx q[0],q[1];\n",
        "rz(3) on q[0]",
        insertion_sort,
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
