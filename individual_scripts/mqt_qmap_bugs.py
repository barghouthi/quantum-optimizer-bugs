#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["mqt.qmap==3.10.0", "mqt.core==3.10.1", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""An MQT QMAP 3.10.0 bug. Map the OpenQASM 2 input and print the output QASM.

Usage: uv run mqt_qmap_bugs.py
"""
from mqt.core import load
from mqt.qmap import sc

from equivalence import check  # exact 100-digit check, in this folder

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def map_heuristic(source):
    config = sc.Configuration()
    config.add_measurements_to_mapped_circuit = False
    config.add_barriers_between_layers = False
    output, _ = sc.map_(load(source), sc.Architecture(2, {(0, 1), (1, 0)}), config)
    return output.qasm2_str()


EXAMPLES = {
    # The heuristic mapper emits a CNOT for every controlled gate.
    # crz(2*pi) is Z on its control, but the output is cx.
    "controlled rotation": (
        HEADER + "qreg q[2];\ncrz(2*pi) q[0],q[1];\n",
        "Z on q[0]",
        map_heuristic,
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
