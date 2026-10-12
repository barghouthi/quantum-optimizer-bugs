#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["qiskit==2.5.1"]
# ///
"""A quilc 1.26.0 angle bug. Compile the OpenQASM 2 input with quilc and print the output Quil.

quilc reads and writes Quil. The script parses the input with Qiskit, writes each
rz as a Quil RZ with the exact float angle, and runs quilc in Docker.

Usage: uv run quilc_bugs.py
"""
import subprocess

from qiskit import qasm2

# quilc 1.26.0, the newest official Docker image, pinned by digest.
QUILC_IMAGE = "rigetti/quilc@sha256:0e65aefeedd467f6dfb4727fc12865f4ab9316d453602a9f7f097ba3b749a424"
HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def quilc(source):
    quil = ""
    for i in qasm2.loads(source).data:
        assert i.operation.name == "rz"
        quil += f"RZ({i.operation.params[0]!r}) 0\n"  # repr keeps the exact float
    result = subprocess.run(["docker", "run", "--rm", "-i", QUILC_IMAGE], input=quil, capture_output=True, text=True, check=True)
    return f"# Input Quil:\n{quil}# Output Quil:\n{result.stdout}"


EXAMPLES = {
    # Large angles: the RZ merge rule rounds 1e17 + pi back to 1e17.
    # The full-rotation rule also removes each 1e17 rotation on its own.
    "large angles": (
        HEADER + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
        "RZ(pi), equivalent to Z",
        quilc,
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
