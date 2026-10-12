#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["qiskit==2.5.1", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""A quilc 1.26.0 angle bug. Compile the OpenQASM 2 input with quilc and print the output Quil.

quilc reads and writes Quil. The script parses the input with Qiskit, writes each
rz as a Quil RZ with the exact float angle, and runs quilc in Docker.

Usage: uv run quilc_bugs.py
"""
import re
import subprocess

from qiskit import qasm2

from equivalence import check  # exact 100-digit check, in this folder

# quilc 1.26.0, the newest official Docker image, pinned by digest.
QUILC_IMAGE = "rigetti/quilc@sha256:0e65aefeedd467f6dfb4727fc12865f4ab9316d453602a9f7f097ba3b749a424"
HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def quilc(source):
    quil = ""
    for i in qasm2.loads(source).data:
        assert i.operation.name == "rz"
        quil += f"RZ({i.operation.params[0]!r}) 0\n"  # repr keeps the exact float
    result = subprocess.run(["docker", "run", "--rm", "-i", QUILC_IMAGE], input=quil, capture_output=True, text=True, check=True)
    # The same gates as OpenQASM 2, for the exact check.
    gates = []
    for line in result.stdout.splitlines():
        line = line.split("#")[0].strip()
        if line and line != "HALT":
            gate = re.fullmatch(r"(RX|RY|RZ)\((.*)\) (\d+)", line)
            assert gate, f"Unsupported quilc output: {line}"
            gates.append(f"{gate[1].lower()}({gate[2]}) q[{gate[3]}];")
    return f"# Input Quil:\n{quil}# Output Quil:\n{result.stdout}", HEADER + "qreg q[1];\n" + "\n".join(gates)


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
    print("\n".join(lines))


for name, (source, expected, optimize) in EXAMPLES.items():
    print(f"// {name}: {gate_count(source)} input gates; expected output {expected}", flush=True)
    output = optimize(source)
    # optimize() returns the output text, or (text, circuit to check) when the
    # tool's output format is not OpenQASM.
    text, circuit = output if isinstance(output, tuple) else (output, output)
    show(text)
    print(check(source, circuit) + "\n", flush=True)
