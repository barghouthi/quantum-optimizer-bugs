#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["cuda-quantum-cu13==0.16.0", "qiskit==2.5.1"]
# ///
"""A CUDA-Q 0.16.0 angle bug. Compile the OpenQASM 2 input and print the compiled QIR.

CUDA-Q does not read OpenQASM. The script parses the input with Qiskit and builds
a kernel with the same gates. Its OpenQASM export shows the circuit before
optimization, so the script prints the compiled QIR instead.

Usage: uv run cudaq_bugs.py
"""
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

import cudaq  # noqa: E402
from qiskit import qasm2  # noqa: E402

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def compile_kernel(source):
    kernel = cudaq.make_kernel()
    qubit = kernel.qalloc(1)
    # Prepare |+>. Without it, rz on |0> is only a global phase.
    kernel.h(qubit[0])
    for i in qasm2.loads(source).data:
        assert i.operation.name == "rz"
        kernel.rz(i.operation.params[0], qubit[0])
    qir = cudaq.translate(kernel, format="qir")
    return "\n".join(line.strip() for line in qir.splitlines() if "call void @__quantum__qis__" in line)


EXAMPLES = {
    # Large angles: the default pipeline rounds 1e17 + pi back to 1e17.
    "large angles": (
        HEADER + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
        "h, then rz(pi)",
        compile_kernel,
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
