#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["cuda-quantum-cu13==0.16.0", "qiskit==2.5.1", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""A CUDA-Q 0.16.0 angle bug. Compile the OpenQASM 2 input and print the compiled QIR.

CUDA-Q does not read OpenQASM. The script parses the input with Qiskit and builds
a kernel with the same gates. Its OpenQASM export shows the circuit before
optimization, so the script prints the compiled QIR instead.

Usage: uv run cudaq_bugs.py
"""
import re
import struct
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

import cudaq  # noqa: E402
from qiskit import qasm2  # noqa: E402

from equivalence import check  # exact 100-digit check, in this folder

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
    text = "\n".join(line.strip() for line in qir.splitlines() if "call void @__quantum__qis__" in line)
    # The same gates as OpenQASM 2, for the exact check. LLVM may print a double in hex.
    gates = []
    for name, value in re.findall(r"call void @__quantum__qis__(\w+?)(?:__body)?\((?:double ([^,]+), )?ptr", qir):
        if value.startswith("0x"):
            value = repr(struct.unpack(">d", bytes.fromhex(value[2:]))[0])
        gates.append(f"{name}({value}) q[0];" if value else f"{name} q[0];")
    assert gates[0] == "h q[0];"  # the |+> preparation; the check compares the rest
    return text, HEADER + "qreg q[1];\n" + "\n".join(gates[1:])


EXAMPLES = {
    # Large angles: the default pipeline rounds 1e17 + pi back to 1e17.
    "large angles": (
        HEADER + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
        "h (prepares |+>), then rz(pi); the check compares the gates after h",
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
    print("\n".join(lines))


for name, (source, expected, optimize) in EXAMPLES.items():
    print(f"// {name}: {gate_count(source)} input gates; expected output {expected}", flush=True)
    output = optimize(source)
    # optimize() returns the output text, or (text, circuit to check) when the
    # tool's output format is not OpenQASM.
    text, circuit = output if isinstance(output, tuple) else (output, output)
    show(text)
    print(check(source, circuit) + "\n", flush=True)
