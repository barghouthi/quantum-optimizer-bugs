#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["mqt.core==4.0.0", "qiskit==2.5.1", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""An MQT Core 4.0.0 angle bug. Compile the OpenQASM 3 input and print the optimized IR.

MQT Core 4 compiles OpenQASM 3 to MLIR. The script prints the optimized MLIR and
the same circuit as OpenQASM 2, through MQT's own conversion to Qiskit.
MQT Core 4 needs its own environment, separate from MQT QMAP 3.10.

Usage: uv run mqt_core4_bugs.py
"""
from mqt.core import mlir
from qiskit import qasm2

from equivalence import check  # exact 100-digit check, in this folder


def optimize(source):
    result = mlir.compile_program(source, output=mlir.OutputFormat.QCO_OPTIMIZED)
    qasm = qasm2.dumps(result.to_qiskit())
    return f"// Optimized MLIR:\n{result.ir}\n// As OpenQASM 2:\n{qasm}", qasm


EXAMPLES = {
    # Large angles: the global-phase normalization multiplies the phase by K and
    # reduces it modulo the binary64 value of 2*pi. The small error in 2*pi adds up
    # to about pi over 1.3e16 periods. The pass skips MQT's 1e4-radian angle limit.
    "large angles": (
        'OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[1] q;\nctrl @ pow(80599999999716496) @ gphase(1) q[0];\n',
        "p(80599999999716496), almost identity",
        optimize,
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
