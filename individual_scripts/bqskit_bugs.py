#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["bqskit==1.2.1", "bqskitrs==0.4.1", "mpmath==1.3.0", "sympy==1.13.0"]
# ///
"""Two BQSKit 1.2.1 bugs. Compile each OpenQASM 2 input and print the output.

The first output contains raw unitary blocks, which OpenQASM 2 cannot express.
The script prints those blocks as matrices.

Usage: uv run bqskit_bugs.py
"""
import mpmath as mp
import numpy as np
from bqskit import compile
from bqskit.compiler import Compiler
from bqskit.ir.lang.qasm2 import OPENQASM2Language
from bqskit.passes import BlockConversionPass
from bqskit.passes.processing.extract_diagonal import ExtractDiagonalPass

from equivalence import check  # exact 100-digit check, in this folder

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def extract_diagonal(source):
    circuit = OPENQASM2Language().decode(source)
    with Compiler(num_workers=1) as compiler:
        output = compiler.compile(
            circuit,
            [
                BlockConversionPass("variable"),
                # Seed 75 gives the largest error. Other seeds give errors from 0 to 99 pp.
                ExtractDiagonalPass(success_threshold=1e-12, instantiate_options={"seed": 75, "multistarts": 1}),
            ],
        )
    lines, gates = [], []  # printed text, and the gate list for the exact check
    for op in output:
        qubits = ",".join(f"q[{w}]" for w in op.location)
        if op.gate.name.startswith("VariableUnitaryGate"):
            unitary = np.asarray(op.get_unitary())  # BQSKit's definition of the gate
            lines.append(f"// unitary on {qubits}:")
            for row in unitary:
                lines.append("//   " + "  ".join(f"{z.real:+.6f}{z.imag:+.6f}j" for z in row))
            gates.append(("unitary", [[mp.mpc(complex(z)) for z in row] for row in unitary], list(op.location)))
        else:
            assert op.gate.qasm_name == "cx", op.gate
            lines.append(f"cx {qubits};")
            gates.append(("cx", [], list(op.location)))
    return "\n".join(lines), gates


def compile_level_2(source):
    circuit = OPENQASM2Language().decode(source)
    with Compiler(num_workers=1) as compiler:
        output = compile(circuit, optimization_level=2, synthesis_epsilon=1e-3, seed=1, compiler=compiler)
    return output.to("qasm")


EXAMPLES = {
    # Wrong qubit: ExtractDiagonalPass merges a diagonal from the second block
    # into the first without remapping their different qubit orders.
    "wrong qubit": (
        HEADER + "gate g(a) x,y { rz(a) x; id y; }\nqreg q[2];\ng(0) q[1],q[0];\ng(pi) q[0],q[1];\n",
        "Z on q[0]",
        extract_diagonal,
    ),
    # Small angles: compile passes synthesis_epsilon, documented as a maximum
    # distance, to a gate-removal pass that compares it with a squared cost.
    # Each block loses its rz(pi/36). The 36 rotations add up to pi.
    "small angles": (
        HEADER + "qreg q[4];\n" + "rz(pi/36) q[0];\ncx q[0],q[1];\ncx q[0],q[2];\ncx q[0],q[3];\n" * 36,
        "rz(pi) on q[0], equivalent to Z",
        compile_level_2,
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
