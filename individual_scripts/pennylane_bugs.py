#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["pennylane==0.45.1", "qiskit==2.5.1"]
# ///
"""Two PennyLane 0.45.1 angle bugs. Transform each OpenQASM 2 input and print the output QASM.

PennyLane reads OpenQASM 2 only through the pennylane-qiskit plugin, which uses
Qiskit's parser. This script uses Qiskit's parser directly.

Usage: uv run pennylane_bugs.py
"""
import pennylane as qml
from qiskit import qasm2
from qiskit.quantum_info import Operator

HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def to_qasm(tape, qubits):
    text = qml.to_openqasm(tape, wires=range(qubits), measure_all=False)
    # PennyLane leaves out the register when the circuit has no gates.
    return text if "qreg" in text else text + f"qreg q[{qubits}];\n"


def merge_rotations(source):
    parsed = qasm2.loads(source)
    tape = qml.tape.QuantumScript([qml.RZ(i.operation.params[0], 0) for i in parsed.data])
    output = qml.transforms.merge_rotations(tape, atol=0)[0][0]
    return to_qasm(output, 1)


def unitary_to_rot(source):
    # Like pennylane-qiskit, convert cu1 to a QubitUnitary.
    parsed = qasm2.loads(source)
    ops = []
    for i in parsed.data:
        wires = [parsed.find_bit(q).index for q in i.qubits]
        if i.operation.name == "h":
            ops.append(qml.Hadamard(wires))
        elif i.operation.name == "x":
            ops.append(qml.PauliX(wires))
        else:
            ops.append(qml.QubitUnitary(Operator(i.operation).reverse_qargs().data, wires))
    output = qml.transforms.unitary_to_rot(qml.tape.QuantumScript(ops))[0][0]
    return to_qasm(output, 2)


EXAMPLES = {
    # Large angles: merge_rotations rounds 1e17 + pi back to 1e17.
    "large angles": (
        HEADER + "qreg q[1];\nrz(100000000000000000) q[0];\nrz(pi) q[0];\nrz(-100000000000000000) q[0];\n",
        "rz(pi), equivalent to Z",
        merge_rotations,
    ),
    # Small angles: a tolerance 400 times larger than intended drops each
    # controlled phase. The 351 cu1(pi/351) gates add up to Z on q[0].
    "small angles": (
        HEADER + "qreg q[2];\nh q[1];\nx q[1];\n" + "cu1(pi/351) q[0],q[1];\n" * 351,
        "Z on q[0]",
        unitary_to_rot,
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
