#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["qbraid==0.12.1", "pyqasm==1.0.4"]
# ///
"""Two qBraid 0.12.1 angle bugs. Convert each OpenQASM 2 input to IonQ JSON and print it.

IonQ JSON is qBraid's output format for IonQ devices. The script prints one gate per line.

Usage: uv run qbraid_bugs.py
"""
import json
import warnings

from qbraid.runtime.ionq.device import IonQDevice
from qbraid.runtime.ionq.provider import IonQProvider
from qbraid.transpiler.conversions.qasm2.qasm2_to_qasm3 import qasm2_to_qasm3
from qbraid.transpiler.conversions.qasm3.qasm3_to_ionq import qasm3_to_ionq

warnings.filterwarnings("ignore", module="qbraid")
HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


def ionq_lines(program):
    lines = [json.dumps({k: v for k, v in program.items() if k != "circuit"})]
    return "\n".join(lines + [json.dumps(gate) for gate in program["circuit"]])


def qasm_to_ionq(source):
    # The QASM 3 path lowers u3 to IonQ rotations with PyQASM.
    return ionq_lines(qasm3_to_ionq(qasm2_to_qasm3(source)))


def ionq_device(source):
    # IonQDevice.run applies transform() and prepare() before it submits the job.
    # Run these two steps offline on the IonQ simulator profile.
    provider = IonQProvider(api_key="offline")
    device = IonQDevice(provider._build_profile({"backend": "simulator", "qubits": 29}), provider.session)
    return ionq_lines(device.prepare(device.transform(source)))


EXAMPLES = {
    # Large angles: PyQASM lowers u3(0,phi,0) with phi + pi, which rounds back to phi.
    "large angles": (
        HEADER + "qreg q[1];\nu3(0,100000000000000000,0) q[0];\nrz(-100000000000000000) q[0];\n",
        "identity (up to global phase)",
        qasm_to_ionq,
    ),
    # Small angles: the angle rewriting ignores operator precedence, so
    # 1/408*pi becomes 1/(408*pi). The correct angles add up to pi.
    "small angles": (
        HEADER + "qreg q[1];\n" + "rz(1/408*pi) q[0];\n" * 454 + "rz(-pi/408) q[0];\n" * 46,
        "rz(pi), equivalent to Z",
        ionq_device,
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
