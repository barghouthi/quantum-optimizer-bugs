#!/usr/bin/env -S uv run --locked --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["mqt.core==4.0.0", "qiskit==2.5.1", "mpmath==1.3.0"]
# ///
"""Additional MQT Core 4 reproducer. Prepare q[0] in |+>; measure X on q[0].

MQT QMAP 3.10.0 requires MQT Core 3.10, so this example has its own script.
"""

import mpmath as mp
import numpy as np
from mqt.core import mlir
from qiskit.quantum_info import Operator

K = 80599999999716496  # an exact 64-bit floating-point integer
SOURCE = f"""OPENQASM 3.0;
include "stdgates.inc";
qubit[1] q;
ctrl @ pow({K}) @ gphase(1) q[0];
"""


def main():
    # ctrl @ pow(K) @ gphase(1) is p(K): it applies phase e^(iK) when q[0] is 1.
    mp.mp.dps = 100
    expected = (1 + mp.cos(mp.mpf(K))) / 2

    imported = mlir.compile_program(SOURCE, output=mlir.OutputFormat.QC_IMPORT).ir
    assert "0x4371E5945E0E12C9" in imported  # the binary64 bits of K

    output = mlir.compile_program(
        SOURCE, output=mlir.OutputFormat.QCO_OPTIMIZED
    ).to_qiskit()
    state = Operator(output).data @ np.array([1, 1]) / np.sqrt(2)
    actual = abs(state[0] + state[1]) ** 2 / 2

    print("mqt core 4")
    print(f"  Input:    ctrl @ pow({K}) @ gphase(1) q[0]  (= p({K}))")
    print(f"  Output:   {[(i.operation.name, i.operation.params) for i in output.data]}")
    print(f"  Measure X: expected {float(100 * expected):.6f}% -> output {100 * actual:.6f}%")
    error = float(100 * abs(expected - actual))
    print(f"  Error:    {error:.6f} percentage points")
    print(f"  Reproduced (>99 pp): {'yes' if error > 99 else 'no'}")
    if error <= 99:
        raise SystemExit("The reported failure did not reproduce.")


if __name__ == "__main__":
    main()
