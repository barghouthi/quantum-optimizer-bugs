# quantum-optimizer-bugs

This script runs six OpenQASM circuits through quantum circuit optimizers.
Each optimizer removes a rotation that affects the measurement result.
The errors range from **99.5 to 100 percentage points**.

Run all six examples:

```sh
uv run bugs.py
```

## How bad are the outputs?

For each example, prepare `|+>`.
Apply the circuit. Then measure in the X basis.
Compare the expected probability of the **+1** outcome with the output probability.

The error is the absolute probability difference, in **percentage points (pp)**:

| Tool / tested version | Expected +1 probability | Output +1 probability | Error | Cause |
|---|---:|---:|---:|---|
| Qiskit 2.5.1 | 0% | 100% | **100 pp** | Floating-point addition |
| PennyLane 0.45.1 | 0.500375% | 100% | **99.499625 pp** | Floating-point addition |
| Cirq 1.7.0 | 0% | 100% | **100 pp** | Floating-point addition |
| BQSKit 1.2.1 | 0% | 100% | **100 pp** | Repeated small approximations |
| QuiZX 0.3.0 | ≈0% | 100% | **≈100 pp** | Integer overflow |
| tket 2.18.5 | 0% | 100% | **100 pp** | Rotation removed during optimization |

## Input and output for each tool

The input gates below execute from left to right.
`RZ(a)` is a Z rotation by angle `a`, in radians.

All six output circuits act as **identity**. They leave the prepared state unchanged.
Gate comparisons below ignore global phase, which does not affect measurement probabilities.

### Qiskit: 100 pp

```text
RZ(1e17); RZ(pi); RZ(-1e17)
```

**Expected:** RZ(pi), equivalent to Z.

**Output:** identity after `CommutativeCancellation(basis_gates=["rz"])`.

**Cause:** floating-point addition loses the middle rotation.
Both large angles are exact 64-bit floating-point values.
QASM import preserves the middle rotation.
The pass rounds `1e17 + pi` back to `1e17`.
It then adds `-1e17` and gets zero.

**Measurement:** the +1 probability changes from **0% to 100%**.

Standard `transpile` also produces this error at optimization levels 2 and 3.
These checks use `basis_gates=["rz", "sx", "x", "cx"]`.
Level 1 preserves the Z rotation.

### PennyLane: 99.499625 pp

```text
RZ(1e16); RZ(1); RZ(1); RZ(1); RZ(-1e16)
```

**Expected:** RZ(3), because the large positive and negative angles cancel.

**Output:** identity after `merge_rotations(atol=0)`.

**Cause:** floating-point addition loses each of the three one-radian rotations.
The large positive and negative angles then cancel.

**Measurement:** the +1 probability changes from **0.500375% to 100%**.

Qiskit parses the QASM input.
The script converts the parsed gates to PennyLane RZ gates.
PennyLane then runs the optimization.

### Cirq: 100 pp

```text
RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)
```

**Expected:** RZ(pi), equivalent to Z.

**Output:** identity after `eject_z(atol=0)`.

**Cause:** the pass combines large rotation exponents with floating-point addition.
It loses the middle rotation. QASM import preserves that rotation.

**Measurement:** the +1 probability changes from **0% to 100%**.

### BQSKit: 100 pp

```text
[RZ(pi/12000)] × 12,000
```

Each rotation is in its own one-gate block.
The rotations add up to pi.

**Expected:** RZ(pi), equivalent to Z.

**Output:** identity after `ForEachBlockPass(ScanningGateRemovalPass())`, then `UnfoldPass`.

**Cause:** the selected passes permit small errors in each block.
`ScanningGateRemovalPass` uses its default threshold of `1e-8`.
It removes each small rotation separately.
Across 12,000 blocks, these errors add up to a large error.
This example uses an optimization workflow that permits approximation.

**Measurement:** the +1 probability changes from **0% to 100%**.

### QuiZX: approximately 100 pp

```text
RZ(2147483647*pi/4294967295); RZ(2147483648*pi/4294967297)
```

**Expected:** RZ(pi × (1 − 2³²/(2⁶⁴ − 1))), almost exactly Z.
Both input angles are close to pi/2.

**Output:** identity after `full_simp`.

**Cause:** integer overflow during phase addition.
QASM import preserves both rational phases exactly.
The script checks them before optimization.
QuiZX adds the phases with signed 64-bit rational arithmetic.
The two denominators are coprime. Their product is 2⁶⁴ − 1.
This value exceeds the signed 64-bit integer limit.
The addition overflows and returns zero, so the rotation disappears.

**Measurement:** the +1 probability changes from **≈1.33758 × 10⁻¹⁷% to 100%**.
The error is approximately **100 pp**.

**Other tools, same input:** Qiskit, PennyLane, Cirq, BQSKit, and tket preserve this rotation.
Each error is below **2 × 10⁻²³ pp**.
These checks use the same QASM and the passes in this script.
Qiskit's `Optimize1qGatesDecomposition` and tket's `SquashTK1` also preserve the rotation.

PyZX 0.9.0 rounds both imported phases to 1/2.
Its `full_reduce` pass produces Z.
This import approximation causes an error of **≈1.34 × 10⁻¹⁷ pp**.
Only QuiZX produces identity for this input in these checks.

### tket: 100 pp

```text
RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)
```

**Expected:** RZ(pi), equivalent to Z.

**Output:** identity after `RemoveRedundancies`.

**Cause:** the optimization pass removes the middle rotation.
QASM import preserves that rotation.

**Measurement:** the +1 probability changes from **0% to 100%**.

## Run the examples

[uv](https://docs.astral.sh/uv/) uses Python 3.12 and installs the packages listed in `bugs.py`.
Later runs reuse the installed packages.
Every example starts with an OpenQASM 2.0 circuit.

```sh
# All six examples
uv run bugs.py

# Select tools
uv run bugs.py --tools pennylane cirq quizx tket
```

The BQSKit example takes the most time because it optimizes 12,000 separate blocks.
The other examples use at most five gates.

## How the expected probabilities are checked

The reference calculation uses exact fractions for numeric angles and mathematical π.
It derives the expected rotation from the input gates.
The script calculates expected probabilities at 100 decimal digits.
It calculates output probabilities from each tool's output matrix.
It rounds the percentages for display.
