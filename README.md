# Quantum optimizer bugs

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![uv](https://img.shields.io/badge/uv-run-DE5FE9?logo=uv&logoColor=2C1A47)](https://docs.astral.sh/uv/)

> [!WARNING]
> **All popular quantum computing SDKs and optimization tools we have tested do not correctly handle angles, resulting in unsound optimization of quantum circuits.**

Across several quantum SDKs and compilation tools, we find bugs that change expected measurement outcomes.
This repository provides seven OpenQASM examples, and one additional MQT Core 4 example, that reproduce these errors.
Each example shows an almost completely flipped measurement outcome.

Run all examples with the pinned dependencies:

```sh
uv run bugs.py
```

## How bad are the outputs?

For each example, prepare every qubit in `|+>`.
Apply the circuit. Then measure `q[0]` in the X basis.
Compare the expected probability of the **+1** outcome with the output probability.

The error is the absolute probability difference, in **percentage points (pp)**:

| Tool / tested version | Expected +1 probability | Output +1 probability | Error | Cause |
|---|---:|---:|---:|---|
| Qiskit 2.5.1 | 0% | 100% | **100 pp** | Floating-point addition |
| PennyLane 0.45.1 | 0% | 100% | **100 pp** | Floating-point addition |
| Cirq 1.7.0 | 0% | 100% | **100 pp** | Floating-point addition |
| QuiZX 0.3.0 | ≈0% | 100% | **≈100 pp** | Integer overflow |
| tket 2.18.5 | 0% | 100% | **100 pp** | Floating-point addition |
| BQSKit 1.2.1 | 0% | 99.998558% | **99.998558 pp** | Angle applied to the wrong qubit |
| MQT QMAP 3.10.0 / Core 3.10.1 | 0% | 100% | **100 pp** | Angle silently dropped |

## Input and output for each tool

The input gates below execute from left to right.
`rz(a)` is a Z rotation by angle `a`, in radians.

The optimized circuits act as **identity** in five examples.
BQSKit and MQT QMAP produce different incorrect circuits.
None of the input circuits is equivalent to identity.
Gate comparisons below ignore global phase, which does not affect measurement probabilities.

### Qiskit: 100 pp

```text
rz(1e17); rz(pi); rz(-1e17)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after `CommutativeCancellation(basis_gates=["rz"])`.

**Cause:** floating-point addition loses the middle rotation.
Both large angles are exact 64-bit floating-point values.
QASM import preserves the middle rotation.
The pass rounds `1e17 + pi` back to `1e17`.
It then adds `-1e17` and gets zero.

**Measurement:** the +1 probability changes from **0% to 100%**.

### PennyLane: 100 pp

```text
rz(1e17); rz(pi); rz(-1e17)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after `merge_rotations(atol=0)`.

**Cause:** floating-point addition rounds `1e17 + pi` back to `1e17`.
Adding `-1e17` then gives zero, so the middle rotation disappears.

**Measurement:** the +1 probability changes from **0% to 100%**.

Qiskit parses the QASM input.
The script converts the parsed gates to PennyLane RZ gates.
PennyLane then runs the optimization.

### Cirq: 100 pp

```text
rz(1e16*pi); rz(pi); rz(-1e16*pi)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after `eject_z(atol=0)`.

**Cause:** the pass combines large rotation exponents with floating-point addition.
It loses the middle rotation. QASM import preserves that rotation.

**Measurement:** the +1 probability changes from **0% to 100%**.

### QuiZX: approximately 100 pp

```text
rz(2147483647*pi/4294967295); rz(2147483648*pi/4294967297)
```

**Expected:** rz(pi × (1 − 2³²/(2⁶⁴ − 1))), almost exactly Z.
Both input angles are close to pi/2.

**Output:** identity after `full_simp`.

**Cause:** integer overflow during phase addition.
QASM import preserves both rational phases exactly.
The script checks them before optimization.
QuiZX adds the phases with signed 64-bit rational arithmetic.
The two denominators are coprime. Their product is 2⁶⁴ − 1.
This value exceeds the signed 64-bit integer limit.
The addition overflows and returns zero, so the rotation disappears.

**Measurement:** the +1 probability changes from **almost 0% to 100%**.
The expected probability is approximately `1.33758e-19`.
The error is approximately **100 pp**.

### tket: 100 pp

```text
rz(1e16*pi); rz(pi); rz(-1e16*pi)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after `RemoveRedundancies`.

**Cause:** floating-point addition loses the middle rotation.
tket stores these angles in half-turns: `1e16`, `1`, and `-1e16`.
Adding `1e16 + 1` rounds back to `1e16`. Adding `-1e16` then gives zero.
QASM import preserves that rotation.

**Measurement:** the +1 probability changes from **0% to 100%**.

### BQSKit: 99.998558 pp

```qasm
gate g(a) x,y { rz(a) x; id y; }
qreg q[2];
g(0) q[1],q[0];
g(pi) q[0],q[1];
```

The first block is identity. The second applies Z to `q[0]`.

**Expected:** Z on `q[0]` and identity on `q[1]`.

**Output:** an incorrect two-qubit circuit after `BlockConversionPass("variable")`
and `ExtractDiagonalPass`.

**Cause:** the extraction pass moves a diagonal phase angle from the second block into
the first without remapping their different qubit orderings, so the angle lands on the wrong qubit.
The script uses a local cost threshold of `1e-12`, seed `75`, and one synthesis start.
The large error appears when the phase is merged into the other block.
Using the same qubit ordering in both blocks preserves the result.
The incorrect merge is in [ExtractDiagonalPass](https://github.com/BQSKit/bqskit/blob/1.2.1/bqskit/passes/processing/extract_diagonal.py).

**Measurement:** the +1 probability changes from **0% to 99.998558%**.

### MQT QMAP: 100 pp

```text
crz(2*pi) q[0],q[1]
```

**Expected:** Z on the control `q[0]` and identity on the target `q[1]`.
A full 2*pi RZ rotation is −identity. Controlling it makes that minus sign a Z gate on the control.

**Output:** `cx q[0],q[1]` after native `mqt.qmap.sc.map_` on two connected qubits.
CNOT leaves `|++>` unchanged.

**Cause:** the native mapper drops the rotation angle. It assumes that every controlled gate is a CNOT.
It accepts the unsupported CRZ gate and silently emits a CNOT in the release build.
The input retains the correct CRZ semantics in MQT Core before mapping.
This is an input-validation failure in the native mapper; it also occurs with its optional optimizations disabled.
The native mapper expects two-qubit gates to be lowered to CNOTs first.
The bug is that it returns a wrong circuit instead of rejecting the unsupported gate.
The Qiskit `compile_` wrapper in this version also accepts it and emits a CNOT.
This example does not test MQT Core 4's compiler pipeline.
See the related [upstream limitation](https://github.com/munich-quantum-toolkit/qmap/issues/843)
and [mapper source](https://github.com/munich-quantum-toolkit/qmap/blob/v3.10.0/src/sc/Mapper.cpp).

**Measurement:** the +1 probability changes from **0% to 100%**.

#### Additional example: MQT Core 4.0.0, approximately 100 pp

```qasm
ctrl @ pow(80599999999716496) @ gphase(1) q[0];
```

This OpenQASM 3 gate is `p(K)` on `q[0]`, with `K = 80599999999716496`.
It applies the phase `e^(iK)` when `q[0]` is 1.

**Expected:** `p(K)`. `K` mod 2π is approximately −9×10⁻⁵, so the gate is almost identity.

**Output:** `p(-3.1413486180563552)`, almost Z, after `compile_program` with `OutputFormat.QCO_OPTIMIZED`.

**Cause:** the global-phase normalization pass multiplies the phase by `K`.
It then reduces the result modulo the 64-bit floating-point value of 2π.
That value is approximately 2.4×10⁻¹⁶ less than the real 2π.
Over approximately 1.3×10¹⁶ periods, this difference adds up to approximately π.
QASM import preserves `K` exactly.
See [`normalizeAngle`](https://github.com/munich-quantum-toolkit/core/blob/v4.0.0/mlir/lib/Dialect/MQT/Utils/Angles.cpp)
and [the normalization pass](https://github.com/munich-quantum-toolkit/core/blob/v4.0.0/mlir/lib/Dialect/MQT/Transforms/NormalizeGlobalPhases.cpp).

**Measurement:** the +1 probability changes from **almost 100% to almost 0%**.
MQT QMAP 3.10.0 requires MQT Core 3.10, so `bugs.py` runs this example in a separate environment.

## How the expected probabilities are checked

The script calculates each expected probability from the QASM input with exact arithmetic and 100-digit precision.
It does not use the tool under test for this calculation, and it rejects unsupported gates.
It then compares this value with the probability from the tool's output circuit.

Released under the [MIT license](LICENSE).
