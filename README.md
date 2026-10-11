# Quantum optimizer bugs

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![uv](https://img.shields.io/badge/uv-run-DE5FE9?logo=uv&logoColor=2C1A47)](https://docs.astral.sh/uv/)

> [!WARNING]
> **Optimization and mapping change a certain or almost-certain −1 outcome into a certain or almost-certain +1 outcome in these examples.**

Across several quantum SDKs and compilation tools, we find bugs that change expected measurement outcomes.
This repository provides seven OpenQASM examples that reproduce these errors.
Each example shows an almost completely flipped measurement outcome.

Run all seven examples with the locked dependencies:

```sh
uv run --locked bugs.py
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
| BQSKit 1.2.1 | 0% | 99.998558% | **99.998558 pp** | Incorrect qubit ordering |
| QuiZX 0.3.0 | ≈0% | 100% | **≈100 pp** | Integer overflow |
| tket 2.18.5 | 0% | 100% | **100 pp** | Floating-point addition |
| MQT QMAP 3.10.0 / Core 3.10.1 | 0% | 100% | **100 pp** | Unsupported gate silently replaced |

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

**Cause:** the extraction pass moves a diagonal phase from the second block into
the first without remapping their different qubit orderings.
The script uses a local cost threshold of `1e-12`, seed `75`, and one synthesis start.
The large error appears when the phase is merged into the other block.
Using the same qubit ordering in both blocks preserves the result.
The incorrect merge is in [ExtractDiagonalPass](https://github.com/BQSKit/bqskit/blob/1.2.1/bqskit/passes/processing/extract_diagonal.py).

**Measurement:** the +1 probability changes from **0% to 99.998558%**.

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

### MQT QMAP: 100 pp

```text
crz(2*pi) q[0],q[1]
```

**Expected:** Z on the control `q[0]` and identity on the target `q[1]`.
A full 2*pi RZ rotation is −identity. Controlling it makes that minus sign a Z gate on the control.

**Output:** `cx q[0],q[1]` after native `mqt.qmap.sc.map_` on two connected qubits.
CNOT leaves `|++>` unchanged.

**Cause:** the native mapper assumes that every controlled gate is a CNOT.
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

## How the expected probabilities are checked

The reference reads each QASM input independently of the tool being tested.
It expands the BQSKit block and computes the diagonal phases with exact rational numbers and mathematical π.
It supports only the gates in these fixtures and rejects other gates.
It then calculates the expected X measurement probability at 100 decimal digits.
It calculates output probabilities from each tool's output matrix.
It accounts for the tools' different qubit orderings and rounds percentages for display.

The script prints whether each error exceeds 99 pp and exits with a failure if any example does not reproduce.
`bugs.py.lock` pins both direct and indirect dependencies.

The large-angle examples use finite, legal QASM angles to expose lost rotations.
An optimizer should preserve their measurement probabilities or reject inputs it cannot handle reliably.

Released under the [MIT license](LICENSE).
