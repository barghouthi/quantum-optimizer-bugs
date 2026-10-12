# Quantum optimizer bugs

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![uv](https://img.shields.io/badge/uv-run-DE5FE9?logo=uv&logoColor=2C1A47)](https://docs.astral.sh/uv/)

> [!CAUTION]
> **All popular quantum computing SDKs and optimization tools we have tested do not correctly handle angles, resulting in unsound optimization of quantum circuits.**

Across many quantum SDKs and compilation tools, we find many bugs that change expected measurement outcomes.
This repository provides OpenQASM examples that reproduce these errors.
Each example shows an almost completely flipped measurement outcome!

This page focuses on bugs that are fast to reproduce.
[SMALL_ANGLES.md](SMALL_ANGLES.md) covers the slower examples, where many small rotations add up to an almost completely flipped outcome.

Run the fast examples with the pinned dependencies:

```sh
uv run bugs.py
```

Run the slower small-angle examples with the pinned dependencies:

```sh
uv run slowbugs.py
```

The [individual_scripts](individual_scripts/) folder has one short script for each tool.
Each script reads the OpenQASM inputs of that tool's examples, runs the tool, and prints the output circuit.
Then it checks whether the output is equivalent to the input. For example:

```sh
uv run individual_scripts/qrisp_bugs.py
```

Most scripts print OpenQASM 2. Some tools use other formats.
The qBraid script prints IonQ JSON, the CUDA-Q script prints compiled QIR, and the quilc script prints Quil.
The MQT Core 4 script prints MLIR and OpenQASM 2. The BQSKit script prints raw unitary blocks as matrices.
The quilc script needs Docker.

The equivalence check is in [equivalence.py](individual_scripts/equivalence.py). It does not use any tool under test.
It reads each angle exactly and calculates both unitaries with 100-digit arithmetic.
Then it compares the unitaries up to global phase.
Run `uv run individual_scripts/equivalence.py` to test the check on circuit pairs with known answers.

## How bad are the bugs?

For each example, prepare every qubit in `|+>`, apply the circuit, and measure `q[0]` in the X basis.
Compare the ideal probability of the **+1** measurement outcome for the original input circuit with the probability for the optimized output circuit.
In every example below, optimization almost completely flips this probability: from approximately 0% to 100%, or vice versa.

The error is the absolute probability difference, in **percentage points (pp)**:

| Tool / tested version | Ideal input +1 probability | Optimized output +1 probability | Error | Cause |
|---|---:|---:|---:|---|
| Qiskit 2.5.1 | 0% | 100% | **100 pp** | Floating-point addition |
| PennyLane 0.45.1 | 0% | 100% | **100 pp** | Floating-point addition |
| Cirq 1.7.0 | 0% | 100% | **100 pp** | Floating-point addition |
| QuiZX 0.3.0 | ≈0% | 100% | **≈100 pp** | Integer overflow |
| tket 2.18.5 | 0% | 100% | **100 pp** | Floating-point addition |
| BQSKit 1.2.1 | 0% | 99.998558% | **99.998558 pp** | Angle applied to the wrong qubit |
| Qrisp 0.9.9 | 0% | 100% | **100 pp** | Floating-point addition |
| qBraid 0.12.1 / PyQASM 1.0.4 | 100% | 0% | **100 pp** | Floating-point addition in gate decomposition |
| CUDA-Q 0.16.0 | 0% | 100% | **100 pp** | Floating-point addition |
| quilc 1.26.0 (pyQuil's compiler) | 0% | 100% | **100 pp** | Floating-point addition |
| MQT QMAP 3.10.0 / Core 3.10.1 | 0% | 100% | **100 pp** | Angle silently dropped |

## Input and output circuits for each tool

The input gates below execute from left to right.
`rz(a)` is a Z rotation by angle `a`, in radians.

The optimized circuits act as **identity** in eight examples.
BQSKit, qBraid, and MQT QMAP produce different incorrect circuits.
The qBraid input is equivalent to identity; the other inputs are not.
Gate comparisons below ignore global phase, which does not affect measurement probabilities.

### Qiskit: 100 pp

```text
rz(1e17); rz(pi); rz(-1e17)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after `transpile()` with default settings.

**Cause:** floating-point addition loses the middle rotation.
Both large angles are exact 64-bit floating-point values.
QASM import preserves the middle rotation.
The default optimization level 2 runs `CommutativeCancellation`.
This pass rounds `1e17 + pi` back to `1e17`.
It then adds `-1e17` and gets zero.
Optimization level 1 does not run this pass and gives the correct result.
See [the angle sum](https://github.com/Qiskit/qiskit/blob/2.5.1/crates/transpiler/src/passes/commutation_cancellation.rs#L236-L276).

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
rz(pi); rz(-1e16*pi)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after `RemoveRedundancies`.

**Cause:** floating-point addition loses the first rotation.
tket stores these angles in half-turns: `1` and `-1e16`.
The pass adds the two angles. `1 + (-1e16)` rounds to `-1e16`.
A rotation by `-1e16` half-turns is identity, so the pass removes it.
QASM import preserves both angles.
`rz(-1e16*pi)` alone is identity, and the pass removes it correctly.
See [the angle sum](https://github.com/CQCL/tket/blob/v2.18.5/tket/src/Transformations/RedundancyRemoval.cpp#L168-L171).

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

### Qrisp: 100 pp

```text
rz(1e17); rz(pi); rz(-1e17)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after the default compiler, `QuantumSession.compile()`.

**Cause:** floating-point addition loses the middle rotation.
QASM import preserves all three angles.
The gate-fusion pass rounds `1e17 + pi` back to `1e17`.
Adding `-1e17` then gives zero. The pass then removes the gate.
See [`fuse_adjacents`](https://github.com/eclipse-qrisp/Qrisp/blob/c30b323ad93507bd58b25e75a705e946a29d790b/src/qrisp/circuit/pass_management/passes/fuse_adjacents.py#L223).
Qrisp has no 0.9.9 tag. This link uses the commit that sets the version to 0.9.9.

**Measurement:** the +1 probability changes from **0% to 100%**.

The script prepares `|+>` with an `h` gate in the Qrisp session.
Thus the compiler cannot use the fact that a new qubit starts in `|0>`.
The script removes this H from the output before the comparison.

### qBraid: 100 pp

```text
u3(0,1e17,0); rz(-1e17)
```

**Expected:** identity up to global phase. With a zero first parameter,
`u3(0,a,0)` is `diag(1,exp(i*a))`, equivalent to `rz(a)` up to global phase.
The following `rz(-a)` cancels it.

**Output:** six RX/RZ gates equivalent to Z after qBraid's
`qasm3_to_ionq(qasm2_to_qasm3(source))`.
This example tests the SDK's conversion to IonQ format.
The script uses the QASM 3 path to lower U3 to IonQ's supported gates.
It calls the converters directly to avoid discovering unrelated optional backends.

**Cause:** qBraid unrolls the U3 gate through PyQASM before emitting IonQ rotations.
PyQASM's decomposition includes a final `rz(phi + pi)`.
With `phi = 1e17`, floating-point addition rounds `phi + pi` back to `phi`.
Losing that pi rotation changes the circuit from identity to Z.
The input's large angles are exact 64-bit floating-point values.
The script checks that an independent QASM import preserves the angles and the correct outcome.
See [the U3 decomposition](https://github.com/qBraid/pyqasm/blob/v1.0.4/src/pyqasm/maps/gates.py)
and [qBraid's IonQ conversion](https://github.com/qBraid/qBraid/blob/v0.12.1/qbraid/transpiler/conversions/qasm3/qasm3_to_ionq.py).

**Measurement:** the +1 probability changes from **100% to 0%**.

The script simulates the emitted IonQ RX/RZ rotations with Qiskit,
without running a Qiskit optimization pass.

### CUDA-Q: 100 pp

```text
rz(1e17); rz(pi); rz(-1e17)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity after the default compilation pipeline. `cudaq.get_state` runs this pipeline.

**Cause:** floating-point addition loses the middle rotation.
The kernel keeps all three exact angles before optimization.
The `quake-simplify` pass merges the rotations with a 64-bit floating-point addition.
`1e17 + pi` rounds back to `1e17`. Adding `-1e17` then gives zero.
The pass then removes the zero rotation.
See [the rotation merge](https://github.com/NVIDIA/cuda-quantum/blob/0.16.0/cudaq/lib/Optimizer/Transforms/QuakeSimplify.cpp#L363-L379)
and [the default pipeline](https://github.com/NVIDIA/cuda-quantum/blob/0.16.0/cudaq/lib/Optimizer/CodeGen/Pipelines.cpp#L196).

**Measurement:** the +1 probability changes from **0% to 100%**.

CUDA-Q has no QASM importer.
Qiskit parses the QASM input. The script builds the CUDA-Q kernel from the parsed gates.
The kernel first applies H to prepare `|+>`.
If the pipeline removes every gate, the CUDA-Q unitary has no qubit.
Thus the script compares the output state.

### quilc: 100 pp

```text
rz(1e17); rz(pi); rz(-1e17)
```

**Expected:** rz(pi), equivalent to Z.

**Output:** identity (`HALT` only) after quilc. pyQuil uses this compiler.

**Cause:** floating-point addition loses the middle rotation.
quilc's rewriting rule merges adjacent RZ gates and adds their angles in 64-bit floating point.
`1e17 + pi` rounds back to `1e17`. Adding `-1e17` then gives zero.
Another rule then removes the zero rotation.
See [the rewriting rules](https://github.com/quil-lang/quilc/blob/7aef642/src/compressor/rewriting-rules.lisp#L85-L93).

**Measurement:** the +1 probability changes from **0% to 100%**.

The script runs quilc 1.26.0 from the official Docker image. It pins the image by digest.
Without Docker, the script skips this example.
Qiskit parses the QASM input. The script sends the parsed gates to quilc as Quil.
We also built quilc 1.27.0, the latest release, from source. It gives the same result.

quilc also removes a single large rotation. `RZ(1.3e18)` compiles to identity.
But 1.3×10¹⁸ mod 2π is approximately 3.139, so the gate is almost Z.
quilc divides the angle by 2π in floating point and checks whether the result is a whole number.
Every floating-point value above 2⁵² is a whole number.
Thus every angle above approximately 2.8×10¹⁶ passes this check.
See [`full-rotation-p`](https://github.com/quil-lang/quilc/blob/7aef642/src/compressor/rewriting-rules.lisp#L54-L60).

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

## How the results are checked

For each example, the script prints the output gates and does three checks:

1. **Expected probability.** The script calculates the expected probability from the QASM input.
   It uses symbolic phase arithmetic and 100-digit numerical evaluation.
   It does not use the tool under test for this calculation, and it rejects unsupported gates.
2. **Floating-point input angles.** Most angles, such as `pi` or `0.01`, do not have an exact 64-bit floating-point value.
   The script calculates each input angle in 64-bit floating-point arithmetic, as Qiskit's QASM parser does.
   Then it calculates the expected probability again with these angles.
   The result changes by less than 10⁻¹⁴ pp. Thus, input rounding does not cause the errors.
   Some tools keep angles in other forms. For example, tket uses half-turns, and PyZX and QuiZX use fractions of π.
   Before optimization, the script checks the imported angles, exactly or within a small bound.
   For Qiskit, PennyLane, Cirq, tket, and MQT QMAP, it also simulates the imported circuit.
   The result must agree with the expected probability.
   This check uses the angles that the tool keeps, so it replaces the floating-point check for these tools.
   For example, tket keeps `-1e16*pi` exactly as `-1e16` half-turns, but the nearest 64-bit float is off by 0.38 radians.
3. **Output gates.** The script reads the gates of the tool's output circuit. It simulates these gates with the same 100-digit reference.
   This check does not use the tool's own matrix code.
   For CUDA-Q, the script reads the gates from the compiled QIR.

The script stops with an error if the tool's own result and the output-gate simulation do not agree.
It also stops if input rounding changes the expected probability by more than 10⁻¹⁴ pp, or if a matrix or state is not valid.
An example reproduces only if the error is more than 99 pp for the tool's own result and for both checks.

Run the regression tests for these checks with the same pinned dependencies:

```sh
uv run --locked bugs.py --test-harness
```

## What's the solution?

The best solution is to rely on symbolic, arbitrary precision arithmetic instead of floating point. For example, I recently did that for [tzap](https://github.com/qqq-wisc/tzap). Better yet, I believe that quantum software should be formally verified.

Released under the [MIT license](LICENSE).
