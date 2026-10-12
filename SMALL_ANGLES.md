# Small-angle bugs

The examples in [README.md](README.md) mostly use very large angles.
The examples on this page use only small angles.
Tolerance checks, numerical shortcuts, and other angle-handling errors remove or change small rotations.
Many small errors then add up to an almost completely flipped measurement outcome.

Run all small-angle examples with the pinned dependencies:

```sh
uv run slowbugs.py
```

Run one or more examples by tool name, for example `uv run slowbugs.py tket qiskit`.
The Cirq and Qrisp examples take several minutes each.

## How bad are the outputs?

The protocol is the same as in [README.md](README.md#how-bad-are-the-bugs).
Prepare every qubit in `|+>`. Apply the circuit. Then measure `q[0]` in the X basis.
Compare the expected probability of the **+1** outcome with the output probability.

| Tool / tested version | Largest input angle | Input gates | Expected +1 probability | Output +1 probability | Error | Cause |
|---|---:|---:|---:|---:|---:|---|
| Qiskit 2.5.1 | 1.4×10⁻⁴ | 44,432 | 0% | 100% | **100 pp** | Two-qubit synthesis drops small entangling components |
| PennyLane 0.45.1 | 9.0×10⁻³ | 353 | 0% | 100% | **100 pp** | Tolerance is 400 times larger than intended |
| Cirq 1.7.0 | 1×10⁻⁵ | 300,003 | 0.500375% | 100% | **99.499625 pp** | Commutation check with a relative tolerance |
| PyZX 0.9.0 | 0.01 | 163 | 0.010812% | 99.999984% | **99.989172 pp** | Phase teleportation loses the sign of a fused phase |
| QuiZX 0.3.0 | 0.01 | 317 | 0.000486% | 99.402667% | **99.402182 pp** | Integer overflow during phase addition |
| tket 2.18.5 | 6.1×10⁻¹⁰ | 47 | ≈0% | 99.836003% | **99.836003 pp** | Shortcut divides two values near zero |
| BQSKit 1.2.1 | 8.7×10⁻² | 144 | 0% | 100% | **100 pp** | `compile` uses a documented distance as a cost threshold |
| Qrisp 0.9.9 | 7×10⁻⁵ | 44,880 | ≈0% | 100% | **≈100 pp** | Two small rotations are cancelled as inverses |
| qBraid 0.12.1 | 7.7×10⁻³ | 500 | 0% | ≈100% | **≈100 pp** | `a/b*pi` changes to `a/(b*pi)` |

All input angles are in radians.
The script checks that each tool's import preserves every input angle before optimization.
[slowbugs.py](slowbugs.py) calculates each expected probability with a state-vector simulation at 100-digit precision.
This reference does not use the tool under test, and it rejects unsupported gates.
The script also does the two checks in [README.md](README.md#how-the-results-are-checked).
It calculates the expected probability again with input angles calculated in 64-bit floating-point arithmetic.
It also simulates the tool's output gates with the same reference.

## Input and output for each tool

The input gates below execute from left to right.
Gate comparisons ignore global phase, which does not affect measurement probabilities.

### Qiskit: 100 pp

```text
h q[1]; h q[2]
crz(2*pi/44430) q[0],q[1]; crz(2*pi/44430) q[0],q[2]; ...   (44,430 crz gates)
```

Each `crz` angle is approximately 1.4×10⁻⁴.
The `h` gates put both targets in `|0>` for the `|+>` input.
Thus each `crz` applies a small phase to `q[0]`. The 44,430 phases add up to π.

**Expected:** Z on `q[0]`.

**Output:** two single-qubit gates on `q[1]` and `q[2]`, and no gate on `q[0]`,
after `transpile(circuit, basis_gates=["cx", "u"])` with default settings.

**Cause:** two-qubit synthesis drops the entangling component of each `crz`.
The default transpiler resynthesizes each two-qubit block with exact synthesis.
But the Weyl decomposition still uses a fidelity of 1 − 10⁻⁹ to classify the block.
The entangling component of each `crz` is within this fidelity of identity, so the synthesis emits no CNOT.
The local factors survive: for this input angle, an isolated `crz(a)` becomes `rz(a/2)` on its target.
This target rotation contributes only a global phase on the prepared `|0>` state, so the phase on `q[0]` is lost.
The alternating targets prevent Qiskit from merging the `crz` gates first.
See [the default fidelity](https://github.com/Qiskit/qiskit/blob/2.5.1/crates/synthesis/src/two_qubit_decompose/common.rs#L25),
[its use in exact synthesis](https://github.com/Qiskit/qiskit/blob/2.5.1/crates/synthesis/src/two_qubit_decompose/basis_decomposer.rs#L651-L657),
[the identity check](https://github.com/Qiskit/qiskit/blob/2.5.1/crates/synthesis/src/two_qubit_decompose/weyl_decomposition.rs#L646-L668),
and [the retained local factors](https://github.com/Qiskit/qiskit/blob/2.5.1/crates/synthesis/src/two_qubit_decompose/weyl_decomposition.rs#L712-L721).

**Measurement:** the +1 probability changes from **0% to 100%**.

### PennyLane: 100 pp

```text
h q[1]; x q[1]
cu1(pi/351) q[0],q[1]   (351 times)
```

The target `q[1]` is in `|1>`, so the 351 phases on `q[0]` add up to π.

**Expected:** Z on `q[0]`.

**Output:** 2,108 gates and no CNOT after `qml.transforms.unitary_to_rot`.
All rotation angles in the output are 0.
`qml.compile` and `qml.transforms.decompose` also lose the phases when an explicit gate set forces decomposition of `QubitUnitary`:

```python
gate_set = {"RY", "RZ", "CNOT", "Hadamard", "PauliX"}
qml.compile(tape, basis_set=gate_set)
qml.transforms.decompose(tape, gate_set=gate_set)
```

Their default settings preserve `QubitUnitary` and the correct measurement probability.

**Cause:** a tolerance check is 400 times larger than intended.
The decomposition counts the CNOTs that a two-qubit unitary needs.
It selects zero CNOTs if the trace is close to 4, with `allclose(trace, 4, atol=1e-7)`.
The comment says that a tolerance of approximately 10⁻⁷ is intended.
But the default relative tolerance of 10⁻⁵ also applies, so the real tolerance is approximately 4×10⁻⁵.
Each `cu1(pi/351)` passes this check. The zero-CNOT decomposition then loses its phase.
See [the CNOT count](https://github.com/PennyLaneAI/pennylane/blob/v0.45.1/pennylane/ops/op_math/decompositions/unitary_decompositions.py#L518-L521).

**Measurement:** the +1 probability changes from **0% to 100%**.

Qiskit parses the QASM input. The script converts `cu1` to a PennyLane `QubitUnitary`.
PennyLane's own QASM importer (pennylane-qiskit) does the same conversion.
That package requires an older Qiskit, so the script does not use it directly.

### Cirq: 99.499625 pp

```text
h q[1]; cx q[0],q[1]
rz(1e-5) q[1]   (300,000 times)
cx q[0],q[1]
```

For the `|+>` input, the circuit applies `rz(3)` to `q[0]`.

**Expected:** `rz(3)` on `q[0]`. The +1 probability is 0.500375%.

**Output:** all 300,003 gates remain after `cirq.transformers.insertion_sort_transformer`.
The two `cx` gates become adjacent and cancel algebraically, leaving the rotations on `q[1]`.

**Cause:** the commutation check uses a relative tolerance.
The transformer moves the last `cx` to the left past each gate that commutes with it.
Cirq compares the two products with `np.allclose`, which has a default relative tolerance of 10⁻⁵.
Thus `rz(1e-5)` on the target counts as commuting with `cx`, although it does not.
Each move loses 10⁻⁵ of phase on `q[0]`. The 300,000 moves lose 3 radians.
See [the transformer](https://github.com/quantumlib/Cirq/blob/v1.7.0/cirq-core/cirq/transformers/insertion_sort.py#L78)
and [the commutation check](https://github.com/quantumlib/Cirq/blob/v1.7.0/cirq-core/cirq/ops/raw_types.py#L1144).

**Measurement:** the +1 probability changes from **0.500375% to 100%**.

`optimize_for_target_gateset` uses the same transformer with `reorder_operations=True`.

### PyZX: 99.989172 pp

```text
s; x; rz(0.01) (155 times); x; rz(0.01); s; x; sdg; rz(0.01)
```

All gates act on `q[0]`. The two `x` gates around the 155 rotations reverse their direction.

**Expected:** the +1 probability is 0.010812%.

**Output:** seven gates after `zx.teleport_reduce`, which is PyZX's phase-teleportation pass.
The output phase is approximately +1.57. The correct phase at this location is approximately −1.55.
`full_reduce` followed by `extract_circuit` gives the correct result.

**Cause:** phase teleportation loses the sign of a fused phase.
The pass fuses the 155 rotations into one phase of 1.55. The `x` gates make this phase negative.
The pass records the sign as a multiplier.
When the fused phase meets a Clifford phase, the pass sets the multiplier to +1.
The next fusion then adds the phase instead of subtracting it.
See [the phase fusion](https://github.com/zxcalc/pyzx/blob/v0.9.0/pyzx/simplify.py#L276-L290).

**Measurement:** the +1 probability changes from **0.010812% to 99.999984%**.

### QuiZX: 99.402182 pp

```text
rz(0.01) (314 times); rz(0.001); rz(0.002); rz(0.003)
```

The angles add up to 3.146, so the circuit is almost Z.

**Expected:** `rz(3.146)`. The +1 probability is 0.000486%.

**Output:** one rotation of approximately 0.155 after `full_simp`.

**Cause:** integer overflow during phase addition.
This is the same cause as the [QuiZX example](README.md#quizx-approximately-100-pp), but with ordinary decimal angles.
QuiZX stores each phase as a fraction of π with 64-bit integers.
Import converts each decimal angle to a fraction with a large denominator, for example 383/120323 for 0.01.
The conversion error is less than 10⁻⁷ of the angle.
Four different decimal angles are sufficient to make the sum of the fractions overflow.
See [the import conversion](https://github.com/zxcalc/quizx/blob/quizx-py-v0.3.0/quizx/src/circuit.rs#L457-L466)
and [the phase addition](https://github.com/zxcalc/quizx/blob/quizx-py-v0.3.0/quizx/src/phase.rs#L164-L166).

**Measurement:** the +1 probability changes from **0.000486% to 99.402667%**.

### tket: 99.836003 pp

```text
s q[0]
(rz(6.134e-10) q[0]; rx(1.064e-10) q[0]; ry(6.317e-11) q[0]; z q[0];
 rx(8.644e-11) q[0]; y q[0]; cx q[0],q[1]; x q[0]; cx q[0],q[1])   (5 times)
s q[0]
```

Each block acts almost as identity on `q[0]` and applies X to `q[1]`.
The `cx` gates keep the single-qubit gates of adjacent blocks apart.

**Expected:** Z on `q[0]`.

**Output:** two `TK1` gates after `FullPeepholeOptimise()`.
`SynthesiseTket` and `PeepholeOptimise2Q` give the same result.

**Cause:** a shortcut divides two values that are almost zero.
tket squashes each chain of single-qubit gates into one rotation.
A shortcut applies when a test value is within 10⁻¹¹ of zero.
The shortcut calculates an angle from the ratio of two small components, each slightly larger than 10⁻¹¹.
This angle is approximately 0.61 radians, but it should be almost zero.
Five blocks add up to approximately π.
See [the shortcut](https://github.com/CQCL/tket/blob/v2.18.5/tket/src/Gate/Rotation.cpp#L146-L160)
and [the tolerance](https://github.com/CQCL/tket/blob/v2.18.5/tket/include/tket/Utils/Constants.hpp#L38).

**Measurement:** the +1 probability changes from **almost 0% to 99.836003%**.

### BQSKit: 100 pp

```text
(rz(pi/36) q[0]; cx q[0],q[1]; cx q[0],q[2]; cx q[0],q[3])   (36 times)
```

For the `|+>` input, the `cx` gates have no effect. The rotations add up to exactly π.

**Expected:** `rz(pi)` on `q[0]`, equivalent to Z.

**Output:** the 108 `cx` gates and no rotations,
after `bqskit.compile` with `optimization_level=2` and `synthesis_epsilon=1e-3`.

**Cause:** `compile` uses a distance limit as a cost limit.
Its documentation says that `synthesis_epsilon` is the maximum distance between the target and the circuit unitary.
At optimization level 2, `compile` divides the circuit into 3-qubit blocks.
It gives `synthesis_epsilon` to `ScanningGateRemovalPass` as `success_threshold`.
This pass compares `success_threshold` with a Hilbert–Schmidt cost, not with a distance.
For a small rotation by θ, this cost is approximately θ²/8. For `rz(pi/36)`, it is 9.5×10⁻⁴, less than 10⁻³.
Thus the pass removes the rotation from each block.
BQSKit's own distance for this change is 0.044, approximately 44 times the documented maximum.
The 36 blocks together lose π.
With the default `synthesis_epsilon=1e-8`, the pass removes rotations up to 2.8×10⁻⁴.
Then the distance is approximately 14,000 times the documented maximum.
See [the `synthesis_epsilon` documentation](https://github.com/BQSKit/bqskit/blob/1.2.1/bqskit/compiler/compile.py#L259-L261),
[the gate-deletion workflow](https://github.com/BQSKit/bqskit/blob/1.2.1/bqskit/compiler/compile.py#L1344-L1368),
[the threshold documentation](https://github.com/BQSKit/bqskit/blob/1.2.1/bqskit/passes/processing/scan.py#L41-L44),
and [the removal check](https://github.com/BQSKit/bqskit/blob/1.2.1/bqskit/passes/processing/scan.py#L131).

**Measurement:** the +1 probability changes from **0% to 100%**.

### Qrisp: approximately 100 pp

```text
u3(0.00007,0,0) q[0]   (44,880 times)
```

Each gate is a Y rotation by 7×10⁻⁵. The rotations add up to 3.1416.

**Expected:** `ry(3.1416)`, almost Y. The +1 probability is approximately 1.3×10⁻¹¹.

**Output:** identity after the default compiler, `QuantumSession.compile()`.

**Cause:** the compiler cancels two small rotations as if they were inverses.
The gate-fusion pass compares the inverse of one gate with the next gate.
It removes both gates if their matrices differ by 10⁻⁴ or less.
`ry(-a)` and `ry(a)` differ by approximately 1.4 × `a`, so each pair of `u3(7e-5,0,0)` gates passes this check.
Each removed pair loses 1.4×10⁻⁴. All 22,440 pairs lose approximately π.
See [the inverse check](https://github.com/eclipse-qrisp/Qrisp/blob/c30b323ad93507bd58b25e75a705e946a29d790b/src/qrisp/circuit/pass_management/passes/fuse_adjacents.py#L434).
Qrisp has no 0.9.9 tag. This link uses the commit that sets the version to 0.9.9.

**Measurement:** the +1 probability changes from **almost 0% to 100%**.

The script prepares `|+>` with an `h` gate in the Qrisp session.
Thus the compiler cannot use the fact that a new qubit starts in `|0>`.
The script removes this H from the output before the comparison.

### qBraid: approximately 100 pp

```text
rz(1/408*pi) q[0]   (454 times)
rz(-pi/408) q[0]    (46 times)
```

All angles are ±π/408, approximately ±0.0077. They add up to π.

**Expected:** rz(pi), equivalent to Z.

**Output:** 454 rotations of `1/(408*pi)` and 46 correct rotations of `-pi/408`.
They add up to approximately 0, so the output is almost identity.
qBraid's IonQ device prepares the job with `transform` and then `prepare` before it submits it.
The script runs these two steps offline.

**Cause:** the angle rewriting ignores operator precedence.
Before conversion to IonQ JSON, qBraid replaces `pi` in each angle with a decimal value.
A regular expression takes only the operand next to `pi`, so `1/408*pi` becomes `1/(408*pi)`.
Each affected angle becomes π² times too small.
This bug does not depend on the size of the angle. For example, `rz(2/3*pi)` becomes `rz(0.2122)`.
The same angle written as `pi/408` stays correct.
See [the regular expression](https://github.com/qBraid/qBraid/blob/v0.12.1/qbraid/passes/qasm/compat.py#L229)
and [its evaluation](https://github.com/qBraid/qBraid/blob/v0.12.1/qbraid/passes/qasm/compat.py#L253-L255).

**Measurement:** the +1 probability changes from **0% to almost 100%**.
