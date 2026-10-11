# quantum-optimizer-bugs

Six OpenQASM circuits whose optimization erases a substantial rotation
and produces the identity circuit. The resulting measurement probabilities change
by **42 to 100 percentage points**. One standalone script reproduces all six:

```sh
uv run bugs.py
```

## How bad are the outputs?

For every example, prepare **|+X> = (|0> + |1>)/√2**, apply the circuit, then
measure the indicated Pauli observable. The table compares the probability of
outcome **+1** for the ideal input circuit and the actual output circuit.

**Percentage points (pp)** are the absolute difference between the probabilities:
0% → 100% is a 100 pp change. These are concrete measurement witnesses, selected
from X and Y measurements, rather than a search over all states and measurements.

| Tool / tested version | Measure | Ideal probability | Output probability | Change | Error occurs during |
|---|---|---:|---:|---:|---|
| Qiskit 2.5.1 | Y | 96.601954% | 50% | **46.601954 pp** | Optimization |
| PennyLane 0.45.1 | X | 0.500375% | 100% | **99.499625 pp** | Optimization |
| Cirq 1.7.0 | X | 0% | 100% | **100 pp** | Optimization |
| BQSKit 1.2.1 | Y | 92.073549% | 50% | **42.073549 pp** | Local approximation |
| QuiZX 0.3.0 | X | ≈0% | 100% | **≈100 pp** | Optimization |
| tket 2.18.5 | X | 0% | 100% | **100 pp** | Optimization |

## Input and output for each tool

Circuits below execute left to right. Angles are in radians; `pi` means
mathematical π. `P(a) = diag(1, exp(i*a))` is a phase gate; `RZ(a)` is a Z rotation.
Expected unitaries are stated up to global phase. All six actual outputs are
the **identity**: the circuit's action has been erased.

### Qiskit: 46.601954 pp

```text
[P(d); P(d); X; P(-d); P(-d); X] × 100,000
d = 3e-6
```

**Expected:** RZ(≈1.2). Each six-gate block implements RZ(4d), so the small
rotations accumulate into a substantial rotation across 600,000 gates.
**Actual:** identity after `CommutativeCancellation` + `InverseCancellation`.
Local commutation cutoffs allow the whole accumulated rotation to disappear.

Measuring Y changes P(+1) from **96.601954% to 50%**.

### PennyLane: 99.499625 pp

```text
RZ(1e16); RZ(1); RZ(1); RZ(1); RZ(-1e16)
```

**Expected:** RZ(3), because the large positive and negative angles cancel.
**Actual:** identity after `merge_rotations(atol=0)`. Floating-point addition
loses the three one-radian contributions while merging the angles.

Measuring X changes P(+1) from **0.500375% to 100%**.
The QASM is parsed with Qiskit and adapted to native PennyLane RZ operations
before running the PennyLane optimization.

### Cirq: 100 pp

```text
RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)
```

**Expected:** RZ(pi), equivalent to Z. **Actual:** identity after `eject_z(atol=0)`.
The importer preserves the rotation; adding the large exponents during
optimization loses the middle Z.

Measuring X changes P(+1) from **0% to 100%**: the two circuits give opposite
outcomes with certainty.

### BQSKit: 42.073549 pp

```text
[RZ(1e-4)] × 10,000
Each rotation is in its own one-gate block.
```

**Expected:** RZ(≈1). **Actual:** identity after
`ForEachBlockPass(ScanningGateRemovalPass())` followed
by `UnfoldPass`. Each block's small rotation is removed independently.

Measuring Y changes P(+1) from **92.073549% to 50%**.
This is an explicitly configured approximation workflow: accepted local errors
accumulate into a large circuit error.

### QuiZX: approximately 100 pp

```text
RZ(2147483647*pi/4294967295); RZ(2147483648*pi/4294967297)
```

**Expected:** RZ(π × (1 − 2³²/(2⁶⁴ − 1))), a rotation almost exactly equal to Z.
Both input angles are close to π/2. **Actual:** identity after `full_simp`.

QASM import preserves both rational phases exactly; the script checks them
before simplifying. Phase fusion adds them using 64-bit rational arithmetic.
Their coprime denominators have least common multiple **2⁶⁴ − 1**, which exceeds
the signed 64-bit range. Overflow causes the sum to become zero and the rotation
to disappear. This is an **integer overflow during optimization**.

Measuring X changes P(+1) from **≈1.33758 × 10⁻¹⁷% to 100%**: approximately
**100 pp**.

### tket: 100 pp

```text
RZ(1e16*pi); RZ(pi); RZ(-1e16*pi)
```

**Expected:** RZ(pi), equivalent to Z. **Actual:** identity after
`RemoveRedundancies`. QASM import retains the middle rotation; the optimization
pass erases it.

Measuring X changes P(+1) from **0% to 100%**: the two circuits give opposite
outcomes with certainty.

## Run the examples

[uv](https://docs.astral.sh/uv/) installs the dependencies declared in `bugs.py`
and uses Python 3.12. The first run installs the quantum frameworks; later runs
reuse its cache. Every tool receives a circuit parsed from OpenQASM 2.0.

```sh
# All six examples
uv run bugs.py

# Select tools
uv run bugs.py --tools pennylane cirq quizx tket
```

Qiskit's 600,000 gates and BQSKit's 10,000 separate block optimizations dominate
runtime. The large circuits demonstrate accumulated error.

## How the expected probabilities are checked

The expected rotations are derived analytically, using exact rational arithmetic
for numeric angles and mathematical π. QASM literals for `3e-6` and `1e-4` spell
out their exact stored binary64 values. The reference also accounts for the X
gates in Qiskit's echo; it does not use a framework's floating-point matrix as
the ideal answer.

Ideal probabilities are evaluated at 100 decimal digits. Output probabilities
are calculated from each tool's native output matrix. The script prints only
the input, expected rotation, output, measurement probabilities, and their
absolute difference in percentage points. The displayed percentages are rounded.

The tests check the ideal probabilities at 100 and 160 digits for agreement
within absolute 10⁻⁷⁰ and compare the Qiskit echo reference with an independent
matrix calculation.

Run the small reference tests with:

```sh
uv run --with mpmath==1.3.0 python -m unittest discover -s tests -v
```
