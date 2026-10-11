# quantum-optimizer-bugs

One standalone [uv](https://docs.astral.sh/uv/) script reproducing six large,
finite discrepancies in quantum circuit importers and optimizers. One example
per tool; no parameter sweeps, overflow example, or dependency on another project.

```sh
uv run bugs.py
```

The script declares its pinned dependencies inline and uses Python 3.12. The first
run installs the quantum frameworks. Subsequent runs reuse uv's cache. All tools
receive circuits parsed from OpenQASM 2.0.

Select examples, export their exact inputs, or save full precision results:

```sh
uv run bugs.py --tools pennylane cirq quizx tket
uv run bugs.py --qasm-dir qasm-inputs --json results.json
```

The Qiskit example expands 600,000 gates. BQSKit optimizes 10,000 individual
blocks with one worker. These two examples dominate runtime; they retain the
large accumulated errors. The implementation uses a specialized rational phase
reference instead of general circuit simulation and has no report UI dependency.

## Examples

Circuits execute left to right. `P(a) = diag(1, exp(i*a))` and
`RZ(a) = diag(exp(-i*a/2), exp(i*a/2))`.

| Tool / version | Input | Expected, up to global phase | Observed output | Unitary error E | Measurement change |
|---|---|---|---|---:|---:|
| Qiskit 2.5.1 | `[P(d); P(d); X; P(-d); P(-d); X]` × 100,000; d=3e-6 | RZ(≈1.2) | Identity | 0.591040413323 | 46.6019543 pp |
| PennyLane 0.45.1 | RZ(1e16); RZ(1); RZ(1); RZ(1); RZ(-1e16) | RZ(3) | Identity | 1.36327752005 | 99.4996248 pp |
| Cirq 1.7.0 | RZ(1e16*pi); RZ(pi); RZ(-1e16*pi) | Z | Identity | 1.41421356237 | 100 pp |
| BQSKit 1.2.1 | RZ(1e-4) × 10,000, one gate per block | RZ(≈1) | Identity | 0.494807918509 | 42.0735492 pp |
| QuiZX 0.3.0 | RZ(10000000000000114 radians) | RZ(≈−3.13309558725) | Identity | 1.41120620719 | 99.9981950 pp |
| tket 2.18.5 | RZ(1e16*pi); RZ(pi); RZ(-1e16*pi) | Z | Identity | 1.41421356237 | 100 pp |

QASM numeric literals for `d` and `1e-4` spell out their exact stored binary64
values. `pi` denotes mathematical π in the reference.

The report checks import and optimization separately. QuiZX's discrepancy occurs
entirely during QASM import: `full_simp` adds zero error. The other five cases
have zero measured import error and a discrepancy introduced by optimization.
Stage errors are separate comparisons; they are not summed to obtain the total.

BQSKit uses an explicitly configured local approximation workflow:
`ForEachBlockPass(ScanningGateRemovalPass(), calculate_error_bound=True)` followed
by `UnfoldPass`. It records a bound but does not enforce it as a global error
limit. This demonstrates accumulated local approximation, rather than a failure
of its default compilation pipeline or a promise of exact semantics.

PennyLane's pinned QASM plugin is incompatible with the pinned Qiskit version.
Its example therefore uses Qiskit's QASM parser and an explicit adapter mapping
the parsed RZ instructions to native PennyLane RZ operations. The report names
this adapter; this is a test of PennyLane optimization, not its own QASM plugin.

## What does the unitary error mean?

The error is the **operator norm distance up to global phase**:

```text
E(U, V) = min over phi in [0, 2*pi) of ||U - exp(i*phi) V||_op
```

`phi` is one real number. `exp(i*phi)` multiplies the entire matrix V by a single
complex scalar of magnitude one. The minimization removes the unobservable
global phase: U and −U represent the same channel, so their error is zero.
It does not allow a separate phase for each basis state or each matrix entry.
The operator/spectral norm is the largest singular value, equivalently the
maximum Euclidean output-vector difference over normalized input states for
that fixed phase. E is unitless; it is not a percentage or an angle.

This is a standard metric on unitary operators modulo global phase. See
Haah, Kothari, O'Donnell, and Tang, *Query-optimal estimation of unitary channels
in diamond distance*, [Definition 1.4, Equation (3), page 4](https://www.cs.cmu.edu/~./odonnell/papers/estimating-unitaries.pdf#page=4)
([arXiv:2302.14066](https://arxiv.org/abs/2302.14066)). Their Proposition 1.6
relates it to diamond-norm distance; it is not itself the diamond norm.

All six circuits have the form `global * diag(1, exp(i*theta))` after their X
gates cancel. If Δ is the difference in relative phases, reduced to [−π, π],
the exact minimization has the simple formula:

```text
E = 2 * sin(abs(Δ) / 4)
```

The optimal global phase centers the two relative eigenphases at ±Δ/2. For
example, identity versus Z gives E=√2. The code computes this formula directly,
without a heuristic phase alignment or a numerical search over phi.

The report also gives a concrete measurement witness: prepare |+X>, measure X
or Y, and compare the probability of outcome +1. Its percentage-point change
has direct physical meaning. The displayed witness is the stronger of these
two probes, not a maximum over every possible state and measurement.

## Independent reference

Native gate angles are accumulated with exact `Fraction` arithmetic, preserving
binary64 values and rational multiples of mathematical π. The analytic ideal is
derived independently from each QASM fixture. The reference handles the X gates
with `X P(a) X = global * P(-a)`; it sums angles before evaluating any trigonometry.
It does not use a framework's numerical matrix as the ideal answer.

Final errors and probabilities are evaluated at 100 and 160 decimal digits. The
script checks agreement within absolute 10⁻⁷⁰. This is an empirical numerical
stability check, not an interval proof or a formal equivalence certificate.
JSON retains 80-digit results; the report rounds them for readability. Native
execution failures and reference failures produce a nonzero process exit status.

Run the small reference tests without installing every quantum framework:

```sh
uv run --with mpmath==1.3.0 python -m unittest discover -s tests -v
```

Tests cover exact residuals, global phase and periodicity, tiny nonzero errors,
stage attribution, and phase accumulation against independent dense matrix
products. To exercise the framework loaders and passes, run `uv run bugs.py`.
