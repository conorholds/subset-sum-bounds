# Tunable checksum moment multipliers

The constants `8.5` and `8.3` in the checksum variance proof can be replaced
by independently chosen rational multipliers. Jointly optimizing the filters,
multipliers, and parameter frontier gives an interval-checked conditional
exponent **`alpha = 0.48972642`**, improving
[PR 13](https://github.com/conorholds/subset-sum-bounds/pull/13)'s `0.48973396`.
Keeping its degree-fourteen filters unchanged already gives `0.48973227`;
changing the filters supplies the larger additional gain.

There is also a bound on parameter headroom: with the predecessor's fixed
filters and at least the same positive reserves, the generalized numerical
model cannot achieve `alpha <= 0.489732159`. This bound permits arbitrary
positive integer block counts and does not assume a balanced parameter search.
It does not bound other filters, smaller reserves, or different algorithms.

The original checker and its receipts are unchanged. The explicit new model
has its own [interval checker](interval_check.py), exact
[candidates](candidates.json), [certificate](CERTIFICATE.json), and
[frontier bound](BOUND.json). These certify numerical conditions and the
localized argument below; the complete parameter-to-algorithm implication
still requires independent mathematical review.

## Predecessor and source

[predecessor.json](predecessor.json) retains the exact `degree14_retuned`
witness from `conorholds/subset-sum-bounds@7d8fa41`, its path, and its source
file hash. The first three candidates use those same coefficients and `kappa`;
the joint candidates change them explicitly.
The balanced cost slice and reserves come from PRs
[12](https://github.com/conorholds/subset-sum-bounds/pull/12) and
[13](https://github.com/conorholds/subset-sum-bounds/pull/13).

The manuscript remains pinned to
[`openai/math@1581d4c`](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026).
Its algorithm, filtering architecture, and original proof are OpenAI's.
The precise source obligations are `eq:checksum-domination-one`,
`eq:checksum-domination-two`, `eq:checksum-unconditional-domination`, and
`eq:checksum-error-exponent` in `checksum.tex`, together with the uniform
finite-circle convergence argument in `filtering.tex`.
[SOURCE.json](../filter-budgets/SOURCE.json) retains their hashes.

## What changes in the proof

Write `Delta = eta-tau`, `Gamma = 2tau-eta`, and `a = Delta/Gamma`.
For a rounded frequency pattern write `U = A+E-D` and `W = 2D-A-E`.
The coordinate-average scores satisfy

```text
|U-ubar| <= 3 epsilon,   |W-vbar| <= 4 epsilon.
```

Choose fixed rational multipliers `lambda1, lambda2` with

```text
lambda1 >= 1,   lambda2 >= 1,
3 lambda1 < 30,   (4a+3) lambda2 < 30.                         (1)
```

The strict upper inequalities keep the source's existing rounding allowance.
They are sufficient conditions for this model, not claimed necessary for every
possible refinement of the estimator.

In this argument `s` denotes the pattern's coordinate-test probability
(the manuscript's local `alpha`); `alpha` here remains the outer exponent.

In the sampling branch with `W <= Gamma`, `s = 1` and `U < Delta`.
Since `Delta-U > 0` and `lambda1 >= 1`,

```text
-U = -Delta + (Delta-U)
   <= -Delta + lambda1 (Delta-U)
   <= -Delta + lambda1 (Delta-ubar) + 3 lambda1 epsilon
   <  -Delta + lambda1 (Delta-ubar) + 30 epsilon.             (2)
```

In the other sampling branch, `W > Gamma`, `s = Gamma/W <= 1`,
and `aW-U > 0`. The identity `s aW = Delta` gives

```text
-s U = -Delta + s (aW-U)
         <= -Delta + lambda2 (aW-U)
         <= -Delta + lambda2 (a vbar-ubar)
                    + lambda2 (4a+3) epsilon
         <  -Delta + lambda2 (a vbar-ubar) + 30 epsilon.     (3)
```

Replace the source's two coordinate exponents by

```text
R1(x) = lambda1 (Delta-u(x))
R2(x) = lambda2 (a v(x)-u(x)).
```

On its good event, the source bounds the match count with `18 epsilon`
of slack. Equations (2) and (3) therefore give exactly its existing
`48 epsilon` pointwise domination. Only after removing the pattern/good-event
indicators from the nonnegative upper bound do we average over the original
independent phase law. The resulting mean squared error exponent is

```text
-tau-Delta+chi+48 epsilon = -eta+chi+48 epsilon < 0,
```

provided both generalized moments are below `2^chi`. No independent-phase
claim is made under the conditioned good-event law. The subexponential
pattern count and polynomial target list are absorbed by this strict negative
linear exponent as before.

The tables, list construction, direct/sampling branch rule, inverse
multiplicities, caps, and arithmetic procedures are unchanged. The multipliers
are analysis parameters. Fixed finite cosine polynomials, positive `kappa`,
and fixed rational multipliers retain the boundedness and uniform convergence
used to pass from continuous to finite-circle moments. Their strict continuous
slack supplies a sufficiently large fixed cutoff; this package does not compute
a useful finite-instance crossover.

The tests retain exact rational examples showing why a multiplier below one
cannot establish these same linear dominations, including the second branch
with `lambda2 < s < 1`. They also test both worst-sign score-rounding
errors and rejection when either allowance in (1) is exhausted.

The resource bounds also retain their general-parameter form. Entropy ordering
gives `D-A,D-E>=0`. With `s=min(1,Gamma/W)` and high-digit density
`H_pi=tau-s(D-A)`, one has `Delta<=H_pi<=tau` and low-entry density
`eta-H_pi+s(D-E)=Delta+sW<=tau`. Hence both digit ranges and both expected
entry lists retain the source's bounds. Its good-event allowances remain
`14 epsilon` for lists and `18 epsilon` for matches. The cap at
`2^ceil((tau+20 epsilon)n)` still bounds every execution; direct patterns fit
because `sU>=Delta`, and sampled patterns still use `2^ceil(tau n)` draws.
Multiplicities and prefix counts keep `O(n)` bits. None of these steps uses
the two multipliers. This imports the source's deterministic score routines,
bounded random draws, and complex arithmetic implementation; it is not a new
implementation of those routines.

## Search and numerical evidence

For fixed scores and a fixed parameter tuple, each log moment is a convex
log-sum-exp in its multiplier. Its derivative is the weighted mean of
`Delta-u` or `a v-u`; its second derivative is `ln(2)` times the weighted
variance. Bisection of the monotone derivative handles interior optima and
both endpoint optima without a random seed.

[multiplier_search.py](multiplier_search.py) uses a 256-point exploratory
grid, at most 60 derivative bisections per multiplier and 45 outer bisections
on `eta in [.5179,.5181]`. It uses the same 17-block cost balance and reserves
as the predecessor. The search retains the full rounded parameter tuple,
exact nine-decimal multipliers, limits, evaluation count, runtime versions,
and input/code hashes in [discovery.json](discovery.json). The floating-point
outer bisection proposes a witness; acceptance comes from the interval replay.

| Candidate | Outer `alpha` | Multipliers | Purpose |
| --- | ---: | --- | --- |
| `predecessor_degree14` | 0.48973396 | 8.5, 8.3 | Reproduce the predecessor through the new path. |
| `multipliers_only` | 0.48973396 | approximately 8.55760, 8.36198 | Isolate the gain without changing any budget. |
| `multipliers_retuned` | 0.48973227 | approximately 8.55737, 8.36171 | Restore every cost and round conservatively. |
| `joint_degree14_retuned` | 0.48972644 | approximately 8.77011, 8.56892 | Optimize the actual feasible parameter frontier jointly. |
| `joint_degree20_retuned` | 0.48972642 | approximately 8.77034, 8.56913 | Add six cosine terms per filter as a bounded extension. |

The exact candidates and generated certificate are authoritative for all
digits and interval endpoints. The improved tuple preserves `epsilon=1e-8`,
`guard=1e-5`, a `0.0001` moment reserve, a `0.00002` deficiency reserve,
`0.000026` coverage reserve, `0.00003` representation reserve, `0.000016`
filtered-list reserve, `0.000010` visited-pair reserve, and `0.00001` assembly
reserve. Conservative rounding adds a small further margin.

The interval checker adapts the frozen checker's direct inner quadratures,
zero-free split-contour argument, and analytic outer remainder bounds. Both
integrands **and both remainder bounds** use the supplied multipliers. It
reuses the original 14 budget checks and adds the two allowances in (1).
The old fixed-multiplier checker rejects the improved tuple; its old
certificate cannot be used for this changed model.

An exploratory extra tilt using the branch inequalities `W<=Gamma` and
`W>Gamma` gave no improvement at the retained design. This is a local search
observation, not an impossibility result for more general moment envelopes.

### Joint objective and controls

The joint search targets the actual parameter boundary. Let `b_i(x)` be the
`eta` at which moment `i` reaches `eta-.0001` on the balanced cost slice,
for filter/multiplier vector `x`. The feasible endpoint is `min(b_1,b_2)`.
Implicit differentiation gives

```text
gradient b_i = - gradient_x log2(I_i) / (derivative_eta log2(I_i)-1).
```

[joint_search.py](joint_search.py) differentiates normalizers, correlations,
scores, multipliers, and this boundary. It optimizes a smoothed maximum of the
two negative endpoints, with temperatures `.01,.001,.0001,.00001` after scaling
the endpoint objective by 1000. This gives the two constraints their actual
parameter sensitivities, rather than minimizing moments at fixed budgets.
Each stage allows 100 iterations, at most 24 backtracking trials per step,
and the whole run allows 4,000 frontier evaluations including the rounded
retained point. Roots use 40 bisections in `[.5175,.51828]`, on 128 nodes.
Coefficients stay within `.35` of their starting values, `kappa` in
`[.95,1.12]`, and multipliers in `[7.5,9.5]` and `[7.5,9]`.

[joint-degree14.json](joint-degree14.json) retains the first run and
[joint-degree20.json](joint-degree20.json) its zero-padded higher-degree
continuation. Their 633 and 882 evaluations include exact-decimal retention.
The degree-twenty extension improves the outer exponent by only `0.00000002`
in this bounded run. This is evidence of limited local headroom for that
extension, not an optimum over cosine polynomials.

A second degree-twenty run starts from the manuscript's original degree-six
filters rather than the improved predecessor. Its
[input](source-start.json) and [725-evaluation result](joint-source-start.json)
reach the same endpoint to about `1e-11` in `eta`. This is an exploratory
initialization control; only the five candidates in `candidates.json` are
claimed to pass the retained interval certificate.

[odd_probe.py](odd_probe.py) tests another local escape route: hold the cosine
coefficients, `kappa`, and multipliers fixed, and allow the first twenty sine
terms in each filter. A 256-node, `.001`-step Hessian probe uses 1,606 frontier
evaluations. The least unfavorable direction of the optimized convex blend
of its two numerical Hessians worsens the endpoint at amplitudes
`.001,.003,.01,.03,.1`; [odd-control.json](odd-control.json) retains that
direction, inputs, and all five outcomes. These are finite exploratory
controls, not an interval proof of local optimality or a bound on asymmetric
filters. The tests replay the outcomes and check reflection symmetry.

The tests compare every implicit derivative at the start and both retained
designs against separate perturbations on 128- and 256-node grids. They also
check zero-padding invariance, evaluation caps, input hashes, and the separate
continuous interval path. No floating-point discovery output is accepted as
the numerical certificate.

## A certified floor for fixed filters and reserves

Here the filters are fixed at the predecessor and the reserve assumptions
mean **at least** the values just listed. In particular,

```text
d-(.5-tau) >= .00002,   J rho delta-d >= .000026,
.5-delta-gamma >= .00003,   eta-chi >= .0001,
guard >= .00001,   alpha-tau >= .00001,
```

together with the two cost reserves in (4). Let
`H=h(1/4)`, `C=.500046`, `A=H-1+.00003 < 0`, and
`B=2H-1.5+.00003 > 0`. Substitution of the deficiency, coverage, and
representation requirements into the two cost inequalities gives, for
arbitrary positive integer `J` and any block density `rho`,

```text
(1+1/J) tau >= .5 + A rho + C/J + .000016
(1+1/J) tau >=  1 + B rho + C/J - eta + .000020.              (4)
```

Weight these inequalities by `B/(B-A)` and `-A/(B-A)` to eliminate `rho`.
For `J=17`, call the resulting affine floor `T(eta)`. Its slope is strictly
negative. This argument is necessary for every density, including densities
that violate other constraints, so relaxing those constraints cannot invalidate
the lower bound.

More explicitly, write the weighted constant without `C/J` as `K(eta)`.
The floor for general `J` is `(J K+C)/(J+1)`. If `K>=C`, every floor is
at least `C>.5`. Otherwise it decreases with `J`, so all `J<=17` have
`tau>=T(eta)`. For `J>=18`, the elementary entropy/coverage bound gives

```text
tau > .5 - (1-H)/(J+1) >= .5 - (1-H)/19 > .49006.
```

Thus any tuple competitive with the retained result must satisfy the
17-block affine floor, regardless of its actual block count.

Set `eta0=.517987190`. At `tau=T(eta0)`, consider

```text
F(lambda) = log2 integral M^2 2^{lambda(eta0-T(eta0)-u)}.
```

For every allowed `1<=lambda<10`, convexity gives
`F(lambda)>=F(lambda0)+(lambda-lambda0)F'(lambda0)`. At 256 nodes and
70 interval digits, [frontier_bound.py](frontier_bound.py) encloses `F` and
`F'`, including analytic errors for the derivative integrand. The tangent
lower bound is strictly above `eta0-.0001` over the entire interval `[1,10]`.
[BOUND.json](BOUND.json) retains the positive infeasibility margin.

Consequently `alpha<T(eta0)+.00001` is impossible in this model. Indeed,
`tau<T(eta0)` and `eta<eta0` contradict the decreasing affine cost floor.
If `eta>=eta0`, the first-moment residual at any fixed `lambda>=1` increases
by

```text
(lambda-1)(eta-eta0) + lambda(T(eta0)-tau) > 0,
```

so the first-moment condition fails. This excludes the entire parameter
region, rather than a finite search grid. The conclusion uses the first
moment only; adding the second moment and other constraints cannot weaken it.
It explains why further parameter tuning with these filters and reserves has
little headroom.

## Replay

Use the existing [requirements](../filter-budgets/requirements.txt). From the
repository root:

```sh
python3 -B research/chernoff-multipliers/test_package.py
python3 -B research/chernoff-multipliers/interval_check.py --check-receipt research/chernoff-multipliers/CERTIFICATE.json
python3 -B research/chernoff-multipliers/frontier_bound.py --output .context/chernoff-bound.json
```

To repeat discovery separately, pass `--output` and `--candidates-output` to
`multiplier_search.py`, writing to temporary paths. The verifier consumes
`candidates.json` beside its code; replace inputs deliberately and regenerate
the receipt before claiming a different candidate passed.

The next question is whether a structural change to the densities or a
sharper pattern-dependent variance envelope can improve the new frontier.
No full algorithm implementation, practical speedup, or new selected project
result is claimed here.
