# Fractional coverage without an integer packing loss

The alias-coverage argument can use a fractional sumset inequality instead of
packing an integer number `J` of small restrictions into a half. This replaces
the sufficient asymptotic condition `J rho delta > d`, `J rho < 1`, with
**`delta > d`**, while preserving the alias extractor and its resource bounds.

With the preceding degree-twenty filters and checksum multipliers **unchanged**,
the refined conditional model admits **`alpha = 0.48972043`**, compared with
the tunable-multiplier proposal's `0.48972642`. The useful-block separation
reserve is slightly larger than before. Retuning the multipliers gives the
same exponent at the retained precision.

For these fixed filters and the stated minimum reserves, a separate bound
rules out `alpha <= 0.489720315962704626`. It ranges over all `eta,tau,rho`
and admissible multipliers, not just the searched balanced slice. It does not
bound other filters or smaller reserves.

These are a proposed proof refinement and interval-checked **conditional**
budgets. The selected result is unchanged; independent mathematical review of
the application and the complete algorithm remains pending.

## Source and predecessor

The manuscript is pinned to
[`openai/math@1581d4c`](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026),
with all seven section hashes in [SOURCE.json](../filter-budgets/SOURCE.json).
The changed step is the 17-way packing argument immediately before
`eq:alias-useful-block` in `aliases.tex`. The preliminary search in `setup.tex`,
alias list/pair accounting, and checksum construction are preserved.

The sumset inequality is classical: see Madiman, Marcus, and Tetali,
*Entropy and set cardinality inequalities for partition-determined functions*,
[Corollary 3.5, arXiv:0901.0055v3](https://arxiv.org/pdf/0901.0055v3).
The contribution here is its application to this coverage interface and the
resulting checked budget, not discovery of that inequality.

[predecessor.json](predecessor.json) pins canonical project main `f25a303`,
the exact tunable-multiplier candidate at `0eedec7`, its complete source file
hash, and the mathematical reference. The predecessor retains the filters
and analytic search approach developed through PRs
[12](https://github.com/conorholds/subset-sum-bounds/pull/12) and
[13](https://github.com/conorholds/subset-sum-bounds/pull/13).
The deficiency-interface audit in
[PR 14](https://github.com/conorholds/subset-sum-bounds/pull/14) supplies the
complementary conditioning and variance context; the proof below is explicit
about its own bounds.

## 1. A fractional restriction inequality

For arbitrary fixed integer weights write

```text
Sigma(U) = {sum_{i in X} c_i : X subset U},
D(U) = |U| - log2 |Sigma(U)|.
```

For a fixed `h`-set `V` and `1<=q<=h`,

```text
average_{U subset V, |U|=q} D(U) <= (q/h) D(V).               (1)
```

Here is a self-contained specialization of the sumset projection proof.
Choose the least binary mask representing each sum on `V`, under one fixed
coordinate order, and call this set of representatives `R`. A restriction of
a least representative to `U` is also least for its sum on `U`: otherwise a
smaller replacement inside `U` would preserve the full sum and improve the
full representative. Thus the sum map is injective on the projection `R_U`,
and `|R_U|<=|Sigma(U)|`.

Let `Z` be uniform on `R`. The entropy chain rule and the fact that conditioning
reduces entropy give

```text
sum_{|U|=q} H(Z_U)
  >= binom(h-1,q-1) H(Z)
   = binom(h-1,q-1) log2 |Sigma(V)|.
```

Indeed, expand each `H(Z_U)` in the fixed coordinate order; conditioning on
the earlier coordinates in `U` gives at least the entropy conditioned on all
earlier coordinates. Every coordinate appears `binom(h-1,q-1)` times.
Combining this with `H(Z_U)<=log2|Sigma(U)|`, then subtracting from `q`,
proves (1). The representative law need not be a product law. Signed, zero,
and repeated weights cause no exception. The case `q=0` is immediate.

Averaging (1) over a uniform `h`-subset of `[n]` gives

```text
E D(V_q) <= (q/h) E D(V_h),     0<=q<=h<=n, h>0.              (2)
```

The two-stage restriction is uniform: every `q`-subset has the same number
of `h`-supersets. Equivalently, expected deficiency per selected coordinate
is nondecreasing with the restriction size.

This does **not** follow from superadditivity alone. The set function
`floor(|U|/3)` satisfies nonnegativity, the unit increment bound,
superadditivity on disjoint sets, and the exchange bound, but has
`E f(V_3)=1 > (3/5) E f(V_5)=3/5`. The sumset structure is essential.

The factor in (2) is sharp for arbitrary integer weights. With `z` zero
weights and the remaining weights distinct powers of two,
`D(U)` is exactly the number of zero positions in `U`; hence
`E D(V_q)=qz/n`, giving equality at every pair of sizes.

## 2. Coverage, floors, and probability laws

Put `h=floor(n/2)`, `m=2 floor(rho n/2)`, and `q=m/2`.
The low-deficiency case is the deterministic input condition
`E D(V_h)<=d n`. As in the source, a uniform flip makes a fixed alias a
uniform subset, independent of the random `m`-block. Conditional on selecting
exactly `q` block positions, those positions `S` are a uniform `q`-subset of
`[n]`; the balance event has probability at least `1/(m+1)`.
Equation (2) gives the finite-size statement

```text
E D(S) <= d n m/(2h).                                       (3)
```

If `delta>d`, then for all sufficiently large `n`, Markov's inequality gives

```text
Pr[D(S)<=delta m | balance] >= 1-d n/(2h delta) > 0.          (4)
```

No independent copies of `S` are asserted. No conditioning of the input's
deficiency expectation on preliminary continuation is used. The high/low
case split and the source's smaller-prime independence remain as before.

There is also a `1-O(1/n)` bound. The source's exchange coupling puts each
conditional martingale value in an interval of length at most one; its
conditional variance is therefore at most `1/4`. Orthogonality gives
`Var D(V_q)<=q/4=m/8`. With

```text
a_n = m (delta-d n/(2h)) > 0,
Pr[D(S)<=delta m | balance] >= a_n^2/(m/8+a_n^2).             (5)
```

This is Cantelli's inequality; apply Markov to `(X+v/t)^2` for a centered
variable of variance `v` to obtain its one-sided form. Write `g=delta-d>0`.
For `n>=2`, (3) and the two floors give

```text
a_n >= rho n g - 2g - 2 rho d.
```

Consequently a sufficient cutoff for the bound in (5) to be at least one half
for every larger integer `n` is

```text
n >= max(2, ceil(4/rho+4d/g), ceil(1/(2 rho g^2))).           (6)
```

Above the first floor cutoff, `a_n>=rho n g/2`; also `m/8<=rho n/8`.
This proves (6), including odd `n`. It is a component cutoff, not a practical
solver size or a complete algorithm cutoff.

In fact the upper tail is exponentially small. The same conditional-range
coupling and Hoeffding's conditional moment lemma give

```text
Pr[D(S)>delta m | balance] <= exp(-4a_n^2/m) <= 2^(-4a_n^2/m). (7)
```

For a centered variable supported in an interval of length one, the log
moment generating function has second derivative at most `1/4`, even under
exponential tilting. Integrating twice gives the bound `t^2/8`. Apply it
conditionally to each of the `q` martingale differences, multiply their
conditional bounds, and use Markov with `t=4a_n/q`. This proves (7), with
no independence of the sampled positions required. The classical background
is [Hoeffding's bounded-sum inequalities](https://doi.org/10.1080/01621459.1963.10500830).
At fixed positive `rho,g`, (7) is `exp(-Omega(n))`; beyond the first floor
cutoff it is at most `2^(-rho*g^2*n)`. The audit retains integer-rounded
dyadic bounds rather than expanding enormous rational probabilities.

On the useful event there are at least `2^((1/2-delta)m)` distinct sums on
`S`. The source's sign translation, sparse distinction event, cardinality
layer, and independent smaller-prime argument apply unchanged when
`gamma<1/2-delta`. They give an inverse-polynomial recovery probability per
trial; the existing polynomial repetitions and unconditional cap bounds
still apply. The lemma is an analysis improvement and adds no algorithmic
work. `J` was a proof packing count, distinct from the checksum's unchanged
`Theta(sqrt(n))` pattern blocks.

### Positive transformed controls and the threshold below `d`

The improvement is not limited to zero-weight examples. Take `z` original
weights equal to one and the other weights `(n+1)2^j`. Choose any legal
isolation power `N0`, set all tie weights to one, and put `C=1+n N0` as in
the source. A restriction with `l` small weights and `k` large weights has

```text
|Sigma(U)| = 2^k(l+1),     D(U)=l-log2(l+1).                 (8)
```

Different large subsets differ by at least `C(n+1)` in their main parts;
all small-weight changes together contribute at most `Cn+n`. Since `C>n`,
these cannot cancel. For a fixed large subset, the `l+1` small multiplicities
give distinct sums. This proves (8) for legal positive transformed inputs.

One finite control uses `n=14,z=7,N0=256,m=4,d=3/25,delta=13/100`.
The exact half-support product proves `E D(V_7)<=14d`, and every two-position
restriction is useful. Nevertheless the largest strictly fitting count is
`J=3`, with old coverage margin `-3/350`. This is a coverage-interface
example, not a full feasible algorithm tuple or an input-length-guard claim
at that small size.

There is an asymptotic obstruction below the new threshold. Fix rational
`0<delta<d<1/2` and `0<rho<1/2`, and use the positive family with `z=2dn`
at compatible sizes where this is integral. The low-deficiency condition
holds because `D(U)<=l` and `E l=(z/n)floor(n/2)<=dn` on a half.
For a uniform `q=m/2` restriction, its small count `L` is hypergeometric with
mean `2dq`. If `D(S)<=delta m`, (8) implies

```text
L <= 2delta q+log2(q+1) <= (d+delta)q
```

once `log2(q+1)<=(d-delta)q`. A direct sampling-without-replacement bound gives

```text
Pr[D(S)<=delta m | balance]
  <= exp(-2q(d-delta)^2) <= 2^(-rho(d-delta)^2 n).            (9)
```

For completeness, the moment generating function of a hypergeometric count
is bounded by that of `q` independent Bernoulli draws with the same mean.
It is an elementary symmetric mean of positive numbers; averaging any two
entries at fixed sum increases that mean, whose maximum is attained when
all entries coincide. For a centered Bernoulli, the log moment generating
function has second derivative at most `1/4`, hence is at most `t^2/8`.
Markov's inequality with `t=4(d-delta)` proves the first bound in (9).
The second uses `2q=rho n` at the compatible sizes and `e>2`.

Thus lowering `delta` a fixed amount below `d` cannot give a uniform
inverse-polynomial useful-block guarantee under this low-deficiency
hypothesis. This is not an impossibility result for the complete algorithm:
preliminary search or a different recovery method could handle such inputs.
Original weights have `O(n)` bits, so the source's exponential input-length
guard holds eventually. [AUDIT.json](AUDIT.json) retains exact support
controls, dyadic cutoffs, compatible sizes, and a positive tail exponent
for `delta=d/2` at the improved density.

## 3. Explicit changed model and reserves

The new schema has no `blocks` field. [budget.py](budget.py) replaces the old
coverage and block-space conditions with `rho*(delta-d)>0`; the other budget
inequalities are unchanged. [interval_check.py](interval_check.py) reuses the
unchanged continuous moment bounds from the tunable-multiplier package and
checks both multiplier rounding allowances. Existing checkers, witnesses,
receipts, and the selected result are preserved.

The retained improved candidates use these minimum reserves:

| Quantity | Reserve |
| --- | ---: |
| `d-(1/2-tau)` | `0.000020` |
| Useful-block separation per `n`: `rho*(delta-d)` | `0.0000016` |
| `1/2-delta-gamma` | `0.000030` |
| Filtered list cost below `tau` | `0.000016` |
| Visited pair cost below `tau` | `0.000010` |
| Bit-length guard | `0.000010` |
| Continuous moment ceiling below `eta` | `0.000100` |
| `alpha-tau` | `0.000010` |

The previous packing reserve was `0.000026` in `J rho delta-d`, giving
useful-block separation `0.000026/17` per `n`, about `0.0000015294`.
The new `0.0000016` is larger. No reserve was silently dropped to obtain
the improvement. The `48 epsilon` checksum slack and `20 epsilon` cap
allowance remain, with `epsilon=1e-8`.

Let `H=h(1/4)`. Balancing the list/pair costs with these reserves gives the
searched slice

```text
rho = (eta-0.500004)/(H-1/2),
tau = [1/2+(H-1/2+0.00005)rho+0.0000016+0.000016]/(1+rho).
```

[search.py](search.py) retains the old filters and bisects `eta` at most 45
times in `[.5175,.51828]` on 256 nodes. The fixed-multiplier control uses
94 moment evaluations. The additional multiplier-tuned run uses bounded
convex derivative bisections and does not improve the rounded exponent.
Both use the predecessor's conservative decimal rounding policy.

| Candidate | Conditional `alpha` | Purpose |
| --- | ---: | --- |
| `predecessor_original_budget` | `0.48972642` | Replay the exact prior tuple after removing its proof-only `blocks` field. |
| `fractional_fixed_multipliers` | `0.48972043` | Isolate the coverage refinement; filters and multipliers unchanged. |
| `fractional_retuned` | `0.48972043` | Additional bounded multiplier control. |

[candidates.json](candidates.json), [CERTIFICATE.json](CERTIFICATE.json), and
the discovery receipts retain all exact digits. At the improved tuple the
largest strictly fitting integer `J` is 17, and its old coverage margin is
negative. Larger `J` fail block space. Thus the old checker cannot certify
this tuple for any positive integer packing count.

## 4. A bound over the full fixed-filter parameter space

This bound assumes exactly the named filters and at least the minimum
reserves above. Put

```text
A=H-1/2+0.00005, B=2H-1+0.00005,
C1=0.5000176, C2=1.0000216.
```

Eliminating `d,delta,gamma` from the coverage, representation, and cost
inequalities gives

```text
tau(1+rho) >= C1+A rho,
tau(1+rho) >= C2+B rho-eta.
```

Any competitive tuple has `H/2<tau<1/2`, using `rho<1/2`, `gamma<1/2`, and
the list bound. Thus `tau>A` and `tau<B`. Solve the two inequalities for
lower/upper bounds on `rho` and eliminate it to obtain

```text
tau >= T(eta) = [A eta+C1 B-A C2]/[eta+C1+B-C2-A].             (10)
```

The denominator is positive on this domain and `T` strictly decreases.
At `eta0=0.517990280`, [frontier_bound.py](frontier_bound.py) evaluates
`tau0=T(eta0)` outward and supports the second log moment by a convex plane.
For fixed scores define

```text
Psi(p,q) = log2 integral M(x)^2 2^(p v(x)-q u(x)) dx.
```

This function is convex. At `p0=lambda0*a0,q0=lambda0`, its gradient is
`(E_tilt[v],-E_tilt[u])`. Writing
`g_lambda=a0 E_tilt[v]-E_tilt[u]` gives the supporting plane

```text
log2 I2(a,lambda)
  >= log2 I2(a0,lambda0)
     +lambda(a-a0) E_tilt[v]+(lambda-lambda0)g_lambda.         (11)
```

The checker bounds `g_lambda` by rigorous left/right convex secants in
`lambda` and bounds `E_tilt[v]` from below by a backward secant in `a`.
Each secant endpoint uses a continuous interval enclosure; these are not
floating-point derivative estimates.

If `alpha<=T(eta0)+0.00001`, (10) and monotonicity force `eta>=eta0` and
`tau<=tau0`. With `a=(eta-tau)/(2tau-eta)`,

```text
a-a0 >= tau0(eta-eta0)/(2tau0-eta0)^2.
```

The certificate bounds `E_tilt[v]*tau0/(2tau0-eta0)^2-1` strictly above zero.
Since admissible multipliers have `1<=lambda<10`, (11) then keeps the
second moment above its ceiling for every such tuple. The certificate uses
an outward upper endpoint for `tau0`, making this test conservative.
[BOUND.json](BOUND.json) retains the supporting point, secants, positive
residual, and all input/code hashes. This bound does not rely on a search
being stationary or on an assumed global optimum of the filters.

## Evidence and replay

[audit.py](audit.py) uses exact integers and fractions. The retained campaign
checks signed, zero, and repeated integer weights, canonical projection
injectivity, exact uniform-cover products, normalized adjacent restriction
inequalities, positive transformed inputs, sharpness controls, the
superadditivity counterexample, exact hypergeometric comparisons, and finite
probability cutoffs. [AUDIT.json](AUDIT.json) records the actual
counts and binds the code, candidates, and predecessor.

From the repository root, with existing requirements installed:

```sh
python3 -B research/fractional-coverage/search.py --author conorholds --fixed-multipliers --output /tmp/fractional-fixed.json
python3 -B research/fractional-coverage/search.py --author conorholds --output /tmp/fractional-retuned.json
python3 -B research/fractional-coverage/interval_check.py --author conorholds --check-receipt research/fractional-coverage/CERTIFICATE.json
python3 -B research/fractional-coverage/audit.py --author conorholds --max-n 8 --output /tmp/fractional-audit.json
python3 -B research/fractional-coverage/frontier_bound.py --author conorholds --output /tmp/fractional-bound.json
python3 -B research/fractional-coverage/test_package.py
make check
git diff --check
```

The incompatible model deliberately uses `candidates.json` and its explicit
checker. The test router includes it in `make check`; the legacy-only
`make verify` command does not discover this model automatically. Hosted
verification remains manual and maintainer-only.

Finite controls do not prove the all-size lemma. Interval checks do not
provide independent review of the full algorithm. The source's prime
estimates, bounded draws, exact exclusion implementation, finite-circle
passage, and final assembly remain separate review dependencies.

## Next question

Can the boundary `delta=d` be treated with a sublinear deficiency allowance
absorbed by the strict representation gap? The current checker deliberately
requires a positive fixed coverage reserve and makes no boundary claim.
