# Checksum sampling without a pattern-count variance loss

The checksum estimator can sum **independent conditional sampling errors**
before averaging over its phases and tables. Combining this with unique
pattern ownership removes the source proof's `Q_n^2` factor from the
mean-squared sampling-error bound. For finite-circle log-moment ceiling `xi`,
the proposed bound is

```text
E[1_G |Vhat_rho - V_rho|^2] <= 8 * 2^[-(eta-xi-48 epsilon)n].       (1)
```

This bound is uniform over the residue input and target. It assumes the
source's coverage and simultaneous list/match-size event `G`, and uses fresh
independent match draws for different patterns. The aggregation makes no
independence assumption about the pattern environments or input residues;
the source's table-construction assumptions still supply `G`. Here `Vhat`
denotes the ideal exact-arithmetic estimate before caps; caps do not interrupt
it on `G`. Coverage/cap failure probabilities and finite arithmetic remain
separate terms. The full running-time exponent is unchanged: the estimator
still constructs its entry lists and processes `Q_n=2^O(sqrt(n))` patterns.
The same calculation also permits a smaller **sampling** exponent; at the
retained `alpha=0.48946580` candidate, it changes from `0.48945580` to
`0.48943092` while preserving exponential decay of the ideal sampling error.

This is a proposed localized proof refinement and implementation audit, not
an independently reviewed theorem or a new Subset Sum running-time record.
All selected metadata and existing certificates are preserved.

## Source, scope and predecessor

Canonical project main is `f25a303525df636c079e6a7c574527cccf56e611`.
The mathematical source is
[`openai/math@1581d4c`](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026).
The exact seven section hashes, sixteen parameter snapshots and source-file
hashes are retained in [INPUTS.json](INPUTS.json). The relevant labels in
`checksum.tex` are:

| Source label | Interface audited here |
| --- | --- |
| `eq:checksum-match-count` | Group products and exact uniform rank indexing. |
| `eq:checksum-multiplicity` | Retesting every indexed block table, retaining repetitions. |
| `eq:checksum-direct-weight`, `eq:checksum-sample-weight` | Unique ownership and inverse-multiplicity weighting. |
| `eq:checksum-conditional-variance` | Conditional complex second moment and independent repeated draws. |
| `eq:checksum-unconditional-domination`, `eq:checksum-error-exponent` | Remove the good-event restriction before averaging under the original phase law. |

This complements the deficiency audit in
[PR #14](https://github.com/conorholds/subset-sum-bounds/pull/14), alias-cost
audit in [PR #15](https://github.com/conorholds/subset-sum-bounds/pull/15), and
filter/coverage work in [PR #16](https://github.com/conorholds/subset-sum-bounds/pull/16).
It does not redo those arguments. The sixteen snapshots cover the four
original witnesses and the tunable-multiplier, fractional-coverage and
asymmetric-filter candidates retained in PR #16.

The exact-integer exclusion used in alias extraction is **not an additional
test on checksum matches**. Alias matches represent pairs of subset records;
checksum matches represent high/low Fourier-frequency digits and table-index
tuples. The source estimates the whole modular-bin checksum and subtracts
retained aliases afterwards in `implementation.tex`. The checksum match
machinery therefore has no subset-pair exclusion obligation.

## 1. Exact compressed counting and weighting

Fix the phases, every table and the deterministic score-approximation rule.
For a pattern `pi`, let `a_nu,b_nu` be its two entry-group sizes for whole
table-index tuple `nu`. Then

```text
L_pi = sum_nu a_nu b_nu.
```

Choose an integer uniformly from `[0,L_pi)`, locate its prefix-sum interval,
then use quotient and remainder within that group. This is a bijection from
ranks to entry pairs. It samples uniformly over **matches**, including
repeated table labels at distinct indices. Uniform selection among groups
would give a different law whenever their products differ.

A frequency has unique high/low digits for that pattern's digit base. Its
number of matches is

```text
mu_pi(k) = product_j |passing_high_j(k) intersection passing_low_j(k)|.
```

Both entry emission and retesting must use the same stored tables, table
indices and deterministic score rule. Labels with identical contents retain
distinct indices. For each sampled or enumerated match, `mu_pi(k)>=1`, even
if coverage of other frequencies failed. Recomputing tests with a different
rounding rule or deduplicating only the multiplicity calculation breaks this
identity.

Let `C_pi` denote the valid frequencies whose deterministic pattern is `pi`.
They partition `[0,P)`. A match for `k>=P` or `k notin C_pi` has value zero.
Conditional on the environment, a single draw has value

```text
X_pi = L_pi J_rho(k)/(P mu_pi(k)) if k in C_pi; zero otherwise.
```

Zero-valued rejected matches still consume a draw. Resampling until a valid
frequency is found changes the law and can bias this weight. An unattained
pattern may have `L_pi>0` and no owned frequencies at all, so such a rejection
loop could also fail to terminate. The fixed raw-draw count avoids that
problem. In the three-point all-star control, conditioning away the extra
digit pair changes the correct estimate mean `1` to `4/3`.

On coverage, `E[X_pi | environment] = v_pi` and

```text
E[|X_pi|^2 | environment]
  = L_pi/P^2 sum_{k in C_pi} |J_rho(k)|^2/mu_pi(k).
```

With `N` independent draws, the exact variance is

```text
Var(Xbar_pi | environment)
  = [L_pi/P^2 sum_{k in C_pi} |J_rho(k)|^2/mu_pi(k) - |v_pi|^2]/N.
                                                                    (2)
```

This uses the complex variance `E|Z-EZ|^2`; conjugation matters. Dropping
`-|v_pi|^2` and using `mu_pi(k)>=1` gives the source bound. If `L_pi=0`,
coverage implies `C_pi` is empty. Without coverage, zero output can be biased;
the proof makes no claim of unbiasedness on failed coverage.

Entry caps and direct-enumeration caps abort before storing or visiting an
extra object. A compressed group product may be much larger than its two
lists; no product expansion is required for indexing or a sampled draw.
The package's grouped oracle and separate full-key reference check this
finite interface, not the source's asymptotic list-construction speed.

## 2. Sum conditional variances before releasing the environment

Write `F` for the phases and tables, and let `G` be the source's event from
`lemma:checksum-tables(iii)`. It is measurable with respect to `F`; it does
not depend on the later match draws. Allocate independent match-draw streams
to different sampled patterns. Conditional on any fixed `F in G`, the pattern
errors have mean zero and independent draws, so

```text
E[|sum_pi E_pi|^2 | F] = sum_sampled_pi E[|E_pi|^2 | F].       (3)
```

Direct patterns have zero error on `G`. Random phases can correlate every
pattern and can make their sampling laws very different; neither fact
affects (3). Independence here concerns only the fresh sampling randomness.
There is no cross-target independence claim or requirement.

Put `Delta=eta-tau`, `Gamma=2tau-eta` and `a=Delta/Gamma`. With tunable
domination slopes `lambda1,lambda2`, define

```text
R1(x)=lambda1(Delta-u(x)),
R2(x)=lambda2(a v(x)-u(x)),
H_r(x)=M(x)^2 2^Rr(x).
```

The source's rounded profile bounds imply its pointwise domination with
these slopes when

```text
lambda1>=1, lambda2>=1,
3 lambda1 <= 30, (4a+3)lambda2 <= 30.                         (4)
```

In the `W_pi<=Gamma` sampling branch, write `t=Delta-U_pi>0`.
Then the difference between `-Delta+lambda1(Delta-ubar)+30epsilon`
and `-U_pi` is at least `(lambda1-1)t+(30-3lambda1)epsilon`.
In the other branch, `alpha W_pi=Gamma` and
`-alpha U_pi=-Delta+alpha(aW_pi-U_pi)`. The gap `aW_pi-U_pi` is
positive; `alpha<=1<=lambda2` and the profile errors cost at most
`lambda2(4a+3)epsilon`. This proves (4)'s sufficiency. The `18epsilon`
match-count allowance plus the `30epsilon` domination allowance gives
`48epsilon`.

For each `k`, there is at most one sampled owning pattern. If it uses branch
`r in {1,2}`, the source's domination gives

```text
1_G sum_sampled_pi 1_{k in C_pi} L_pi |J_rho(k)|^2
  <= 2^[(eta-Delta+48epsilon)n] [product_i H1(x_i(k))
                              + product_i H2(x_i(k))].      (5)
```

This sums over **owning patterns**, not over all patterns that happen to
have a match for `k`. It replaces a possible `Q_n` repetition by two
nonnegative products.

Now average the conditional variance sum (3), bounded using (2) and (5),
under the **original** environment law. At each fixed `k`, the random
`x_i(k)=phi_i+k b_i` are independent uniform residues, whatever the input
`b_i`. If each finite-circle moment satisfies `average H_r<=2^xi`,
each product expectation is at most `2^(xi n)`. Summing over the `P`
frequencies and dividing by `N P^2` gives

```text
E[1_G |sum_pi E_pi|^2]
  <= 2 * 2^[(eta-Delta+48epsilon+xi)n]/(N P).
```

Using `N=2^ceil(tau n)` and `|log2 P-eta n|<=2` proves (1), since
`tau+Delta=eta`. The factor eight is two branches times the factor four
allowed by the modulus interval. In contrast, applying the source's
per-pattern bound and then Cauchy--Schwarz gives a `Q_n^2` factor.

Conditioning the phases themselves on `G` would invalidate the product
expectation. Only the nonnegative left side is masked by `1_G`; (5)'s right
side is averaged without a good-event restriction. This distinction is
essential even on the smallest circle.

## 3. Explicit finite-circle reserve and total failure probability

Suppose a retained continuous certificate bounds both log moments by `chi`,
and put `s=eta-chi-48epsilon>0`. The fixed smooth positive filters have finite
normalizers and scores converging uniformly to their continuous versions,
as in `filtering.tex`. For all sufficiently large `P`, both finite log
moments are bounded by

```text
xi=chi+s/2,
eta-xi-48epsilon=s/2.
```

This spends half the **remaining variance reserve**, rather than silently
treating the continuous ceiling as an already proved finite-circle bound.
It does not require a changed numerical witness or checker. For the retained
asymmetric candidate at `alpha=0.48946580`, the exact values are

```text
s=0.00009952, xi=0.51841070, variance decay=s/2=0.00004976.
```

For `T` targets and ideal absolute-error tolerance `t>0`, Markov and a union
bound give

```text
Pr[G and some target error > t] <= 8 T t^(-2) 2^(-s n/2).    (6)
```

For `t=1/16`, the coefficient is `2048 T`. The overall failure bound also
adds `Pr(G^c)`, bounded-draw coupling failures, and finite arithmetic error.
Those terms are not removed by this refinement. The algorithm still has
`Q_n` patterns in its work and coverage/cap union bounds; only the sampling
error estimate loses the `Q_n^2` factor. No practical crossover size is
established.

### Fewer match draws, with an explicit remaining error reserve

The calculation before substituting `N=2^ceil(tau n)` gives the more general
bound, for any common integer sample count `N>=1`,

```text
E[1_G |sum_pi E_pi|^2] <= 8 * 2^[(xi-Delta+48epsilon)n] / N.   (7)
```

Keep `xi=chi+s/2`, and choose

```text
sigma=tau-s/4 > 0,
N=2^ceil(sigma n).
```

Then `sigma+Delta-xi-48epsilon=s/4`, so (7) is at most
`8*2^(-s n/4)`. At the asymmetric candidate,

```text
sigma=0.48943092,
remaining decay=s/4=0.00002488.
```

This is a factor `2^(0.00002488 n+O(1))` fewer draws per sampled pattern,
paid for by half of the finite-circle variance reserve left in (1).
The source's entry-list construction and direct-pattern cap remain at the
old budgets, so this does not lower `alpha` or the complete algorithm's
leading exponent. It is a proposed subroutine change, not a benchmark of an
implemented full checksum estimator. The sixteen exact snapshots all have
`sigma>0`; the tests check their full reserve equations without floating
point arithmetic.

## 4. A stronger uniform multiplicity event

The source's change-of-measure coverage proof gives, conditional on arbitrary
phases and for a frequency in its own pattern, a per-row acceptance probability
on block `j` of at least

```text
p_j >= (1/2) 2^[-(alpha D_j+2epsilon)d_j].
```

Its `R_j=2^ceil[(alpha D_j+4epsilon)d_j]` independent indexed rows therefore
give a binomial acceptance count `Z_j` with mean

```text
m_j=R_j p_j >= 2^(2epsilon d_j-1).
```

For any binomial count `Z` of mean `m`, exponential Markov with `t=ln 2`
shows

```text
Pr[Z<m/2]
  <= exp[(t/2-(1-exp(-t)))m]
   = exp[-(1-ln 2)m/2]
  <= exp(-m/8).
```

Here `ln 2<3/4`, for example because the first four positive terms of
the series for `exp(3/4)` already exceed two. Hence, uniformly in the phases,

```text
Pr[Z_j < 2^(2epsilon d_j-2)] <= exp[-2^(2epsilon d_j-4)].
```

A union bound over the `q` blocks and `P` valid frequencies gives a failure
probability at most

```text
q P exp[-2^(2epsilon d_min-4)] = o(1),
```

since `d_min=Theta(sqrt(n))`. No union over unattained patterns is needed.
Intersect this event with the original `G` to form `G*`. It is still an
environment event of probability `1-o(1)`, and every valid frequency in its
own pattern now satisfies

```text
mu_pi(k) >= product_j 2^(2epsilon d_j-2)
          = 2^(2epsilon n-2q).                               (8)
```

Retaining the exact `1/mu_pi(k)` factor in (2) instead of dropping it gives,
with the original sample count,

```text
E[1_G* |sum_pi E_pi|^2]
  <= 8 * 2^(2q) * 2^[-(eta-xi-46epsilon)n].                  (9)
```

Thus some of the source's epsilon allowance can also be recovered without
new tables or extra work. Equation (9) preserves an explicit `2^(2q)` factor;
it would be incorrect to erase it at finite size. For all sufficiently large
`n`, it improves the linear error exponent because `q=o(n)`. This is a
strengthening of the analysis, not permission to reduce any retained
checker reserve. The reduced-draw proposal above deliberately uses the
simpler bound (7), with every existing reserve unchanged.

## 5. An effective finite-circle component cutoff

For the exact retained asymmetric filter, [FINITE-CIRCLE.md](FINITE-CIRCLE.md)
gives a self-contained Lipschitz transfer argument. The existing continuous
ceiling `chi=0.51836094` implies the finite ceiling `xi=0.51841070` for every
integer modulus **`P>=2^351`**. The modulus sizing `|log2 P-eta n|<=2` then
makes **`n>=681`** sufficient for this one component.

[FINITE-INPUT.json](FINITE-INPUT.json) binds the exact coefficients and
predecessor; [finite_circle.py](finite_circle.py) calculates only rational
upper bounds and integer powers of two, retained in [FINITE.json](FINITE.json).
A separate calculation with slightly looser rational envelopes recovers the
same modulus power. The argument replaces all finite normalizers, handles
sine terms, and avoids taking a logarithmic derivative of the vanishing
Fourier-amplitude factor.

This cutoff is deliberately conservative and does not make the complete
algorithm valid or practical at 681 inputs. Coverage, cap failures, bounded
draws and finite arithmetic have their own thresholds. The continuous
certificate and independent review of the new transfer proof remain
dependencies.

## Exact evidence and replay

[oracle.py](oracle.py) implements bounded grouped counting, prefix-rank
selection, literal table retesting, and exact character arithmetic on
circles of size 2, 3, 4 and 6. [reference.py](reference.py) independently
enumerates whole table keys and subset expansions; both implementations
share the exact scalar arithmetic and assisted authorship. This is an
internal cross-check, not independent scientific review.

The finite tests use positive rational densities and rational **product**
thresholds to make acceptance exact. They do not implement the source's
entropy-profile grid, approximate logarithm rule, random table-count law or
all-size coverage experiment. Their purpose is to audit the combinatorial
and sampling identities once the acceptance relation is fixed.

From the repository root:

```sh
python3 -B research/checksum-variance/audit.py --author conorholds --output /tmp/checksum-standard.json
python3 -B research/checksum-variance/audit.py --author conorholds --full --output /tmp/checksum-full.json
python3 -B research/checksum-variance/finite_circle.py --author conorholds --output /tmp/checksum-finite.json
python3 -B research/checksum-variance/test_package.py
git diff --check
```

The retained audit includes explicit biased controls for omitted
multiplicity, inconsistent label deduplication and uniform-group sampling;
a shared-rank covariance control; and a coverage-conditioned phase-law
control. Exact replay counts and input/code hashes are in [AUDIT.json](AUDIT.json).
The full bounded campaign passed **163,120 configurations**, **655,120**
Fourier/target identities, **326,240** exact pattern second moments and
**184** independent joint-draw controls. It compared **677,854** literal
matches with compressed indexing, including **6,906** out-of-range matches.
All residue/phase vectors and all pairs of indexed one-coordinate labels
on each of two blocks were exhausted on the four-point circle; the standard
campaign also covers circles of size two, three and six. This finite domain
does not assert exhaustive coverage for larger circles or block dimensions.

Some retained failure values are particularly small:

| Control on the four-point circle | Exact estimate mean, with correct mean `1` |
| --- | --- |
| Omit inverse multiplicity | `9/4 + i/4` |
| Deduplicate labels only when retesting multiplicity | `15/8 - i/8` |
| Choose a table group uniformly, then choose a pair | `15/16 + 3i/16` |

The correct estimator has one-draw variance `17/16` in that fixture. In a
separate two-pattern fixture, independent ranks give total MSE `3/2`, while
reusing the same rank gives `2`. On the two-point circle, a coverage event
of probability `1/4` changes the amplitude second moment from `2` to `4`
and the joint score moment from `1` to `9/4`. This is a finite illustration
of the conditioning hazard, not a counterexample to the source's masked
original-law argument.

In an all-star four-point fixture, target two has zero checksum, unit
multiplicity and one-draw variance exactly `2`, attaining the upper bound
`L/P^2 sum_k |J(k)|^2`. Thus the conditional second-moment inequality used
above can be sharp in this finite interface; a general gain from multiplicity
alone requires additional information. This is not a lower bound on the
source algorithm's asymptotic exponent.

Exact binomial-tail tests independently check the dyadic consequence
`Pr[Z<m/2] <= 2^(-m/8)` for `1<=R<=64` and seven rational probabilities.
They do not establish the all-size concentration theorem; the exponential
Markov argument above supplies that step.

Not established: an independent review of the written proofs, a formalized
theorem, explicit all-size coverage/cap or finite-arithmetic cutoffs,
full finite-arithmetic implementation, or a faster complete
Subset Sum solver. This package introduces no new filtering certificate.
The pinned moments, source coverage/size theorem, bounded random draws,
precision argument and final assembly remain dependencies.

## Credit

OpenAI supplies the source checksum algorithm, table construction and
inverse-multiplicity estimator. The changes here are the conditional
variance aggregation, unique-owner domination and focused exact audit.
Preceding parameter and filter contributions retain their source bindings
in [INPUTS.json](INPUTS.json). Substantial AI assistance: Codex.
Apache-2.0 licensing is retained.

## Next question

Can a quantified coverage/cap failure bound make the sharper error formula
useful at a concrete input size? Start by fixing the source's block partition
and inserting explicit score bounds into its change-of-measure concentration
step, preserving every rounding reserve. The finite-circle component now has
a cutoff, but the other thresholds and a practical crossover remain open.
