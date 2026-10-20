# Effective checksum coverage and cap bounds

The checksum's coverage and simultaneous list/match-size event can be given
an explicit failure bound. A small change to the block partition avoids
requiring square-root-sized blocks before concentration is available.
For the retained asymmetric filter, the sufficient component cutoffs are:

| Partition | Sufficient input size | Failure bound for coverage, multiplicity and caps |
| --- | --- | --- |
| Fixed `q=ceil(sqrt(n))` | `n>=2^122` | `2^-128` |
| Enforce a minimum block size | `n>=2^61` | `2^-128` |

These are conservative sufficient bounds from the argument below, **not
optimal cutoffs or practical solver sizes**. Both are enormous. The second
partition eventually becomes the source's square-root partition and
preserves its asymptotic operation bound. Near its smaller cutoff, its
individual tables can be much larger. No faster complete algorithm or new
leading exponent is claimed.

Composing this event with the companion
[checksum sampling refinement](../checksum-variance/README.md) gives a
conditional **ideal-arithmetic** statement: for `n>=2^61`, at most `n` targets,
and the reduced sampling exponent `0.48943092`, all target errors are at most
`1/16` with failure probability at most `2^-127`. This assumes the pinned
continuous moment ceiling and the companion finite-circle/variance proofs.
The source's bounded-draw coupling and arithmetic errors are quantified in
[IMPLEMENTATION.md](IMPLEMENTATION.md). Allocating `2^-130` to each random
failure term gives checksum error at most `1/8` simultaneously for at most
`n^2` targets, with failure at most `2^-128`, above the same probability
cutoff. At `n>=2^93`, an additional uniform bound pays the table-scan and
pattern overhead with an exponent of at most `epsilon n`.
These remain conditional checksum-component results. The full polynomial
operation count, machine-address accounting and the remaining Subset Sum
components are not quantified here; independent mathematical review is pending.

The explicit term counts also permit smaller deterministic implementation
budgets: **linear phase/contribution precision instead of quadratic precision,
and a linear rejection cap for ordinary uniform draws instead of `n^3` attempts**. Their exact
schedules retain the same error allowances on all capped histories. They
reduce local precision and rejection budgets; no complete runtime gain is
measured or claimed.

## Source and scope

The source is
[`openai/math@1581d4c`](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026),
`filtering.tex` (`eq:filter-scores`, `lemma:finite-filtering`),
`checksum.tex` (`eq:checksum-pattern-order`, `eq:checksum-coverage-probability`,
`eq:checksum-whole-table-count`, `lemma:checksum-tables`), and
`implementation.tex` (`eq:sampling-coupling`, `eq:finite-evaluation-error`).
Canonical project main is `f25a303525df636c079e6a7c574527cccf56e611`.
[INPUTS.json](INPUTS.json) binds the three sections and companion files;
the companion's complete seven-section source manifest is also hash-bound.
The exact filter snapshot comes from the asymmetric candidate in
[PR #16](https://github.com/conorholds/subset-sum-bounds/pull/16).

This resolves the explicit coverage/cap question left by the checksum
sampling package, and quantifies the source's subsequent checksum coupling
and rounding bounds, including smaller precision/rejection schedules.
It supplies no new numerical filtering certificate and
changes neither the selected candidate nor existing reserves. The source's
uniform finite-score evaluator remains a dependency; the bounds calculator
does not implement that evaluator or construct exponentially many tables.
Fix the residue input and target list before drawing the checksum randomness,
as in the source's modular-checksum theorem.

## 1. Score widths give a sufficient block size

All logarithms of densities are base two. Natural `exp` is used for tail
bounds. Let `R,S` be the sums of absolute cosine and sine coefficients of
the real trigonometric polynomials `r,s`; write `B=R+S`.
On every finite circle, cancellation of the normalizers gives

```text
log K_x(w) = r(x+w)+s(w)-log average_w 2^[r(x+w)+s(w)],
log f_*(y) = kappa r(y)-log average_y 2^[kappa r(y)],
log g_*(z) = kappa s(z)-log average_z 2^[kappa s(z)].
```

Their ranges have lengths at most

```text
w_D=2(R+S), w_A=2 kappa R, w_E=2 kappa S.                 (1)
```

Each normalized positive density has values both at most and at least one,
so zero lies in each logarithm's range. Adding the star score zero therefore
does not increase these widths. This matters: using twice an absolute
bound as the width would unnecessarily lose another factor four.
We also have `0<=D<=2B`, `|A|<=w_A`, and `|E|<=w_E`.

Under the source's changed label law, coordinate scores remain independent,
even conditional on arbitrary phases and frequency. Their means are
`alpha A(x)`, `alpha E(x)` and `alpha D(x)`. Only three **one-sided**
deviations need excluding: a low high-side score, a low low-side score,
and a high log likelihood.

For completeness, the bounded-variable estimate used here is the usual
[Hoeffding bound](https://www.tandfonline.com/doi/abs/10.1080/01621459.1963.10500830).
For a variable in an interval of width `w`, the second derivative of its
centered log moment-generating function is a tilted variance, at most
`w^2/4`. Integrating twice gives `log E exp[t(X-EX)]<=t^2 w^2/8`.
Independence, exponential Markov and `t=4epsilon/w^2` give a one-sided
average-deviation probability at most
`exp(-2epsilon^2 d/w^2)`. A constant variable has zero deviation probability.
The coordinates need not have identical laws or share interval endpoints.

The first four terms of `exp(7/10)` sum to `12013/6000>2`.
Consequently it suffices to require

```text
d >= D0 = max(1,
             ceil[7 w_D^2/(10 epsilon^2)],
             ceil[21 w_A^2/(20 epsilon^2)],
             ceil[21 w_E^2/(20 epsilon^2)]).             (2)
```

The likelihood upper tail is then less than `1/4`, and each score lower
tail less than `1/8`. With probability at least `1/2` under the changed
law, the two exact score averages are at least their means minus `epsilon`
and the log likelihood average is at most its mean plus `epsilon`.
The source's profile error and score error are each at most `epsilon`;
thus the same `-3epsilon` tests pass and the block likelihood ratio is at
most `2^[(alpha D_j+2epsilon)d]`. No reserve is reduced.
Changing measure back gives the source's explicit per-row bound

```text
p_j >= (1/2) 2^[-(alpha D_j+2epsilon)d_j].                (3)
```

For the retained filter, exact rational arithmetic gives

```text
D0=2006670879202202015,
individual sufficient sizes:
  2006670879202202015, 809923614872995342, 801885223128037735.
```

The first term controls this estimate. Equation (2) is sufficient, not
necessary; no lower bound of this magnitude is asserted.

## 2. Make the grid and its pattern count explicit

Put `g=epsilon/4`. Obtain deterministic approximations to each block
average with error at most `epsilon/8`. Round the approximation's upper
endpoint upward to a multiple of `g` for `D_j`, and its lower endpoint
downward for `A_j,E_j`. Then

```text
0 <= D_j-Dbar <= epsilon/2,
0 <= Abar-A_j <= epsilon/2, 0 <= Ebar-E_j <= epsilon/2.
```

In particular the source's `epsilon` error allowance holds,
`D_j>=0`, and `A_j,E_j<=D_j`. Choose the fixed grid only within

```text
D_j in [0,2B+epsilon/2],
A_j in [-w_A-epsilon/2,w_A], E_j in [-w_E-epsilon/2,w_E].
```

Every rounded actual profile lies in this grid; enumerate all its ordered
triples, including unattained patterns. Conservative upper counts are

```text
M_D = ceil[(2B+epsilon/2)/g]+1,
M_A = ceil[w_A/g]+ceil[(w_A+epsilon/2)/g]+1,
M_E = ceil[w_E/g]+ceil[(w_E+epsilon/2)/g]+1.
Gbits = ceil(log2 M_D)+ceil(log2 M_A)+ceil(log2 M_E).
Q_n <= 2^(Gbits q).                                      (4)
```

The ceilings deliberately overcount possible boundary indices. For the
retained filter the counts are `6772500498`, `7026149213`, `6991195439`;
each is at most `2^33`, so **`Gbits=99`**. The grid and collection are
fixed independently of the phases. Averaging a joint score after
conditioning on a frequency's pattern would not justify the match bound.

## 3. Coverage, multiplicity and caps with every prefactor retained

Let the balanced partition have `q` blocks and minimum size `dmin>=D0`.
The source draws `R_j=2^ceil[(alpha D_j+4epsilon)d_j]` independent indexed
rows, so (3) gives an acceptance count of mean at least
`2^(2epsilon d_j-1)`. As proved in the companion package, the binomial
lower-tail estimate yields

```text
Pr[some valid frequency/block has fewer than 2^(2epsilon d_j-2)
   accepting indexed rows in its own pattern]
  <= q P exp[-2^(2epsilon dmin-4)].                       (5)
```

There is no union over patterns in (5): after fixing the phases, each
frequency has exactly one owning pattern. The bound is uniform in those
phases, so it can be averaged without changing their original law.
Outside this failure, every frequency is covered and has multiplicity
at least `2^(2epsilon n-2q)` in its own pattern.

For caps, retain the source's parameter choices `alpha,H` and all its
score thresholds. From `P<=4*2^(eta n)`, `H>=0`, and
`2^[(eta-H)n]<=B_pi<=2*2^[(eta-H)n]`, the high digit count is at most
`5*2^(Hn)`, the low count at most `2*2^[(eta-H)n]`, and their product
at most `6*2^(eta n)`. The whole-table count is at most
`2^[(alpha D_pi+4epsilon)n+q]`. Thus the source's marginal and
original-phase joint score moments give

```text
E high entries <= 5*2^[(tau+8epsilon)n+q],
E low entries  <= 2*2^[(tau+8epsilon)n+q],
E matches      <= 6*2^[(eta-alpha U_pi+12epsilon)n+q].
```

The entry bounds even hold conditional on phases; the match bound uses
their original independent uniform law. No independence across patterns
or whole-table tuples is needed. Markov at the source's good-event caps
(`tau+14epsilon` for both entry lists and `eta-alpha U_pi+18epsilon`
for matches), followed by (4), gives

```text
Pr[some good-event cap fails] <= 13*2^[(Gbits+1)q-6epsilon n]. (6)
```

Combining (5) and (6) bounds the complement of an environment event `G*`
contained in the source's coverage/size event. Sampling draws are not used
to define it. On `G*`, operational entry and direct-enumeration caps
`2^ceil[(tau+20epsilon)n]` do not interrupt the ideal estimator.

## 4. Two explicit partitions and bounds for every larger n

Let `L=ceil(log2 n)` and set

```text
t(n)=max(D0, ceil[(L+7)/(2epsilon)], ceil[(Gbits+1)/(2epsilon)]).
q(n)=min(ceil(sqrt(n)), floor[n/t(n)]).                    (7)
```

Use (7) only when `n>=t(n)` and distribute coordinates as evenly as
possible. Then `dmin>=t(n)` and `(Gbits+1)q<=2epsilon n`.
Also `2epsilon dmin-4>=L+3`, so `exp(-x)<=2^-x`, `eta<=1`,
`q<=n`, and `ceil(log2 n)<=n` make (5) at most `2^-n`.
Together,

```text
Pr[(G*)^c] <= 2^-n + 13*2^(-4epsilon n).                 (8)
```

The table-generation tradeoff is explicit. Always
`q<=ceil(sqrt(n))` and `dmax<=max(ceil(sqrt(n)),2t(n))`.
Since `t(n)=O(log n)` for the fixed filter and epsilon, eventually (7)
chooses the source's exact block count. Thus its asymptotic
`2^[(tau+20epsilon)n+O(sqrt(n))] poly(n+T)` form is preserved.
Before that point, larger blocks can increase table scanning sharply.
An explicit per-digit table-scan overhead exponent is

```text
ceil(log2 q)+1+(2B+epsilon/2+4epsilon)dmax,
```

with a further pattern exponent at most `Gbits q`. [BOUNDS.json](BOUNDS.json)
retains these quantities. The smaller cutoff in (8) is a probability
threshold; it is not a finite runtime improvement.

The calculator checks sufficient conditions at `N=2^k` that imply the
bounds for **every integer `n>=N`**, not just powers of two. For (7), these
are `N>=D0`, `2epsilon N>=k+8`, `2epsilon N>=2`,
`2epsilon N>=Gbits+1`, `4epsilon N>=b+5`, and `N>=b+1` for desired
failure `2^-b`. The first four imply `n>=t(n)` uniformly: write `n=Nx`,
use `log2 x<=2(x-1)`, and `ceil(log2 n)<=k+log2 x+1`.
The last two make each term in (8) at most `2^(-b-1)`.

For comparison, with fixed `q=ceil(sqrt(n))`, use
`dmin>=sqrt(n)-2` and `q<=2sqrt(n)`. Sufficient conditions at `N=2^k`,
`u=floor(sqrt(N))`, are

```text
u-2>=D0, u>=2(Gbits+1)/epsilon, 5epsilon N>=b+5,
2epsilon u>=4, 2epsilon(u-2)>=k+7, 2N-2>=b+1.
```

They give `(Gbits+1)q<=epsilon n` and make the cap bound at most
`2^(-b-1)`. For coverage, `log2 x<=4(sqrt(x)-1)` shows
`2epsilon dmin>=ceil(log2 n)+6` uniformly above `N`.
Thus the exponential in (5) is at least `4n`, and (5) is at most
`2^(2-2n)<=2^(-b-1)`. The two elementary logarithm inequalities follow
from `ln 2>1/2`, `ln x<=x-1` and `ln x<=2(sqrt(x)-1)`.

For `b=128`, exact checks first satisfy these particular sufficient
conditions at powers **61** and **122**, respectively. The preceding
dyadic powers fail the concentration condition. These are minima for
these stated conditions, not optimal mathematical thresholds.

## 5. Compose the ideal checksum error

The companion's finite-circle calculation transfers the pinned ceiling
`chi=0.51836094` to `xi=0.51841070` at `P>=2^351`, sufficient for
`n>=681`. Its independent-pattern sampling proof allows
`Ndraw=2^ceil(0.48943092 n)` and gives

```text
E[1_G |Vhat-V|^2] <= 8*2^(-0.00002488 n).
```

Because `G*` is an environment event contained in `G`, the same upper
bound applies with `1_G*`. Markov at error `1/16` and a union over at
most `n^d` targets gives `2048 n^d 2^(-0.00002488 n)`.
At `N=2^61`, writing `delta=0.00002488`, the exact inequalities
`delta N>=d(61+1)+128+11` and `delta N>=2d` hold for `d=1`.
Using `log2 x<=2(x-1)` proves sampling failure at most `2^-128`
for every `n>=N`; they also hold for every integer `1<=d<=1024`.
Adding (8)'s failure gives **`2^-127`** for the ideal stage.
Fresh independent match draws across patterns and the pinned continuous
moment ceiling are required. This composition does not account for
bounded-draw coupling, numerical approximation error in the final estimate,
or the rest of the Subset Sum algorithm.

## Exact evidence and replay

[coverage_bounds.py](coverage_bounds.py) uses only exact rational arithmetic,
integer roots and exponent comparisons; it never constructs a probability
power with an astronomically large exponent. It imports the companion's
exact candidate parser and coefficient bounds. All three receipts bind their
inputs and code. The predecessor inputs remain unchanged.

[coverage_audit.py](coverage_audit.py) provides a separate exact finite
change-of-measure control on `Z/3Z` with `f=(1/2,1/2,2)`, `g=1`, `kappa=1`.
Its dynamic program exhausts the original and tilted score distributions
for `alpha` in `{1/4,1/2,1}` and block sizes `{288,320,384}`.
It checks every likelihood identity, tilted typical mass at least `1/2`,
the worst-case score-approximation allowance, and the transferred acceptance
lower bound. This synthetic law is not the degree-twenty filter.
A one-coordinate control has typical mass zero, illustrating why a block
size hypothesis is needed. Small literal label enumerations independently
check the dynamic program.

From the repository root:

```sh
python3 -B research/checksum-coverage/coverage_bounds.py --author conorholds --output /tmp/coverage-bounds.json
python3 -B research/checksum-coverage/coverage_audit.py --author conorholds --output /tmp/coverage-audit.json
python3 -B research/checksum-coverage/implementation_bounds.py --author conorholds --output /tmp/checksum-implementation.json
python3 -B research/checksum-coverage/test_package.py
git diff --check
```

The tests include outward grid rounding at both approximation-error
endpoints, balanced partition transitions, non-dyadic input sizes,
independent coarser coefficient envelopes, preceding-power failures,
explicit coverage/cap failures, exact binomial tails and malformed inputs.
Written arguments supply the all-size implications. This is an internal
check, not independent scientific review or a formal proof.

## Credit and next question

OpenAI supplies the source algorithm, coverage change of measure and
list/match moment arguments. The new work makes their constants explicit,
chooses a minimum-block partition, composes the checksum error budget and
quantifies table/pattern overhead.
Substantial AI assistance: Codex. Apache-2.0 licensing is retained.

Next: quantify the polynomial operation count of the deterministic finite-score
evaluator and its storage/address requirements, so that the exponential-factor
bound can be turned into a concrete word-RAM resource bound.
