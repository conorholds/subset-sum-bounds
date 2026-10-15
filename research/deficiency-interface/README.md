# Preliminary search and representation coverage: a scoped proof audit

This audit supplies a parameterized argument for one interface in the
conditional model: **preliminary search either resolves the high-deficiency
case with probability `1-O(1/n)`, or the low-deficiency case supplies enough
representations for alias recovery**. It checks the original parameter tuple
and seven other distinct retained tuples, including the proposals in
[PR 12](https://github.com/conorholds/subset-sum-bounds/pull/12) and
[PR 13](https://github.com/conorholds/subset-sum-bounds/pull/13).

Two concrete refinements accompany the argument:

- An explicit random-restriction variance bound, asymptotically sharp at `n/16` for
  halves, sharpens the source's bound `q` for a restriction of size `q`.
- Under the existing strict coverage gap, a balanced representation block is
  useful with probability **`1-O(1/n)`**, strengthening the source's positive
  constant obtained from Markov's inequality.

There is also an exact positive-weight counterexample to a tempting incorrect
inference: a preliminary `CONTINUE` does **not** certify low expected deficiency.
Its uniform-half cap-overflow probability is `1/143` despite high expectation.
A scalable family covers the actual retained thresholds as well.

This is a proposed proof audit with exact finite controls, not a full
algorithmic review, a machine-checked proof, or a new exponent. The selected
result and its pending independent-review status are unchanged.

## Source and scope

Canonical main is pinned at
[`f25a303`](https://github.com/conorholds/subset-sum-bounds/tree/f25a303525df636c079e6a7c574527cccf56e611).
Coordination is in
[Discussion 5](https://github.com/conorholds/subset-sum-bounds/discussions/5#discussioncomment-18950854).
The manuscript is pinned by [SOURCE.json](../filter-budgets/SOURCE.json) to
`openai/math@1581d4cbdd69a0faf34b2089614fe2a814316a6a`.

The relevant source obligations are:

| Section | Labels used | Role in this audit |
| --- | --- | --- |
| [setup.tex](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/build/sections/setup.tex) | `lem:deficiency`, `lem:half-variance`, `prop:preliminary` | Support identities, restriction laws, capped distinct-sum search. |
| [aliases.tex](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/build/sections/aliases.tex) | `eq:alias-half-block`, `eq:alias-useful-block`, `eq:alias-many-weights`, `eq:alias-collision-moment` | Balanced block, useful deficiency, cardinality layer, smaller-modulus hit. |
| [implementation.tex](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/build/sections/implementation.tex) | fixed-input case split after `eq:residual-accuracy` | Combines unconditional success events without conditioning on `CONTINUE`. |

[parameters.json](parameters.json) retains exact tuples, source-package hashes,
and the pinned project commits. The original four-witness package has three
distinct tuples; PR 12 adds three and PR 13 adds two. Filters are immaterial to
this interface, so this package deliberately contains no numerical witness file.

## 1. Deficiency and a sharper variance bound

For fixed arbitrary integer weights define
`Sigma(U)={sum c_i : i in X, X subset U}` and
`D(U)=|U|-log2|Sigma(U)|`. Direct support identities imply:

```text
D(empty)=0;  D(U)>=0;
0 <= D(U union {i})-D(U) <= 1;
D(U union V) >= D(U)+D(V) for disjoint U,V;
|D(U)-D(U')| <= 1 for a one-index exchange at fixed size.
```

These follow respectively from a nonempty support of at most `2^|U|`, a union
with one translate, a sumset of cardinality at most the product of the two
supports, and comparison with the common set after removing the exchanged
index. Coupling consecutive uniform restriction sizes also shows that
`E D(V_q)` is nondecreasing in `q`.

More generally let `f` be any function on `q`-subsets that changes by at most
one under an exchange. Expose a uniform ordered sample without replacement and
its Doob martingale `M_r`. Given a prefix of length `r-1`, couple choices `i,j`
for the next entry by swapping those names in the remainder. If that remainder
contains the other name, the final sets coincide; otherwise they differ by
one exchange. The former event has probability `(q-r)/(n-r)`. Therefore the
conditional values of `M_r` lie in an interval of length at most
`(n-q)/(n-r)`.

A variable in an interval of length `L` has variance at most `L^2/4`:
`E[(X-a)(b-X)]>=0` gives `Var(X)<=(E X-a)(b-E X)<=L^2/4`.
Orthogonality of martingale differences now gives

```text
Var f(V_q) <= C(n,q)
C(n,q) = (n-q)^2/4 * sum_{k=n-q}^{n-1} 1/k^2.                 (1)
```

The endpoint cases `q=0,n` have variance zero. The same argument applied to
`g(U)=f([n] minus U)` gives the complementary bound `C(n,n-q)`; it does not
assume that `g` is a deficiency function. For `1<=q<n`, integration of `1/x^2`
after its first term gives

```text
C(n,q) <= A(n,q) = 1/4 + (n-q)(q-1)/(4(n-1)).                (2)
V(n,q) = min(A(n,q), A(n,n-q)) <= (n+3)/16.
```

Use `V(n,0)=V(n,n)=0`. The last inequality follows from
`(n-q)(q-1)<=(n-1)^2/4`; also `A(n,q)<=q/4`. In particular the half-restriction bound is
`n/16+O(1)`, rather than `n/2+O(1)`. These bounds apply to `D` by the exchange
identity above. The computation uses the closed rational bound `V` for large
`n`, and the sharper finite sums for small controls.

The leading constant is sharp for arbitrary integer weights. For even `n`,
take `n/2` zero weights and `n/2` distinct powers of two. Every restriction
has deficiency equal to its number of zero positions, since its other subset
sums are all distinct. On a uniform half this hypergeometric count has
variance `n^2/(16(n-1))`: each zero-position indicator has mean `1/2`, variance
`1/4`, and distinct indicators have covariance `-1/(4(n-1))`. Thus no uniform
bound `c n+O(1)` with `c<1/16` can hold in this class.

## 2. The high-deficiency preliminary case

Let the preliminary exponent be `p`, the expectation threshold be `d n`, and
assume `g=p-1/2+d>0`. Cap each distinct-sum list at `2^ceil(p n)` as in the
source. Intermediate supports never exceed a side's final support, so the
procedure overflows precisely when a side's final support exceeds that cap.
On a completed search the answer is exact; its all-history cost is
`2^(p n) poly(n+b)`.

Suppose `E D(V_floor(n/2))>d n`. Monotonicity gives this lower expectation
bound for both marginal side sizes `q=floor(n/2),ceil(n/2)`. Put

```text
t_q = d n-q+ceil(p n) >= g n-1/2.
```

For positive `t_q`, overflow implies a downward deviation from the mean of
more than `t_q`. Cantelli's inequality bounds its probability by
`V(n,q)/(V(n,q)+t_q^2)`. For completeness, for a centered variable of variance
`v`, Markov applied to `(X+v/t)^2` bounds `Pr(X>=t)` by `v/(v+t^2)`; apply the
same fact to `-X` for the downward tail.

Union over the two sides and add the source's capped uniform-sampler failure
`n*2^(-n^2)`. This is `O(1/n)` uniformly in the fixed weights because
`V=O(n)` and `t_q` grows linearly. The exact replay conservatively replaces
the sampler term by `1/n`, and returns the vacuous upper bound one if a
finite-size separation is nonpositive.

## 3. The low-deficiency representation case

Now suppose `E D(V_floor(n/2))<=d n`. Let `m=2 floor(rho n/2)`, `q=m/2`,
and let `J` be a positive integer with `J rho<1`. A uniform half can be
partitioned into `J` disjoint uniform-marginal `q`-subsets and a remainder;
indeed `J q<=floor(n/2)` follows directly from the floor and `J rho<1`.
Superadditivity and nonnegativity imply

```text
E D(S) <= d n/J                                                 (3)
```

for a uniform `q`-subset `S`; independence between the partition pieces is
unnecessary.

Fix the large prime and an alias `I`. The trial's independent uniform flip
`F` makes `I triangle F` a uniform subset, regardless of `I`. For independent
uniform `m`-subset `M`, the balance event has probability
`binom(m,m/2)/2^m >= 1/(m+1)`. Conditional on balance,
`S=(I triangle F) intersect M` is uniform over all `q`-subsets of `[n]`.
For any specified `S`, both the number of extensions to `M` and the number of
compatible flips are the same; this also proves the conditional law directly.

Assume the existing strict coverage inequality `J rho delta>d`. With
`a_n=delta m-d n/J`, equations (1)--(3) and Cantelli give, when `a_n>0`,

```text
Pr[D(S)<=delta m | balance] >= a_n^2/(V(n,q)+a_n^2).             (4)
```

Here `a_n=(delta rho-d/J)n-O(1)`, so the lower bound is **`1-O(1/n)`**.
The source's Markov bound `max(0,1-d n/(J delta m))` remains independently
valid; the replay takes the larger of the two. This strengthening uses the
same strict inequality and requires no change to the numerical checker.

An explicit sufficient cutoff for probability at least one half is

```text
G = J rho delta-d > 0
n >= max(ceil(4 J delta/G), ceil(rho J^2/(2 G^2))).              (5)
```

Indeed `m>=rho n-2` then gives `a_n>=G n/(2J)`, and
`V(n,q)<=q/4<=rho n/8` makes `a_n^2>=V` at this cutoff.

## 4. Distinct masks, smaller moduli, and probability laws

On the useful-block event the support has at least `2^((1/2-delta)m)` values.
For every realized flip, signed support is an exact translate of unsigned
support: replacing `X` by `X triangle (F intersect S)` translates its sum by
`-sum_{i in F intersect S}c_i`. This is a pointwise identity and needs no
independence of `S` and `F` after conditioning on balance.

Condition now on the source's sparse-distinction property of the large prime,
which holds simultaneously for all blocks and signs. Distinct signed integer
sums remain distinct modulo that prime, hence also as signed reduced integer
sums. One cardinality layer therefore has at least

```text
B >= 2^((1/2-delta)m)/(q+1) >= 2^((1/2-delta)m)/(m+1)
```

distinct left weights, each with its disjoint complementary right mask.
The strict representation gap `r=1/2-delta-gamma>0` ensures `B>=2^(gamma m)`
for sufficiently large `n`. The replay supplies a conservative effective
cutoff: find a dyadic `m0=2^k` with `r m0>=k+1`, then any `m>=m0` satisfies
`r m>=log2(m+1)`. Doubling preserves the inequality, and it holds throughout
each intervening dyadic interval. Taking `n>=ceil(m0/rho)` ensures `m>=m0`.

The useful layer and these weights depend on `P,M,F,I`, but not on that
layer's independent smaller prime and residue. Their distinct nonzero
differences have `O(n)` bits. Importing the source's random-prime collision
lemma, the expected squared residue occupancies are at most
`B+K n^2 B^2 2^(-gamma m)`. Absorb the first term using `B>=2^(gamma m)`.
Markov and Cauchy--Schwarz give probability at least one half of at least
`2^(gamma m)/(2 K' n^2)` occupied residues; a uniform smaller-modulus residue
hits these with probability at least `1/(4 K' n^2)`. The precise absolute
constant is not computed by this audit.

Thus, above the component cutoffs, a trial recovers a fixed alias with
probability at least `1/(16 K' (m+1)n^2)=Omega(n^-3)`, conditional on a
sparse-distinction prime. Every layer is tried, so the algorithm need not
know the useful layer. Independent repetitions `n^10` still suffice to cover
all at most `2^n` aliases and polynomially many targets with probability
`1-o(1)`. This argument uses the uncapped experiment. Add sparse-prime and
cap-failure probabilities by a union bound, using the source's **unconditional**
list/pair bounds; do not restrict that cost average to good primes.

Finally, the high/low split concerns a deterministic expectation for each fixed
transformed input. In the high case, failure is contained in preliminary
continuation. In the low case, draw the independent main randomness regardless
of whether the preliminary procedure finishes: a bad final result is contained
in the unconditional bad-main event. Conditioning the main analysis on
`CONTINUE`, or treating continuation as a low-deficiency certificate, is
unnecessary and unjustified.

## 5. An exact continuation counterexample

Take 16 positive original weights: ten copies of `1` and the six weights
`17*2^j`, `j=0,...,5`. Use the legal isolation choice `N0=256` and all tie
weights equal to `1`, so `c_i=4097 a_i+1`. A half with `k` large weights has
exactly `2^k(9-k)` distinct sums: different large subsets have different
integer parts in units of `17`, and the small-weight multiplicity ranges
independently from zero to `8-k`.

For `p=2/5`, the cap is `128`. It overflows precisely when one half contains
all six large weights, giving uniform-half probability
`2 binom(10,2)/binom(16,8)=1/143`. However, with `E k=3`, concavity of `log`
gives `E D(V)>=5-log2(6)>12/5=(3/20)*16`, since `6^5<2^13`.
Thus the input is in the high case for `d=3/20` and has a positive
preliminary gap `p-1/2+d=1/20`, yet it can continue. The probability `1/143`
is conditional on a successful uniform-half draw; sampler failure adds
another possible continuation.

This is not solely a small-size phenomenon. For even `n`, use `n/2` copies
of `1` and `n/2` weights `(n+1)2^j`, with all tie weights one. The same support
formula gives `E D(V)>=n/4-log2(n/4+1)`. Any fixed `d<1/4` is therefore in
the high case for sufficiently large `n`. Yet a half containing all large
weights exceeds any fixed `p<1/2` cap eventually, with probability at least
`2/binom(n,n/2)>0`. Original input bit length is `O(n)`, so the source's
exponential bit-length guard also holds eventually. At `n=128` the replay
already checks the mean/cap comparison for all eight retained `p,d` pairs;
that finite illustration is of the preliminary interface, not the guard.

## Evidence and replay

[audit.py](audit.py) uses only standard-library integers and rational arithmetic.
[AUDIT.json](AUDIT.json) binds its code, parameter snapshot, and source manifest,
and records:

- 462 integer-weight multisets of sizes zero through six, with weights in
  `{-2,-1,0,1,2}`; 18,943 supports, 190,474 disjoint unions, 53,125 increments,
  124,770 exchanges, and 18,943 sign translations.
- All 300 toy `(n,m,alias)` sampling cases through `n=6`, enumerating every
  block and flip; 501 exact variance controls through `n=10`.
- The legal positive-weight continuation example, eight exact parameter tuples,
  component cutoffs, and rational probability bounds at three sizes.

For the proposed degree-fourteen tuple, the limiting Markov lower bound is
only about `0.002519`. Equation (4) instead guarantees at least `0.997011325904`
at `n=10^12`, and equation (5) guarantees at least one half for every
`n>=12,347,133,188`. The cardinality-loss cutoff from the dyadic construction
is `18,153,213`. These are conservative cutoffs for individual proof steps;
they are not practical solver sizes or complete algorithm cutoffs.

From the repository root:

```sh
mkdir -p .context
python3 -B research/deficiency-interface/audit.py --output .context/deficiency-audit.json
python3 -B research/deficiency-interface/test_package.py
make check
git diff --check
```

The controls verify exact finite statements and reproduce the receipt, not
the written argument for every `n`. The source's sparse-prime estimate,
smaller-prime theorem, capped cost bounds, checksum analysis, finite-circle
limits, and final assembly still require their own review. Cached manuscript
sections were checked against all seven hashes in the frozen source manifest;
the source itself is not duplicated in this package.

## Next question

Can the equality boundary `J rho delta=d` be treated with a sublinear allowance
in the useful-block threshold, absorbed by the strict representation gap?
That would change a hypothesis of the current model and needs its own argument
and verification path; no boundary relaxation is claimed here.
