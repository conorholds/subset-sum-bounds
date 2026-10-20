# Alias costs and the exact-integer exclusion

This audit checks the generalized alias-list and visited-pair cost formulas,
complementing the deficiency/coverage argument in
[PR 14](https://github.com/conorholds/subset-sum-bounds/pull/14). It also proves
a concrete obstruction to delaying the source's exact-integer exclusion:
**equal positive weights can produce exponentially many overlapping exact
pairs, even though there are no modular aliases**.

For the proposed degree-fourteen parameters, the obstruction has expected
uncapped pair exponent at least `1.0076880743844`, up to polynomial factors,
on an infinite sequence of compatible sizes. The correctly excluded nonexact
pair estimate has exponent below `0.489713855`. The difference explains why
the large-prime collision saving cannot be charged to exact pairs.

The contribution is a scoped proof and a finite execution oracle. It changes
no selected exponent, checker, or model hypothesis. The proofs still require
independent review; no full Subset Sum solver or complete algorithmic audit
is supplied.

## Starting point and source obligations

Canonical main is pinned at
[`f25a303`](https://github.com/conorholds/subset-sum-bounds/tree/f25a303525df636c079e6a7c574527cccf56e611).
The work is coordinated in
[Discussion 5](https://github.com/conorholds/subset-sum-bounds/discussions/5#discussioncomment-18974325).
The manuscript pin is `openai/math@1581d4cbdd69a0faf34b2089614fe2a814316a6a`,
as retained in [SOURCE.json](../filter-budgets/SOURCE.json).

| Source | Labels | Obligation checked here |
| --- | --- | --- |
| [aliases.tex](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/build/sections/aliases.tex) | `eq:alias-raw-sizes`, `eq:alias-generation` | Balanced partitions and output-sensitive list setup. |
| Same section | exact-integer exclusion paragraph, `eq:alias-fixed-difference` | Skip forbidden intervals without visiting their records, including overlaps. |
| Same section | `eq:alias-list-bound`, `eq:alias-large-modulus-bound`, `eq:alias-pair-bound`, `eq:alias-cap` | General exponents, probability laws, and all-history caps. |
| [setup.tex](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/build/sections/setup.tex) | `eq:transformed`, `lem:prime-collision` | Legal transformed input and the nonzero integer collision estimate. |

[parameters.json](parameters.json) retains the eight distinct parameter tuples
from main, [PR 12](https://github.com/conorholds/subset-sum-bounds/pull/12), and
[PR 13](https://github.com/conorholds/subset-sum-bounds/pull/13), with exact
source-package hashes and commits. No filter coefficients enter this argument.

## 1. Balanced raw lists and their generation

Fix a block `M` of even size `m=2 floor(rho n/2)` and a layer `0<=k<=m/2`.
The ordinary coordinates split into `L,R` of sizes `l,r`, with `l+r=n-m`.
The raw record counts are

```text
N_L = 2^l binom(m,k)
N_R = 2^r binom(m,m/2-k)
N_L N_R = 2^(n-m) binom(m,k) binom(m,m/2-k).                    (1)
```

Incrementing `l` doubles the left count and halves the right count. Rounding
their real-valued equality point to a nearest integer therefore balances them
within a factor two. The equality point lies in `[(n-2m)/2,n/2]`, since each
binomial coefficient is between one and `2^m`. Thus the partition exists for
`m<n/2`, and `max(l,r)<=n/2+O(1)`. Exact integer comparisons suffice to choose
it; neither the prime nor the weights influence this partition.

To generate a filtered side, split its ordinary positions in half. On one
half attach all permitted block masks, and on the other enumerate only ordinary
subsets. Sort by weight modulo the smaller prime and emit complementary
residue pairs. This takes a polynomial factor times setup plus emitted records,
with setup at most

```text
2^(m+ceil(max(l,r)/2)) <= 2^((rho+1/4)n+O(1)).                  (2)
```

Sorting, residue searches, weights of bit length `O(b+log n)`, and masks of
length `n` all add polynomial factors. A cap can interrupt before the next
filtered record is emitted. The finite oracle below deliberately enumerates
raw records exhaustively; it checks the identities and cap semantics, not the
asymptotic speed of this generation procedure.

## 2. One residue, not two independent filters

Conditional on every choice except the uniform residue `j' mod P'`, a left
record survives for exactly one value of `j'`. A right record likewise survives
for exactly one value, using the complementary residue. Thus each expected
list length is its raw length divided by `P'`.

A pair whose reduced integer weights sum to the required target survives
both filters for exactly **one common residue**. Its survival probability is
`1/P'`, rather than `1/(P')^2`; no independence between the two filters is
available. The exact replay checks this identity by summing over every residue.

Write `H=h(1/4)`. Balancing (1), applying the binomial entropy bound and
`h(k/m)+h(1/2-k/m)<=2H`, and using `1/P'<=2^(1-gamma m)` give

```text
E[one filtered list length]
  <= 2^(n/2+(H-1/2-gamma)m+O(1)).                              (3)
```

The entropy inequality follows from concavity, with equality at `k/m=1/4`.
Together with (2), this gives the model's setup and list conditions
`rho+1/4<tau` and `1/2+(H-1/2-gamma)rho<tau`.

## 3. The nonzero discrepancy and original probability laws

Fix a transformed input, target index, carry index, layer, block, flip, and two
raw records before averaging over the large prime. All of these choices and
the ordinary partition can be fixed independently of that prime. The carry
index ranges over the predetermined set `0,...,n`; it is not selected from a
prime-dependent alias during this cost analysis.

Let `mu_i` be the combined record multiplicity, possibly two on overlapping
masks, and define the **unreduced** integer

```text
D = sum_i (-1)^F_i mu_i c_i - t_s + sum_{i in F} c_i.           (4)
```

This discrepancy is independent of the large prime. The reduced integer
matching test implies `P divides D`, because reduced weights are congruent to
the original weights and the carry target is congruent to `t_s`. The exclusion
requires `D!=0` before any matching pair is visited. For a nonzero `D`, the
source's random-prime lemma gives

```text
Pr_P[P divides D] <= poly(n)(b+log n) 2^(-eta n).               (5)
```

Multiplicities at most two do not change the bit-length bound. For `D=0`,
exclusion makes the contribution exactly zero; equation (5) is inapplicable.
Taking a union of masks before computing (4) would lose these multiplicities
and change the tested integer, so that substitution is invalid for overlaps.

Conditional on the primes and the raw records, the smaller-residue probability
is at most `2^(1-gamma m)` times the indicator of the two integer tests.
Average over the original large-prime law using (5), then sum over (1). This
does not assume independent collision events for the two moduli and does not
condition the large prime on successful coverage. The result is

```text
E[visited pairs]
 <= poly(n)(b+log n) 2^(n+(2H-1-gamma)m-eta n).                 (6)
```

Under the input guard `b<=2^(zeta n)`, the exponent is
`1+(2H-1-gamma)rho-eta+zeta`. Both (3) and (6) hold without low deficiency.
Fixing a prime-dependent alias first would not establish these cost bounds;
the proof fixes raw records and sums their indicators instead.

## 4. Exclusion without visiting, caps, and saved aliases

Group right records by their reduced integer weight. Inside each group sort
by unreduced integer weight, retaining duplicate records. For a left record,
select its complementary reduced group. Every forbidden right record has
the same unreduced weight, `t_s^F-left_integer`. Two binary searches identify
its entire interval, and iteration over the prefix and suffix visits exactly
the permitted pairs. Equal keys and overlapping masks cause no exception.

Check mask disjointness only for visited pairs. On disjoint masks, undoing
the flip gives a genuine subset whose reduced sum matches the target and
whose unreduced sum differs. Direct verification prevents a false alias from
being saved, and deduplication is by subset separately for each target.
An overlapping pair is rejected even if its multiplicity sum is a nonexact
modular match.

The costs are linear in filtered lists plus permitted visits, up to polynomial
factors, regardless of the number of excluded exact pairs. Caps at
`2^ceil(tau n)` stop before an extra list record or visit is stored. If each
exponent in (2), (3), and (6) is strictly below `tau`, Markov and a union bound
over polynomially many targets, carries, layers, and trials give total cap
failure `o(1)`. The all-history cost is then `2^(tau n) poly(n+b)`, including
histories that abort. This still imports the source's bounded random draws
and fixed bit-arithmetic conventions.

## 5. A homogeneous overlapping-pair obstruction

Choose even `n`, all original weights `a_i=1`, original target `t=n/2`, and
all tie weights `u_i=1`. Any legal isolation choice `12n<=N0<24n` gives
`c_i=C+1=:c>0`, where `C=1+n N0`. The transformed target with `s=n/2` is
`c n/2`. These are legal positive inputs with polynomially bounded weights.

For sufficiently large `n`, every prime in the prescribed large-prime interval
exceeds `n c`, since `eta>0`. Thus `b_i=c_i`, the relevant carry is zero,
and there are no modular aliases. Nevertheless there are many raw pairs whose
multiplicity weights equal the target. Correct exclusion visits **none** of
these pairs, for every flip and smaller residue.

Take `m` divisible by 32 and layer `k=m/4`. Retain mask pairs whose overlap
has size `z=m/16`; their four membership classes in `M` have proportions

```text
both: 1/16; left only: 3/16; right only: 3/16; neither: 9/16.
```

The number of such ordered mask pairs is the multinomial type-class size `A`.
Its entropy is `2 h(1/4)=2H`, so

```text
A >= 2^(2H m)/(m+1)^4.                                        (7)
```

An elementary justification avoids asymptotic Stirling errors. In the
multinomial distribution with these four probabilities, the integer mean
type is a mode: for each component count `t`, the sequence `t^a/a!` is
maximized at `a=t` (and `t-1`). Ratios of the product masses therefore show
that this type has maximum probability. At most `(m+1)^4` types partition
probability one, proving (7). Its exact rational form is also replayed.

There are `2^(n-m)` ordinary choices per mask pair. Under a uniform flip,
the contribution of an index with multiplicity zero or one is a uniform
choice in `{0,1}` in units of `c`; for multiplicity two it is a uniform
choice in `{-1,2}`. Every coordinate has mean `1/2`. Since both `z` and
`n-z` are even, the event that each group has half its choices on each side
gives exact target sum and probability at least

```text
[binom(z,z/2)/2^z] [binom(n-z,(n-z)/2)/2^(n-z)]
 >= 1/((z+1)(n-z+1)) >= 1/(n+1)^2.                             (8)
```

The central binomial bound follows because the largest coefficient among
`r+1` coefficients is at least their sum `2^r` divided by `r+1`.
For every such exact pair, both smaller filters hold for one common residue,
with probability `1/P'>2^(-1-gamma m)`, independently of the overlap test.
Equations (7)--(8) imply that **delaying exact exclusion until after visiting**
would yield expected overlapping visits at least

```text
2^(n+(2H-1-gamma)m) / [2(m+1)^4(n+1)^2].                      (9)
```

This lower bound averages over the same original flip and prime laws, not a
rare fixed flip. It applies on an infinite sequence of sizes compatible with
the actual rounding: for each retained rational `rho=a/b`, take `n` divisible
by `32b`, so `m=rho n` is divisible by 32. The expected exponent in (9) is
`1+(2H-1-gamma)rho`, with no large-prime saving. It is above one for all eight
retained tuples; the degree-fourteen tuple gives at least `1.0076880743844`.

The equal-weight inputs have high deficiency, so this is an obstruction to
the **unconditional uncapped cost estimate**, not a counterexample to the
source's correct algorithm or a theorem that every alternative extractor
fails. Caps would still bound work, and the preliminary search usually handles
these inputs. The point is that the displayed cost proof cannot charge exact
overlapping pairs to the nonzero-discrepancy estimate (5).

## Evidence and replay

[oracle.py](oracle.py) implements grouped exclusion and a separate exhaustive
full-mask reference. Its toy domain is `4<=n<=10`, even `2<=m<=n/2`, primes
at most 1009, and caps at most 16,384. It is an exhaustive finite oracle,
not the fast list generator used in the asymptotic argument.

[audit.py](audit.py) and [AUDIT.json](AUDIT.json) retain code/input hashes,
source bindings, and these bounded controls:

- 23 balanced-partition identities and 9,960 layer comparisons across five
  weight families, including negative weights, zero weights, flips, carries,
  both smaller primes, and every residue.
- 352 end-to-end controls recovering exactly the genuine aliases balanced on
  the chosen block when all layers, carries, and residues are tried.
- A legal eight-coordinate homogeneous example (`N0=128`, transformed weight
  `1026`, primes `101` and `3`): across all 768 flip/residue choices,
  exclusion skips 16,128 exact pairs, including 2,688 overlapping ones, and
  visits zero. Delayed exclusion would average `7/2` overlapping visits.
- Exact type-class and flip lower bounds at eight sizes through `n=512`, plus
  interval-enclosed model costs for all eight retained parameter tuples.

The tests also cover list/pair cap failures, skip behavior with a zero visit
cap, one shared residue factor, fixed-discrepancy independence, and explicit
rejection of invalid or larger campaigns. The structural counting controls
use exact integers and fractions. Cost profiles use the unchanged shared
interval arithmetic; none supplies a new filter certificate.

From the repository root, with the existing requirements installed:

```sh
mkdir -p .context
python3 -B research/alias-exclusion/audit.py --author conorholds --output .context/alias-audit.json
python3 -B research/alias-exclusion/test_package.py
make check
git diff --check
```

The finite enumeration does not prove all-size claims. The written argument
still imports the manuscript's prime-density/collision lemma, bounded sampling,
and bit-operation model. Coverage, checksums, finite-circle convergence, and
final assembly are separate obligations. In particular this audit is not an
independent endorsement of the entire selected algorithmic claim.

## Next question

Can the exact-integer exclusion be implemented inside the checksum's implicit
match machinery with the same count and cap guarantees? The finite oracle
here covers the explicit alias lists only; an implicit implementation needs
its own invariants and cost audit.
