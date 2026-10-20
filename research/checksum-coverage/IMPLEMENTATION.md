# Checksum coupling, rounding and exponential overhead

The coverage/cap argument in [README.md](README.md) can be composed with
the source's bounded uniform sampler and complex-summation scheme. This
document makes those **checksum-only** error terms explicit. Prime sampling,
alias extraction, preliminary search and final Subset Sum assembly are
outside its scope. The continuous moment ceiling, companion proofs and
the source's deterministic fixed-accuracy score evaluator remain dependencies.

For the retained filter and minimum-block partition:

- At every `n>=2^61`, at most `n^2` targets have checksum error at most `1/8`
  with failure probability at most `2^-128`, under the stated evaluator and
  arithmetic model.
- At every `n>=2^93`, the table-scan and pattern overhead exponent is also
  at most `epsilon n`. This bounds one exponential factor in the operation
  count; it is not a practical crossover or a complete word-RAM cutoff.

All constants are replayed by [implementation_bounds.py](implementation_bounds.py)
and retained in [IMPLEMENTATION.json](IMPLEMENTATION.json).

## 1. Count ordinary uniform draws on every capped history

Write `g=Gbits=99`, `Q<=2^(gq)`, `C=2^c` with
`c=ceil[(tau+20epsilon)n]`, and

```text
a = 2(R+S)+epsilon/2+4epsilon,
sigma = tau-s/4 = 0.48943092,
s = eta-chi-48epsilon.
```

In any pattern the largest block table has at most
`Rmax=2^(a dmax+1)` indexed rows. Since block sizes sum to `n`, there
are at most `n Rmax` coordinate labels per pattern. Generating each label
uses at most one rational Bernoulli draw and one uniform residue draw.
This upper bound permits drawing the residue even if the label is a star.
There are another `n` initial phase draws and at most `T Q 2^ceil(sigma n)`
match draws, even if separate streams are used for every target.

For `T<=n^d`, integer `d>=1`, the total number of ordinary uniform
calls is therefore at most

```text
n+2n Q Rmax+T Q 2^ceil(sigma n)
 <= 7 n^d 2^[gq+max(a dmax,sigma n)]
 <= 2^H,
H = gq+ceil[max(a dmax,sigma n)]+3+d ceil(log2 n).        (9)
```

This holds before imposing a good event. Aborts only reduce the number of
calls. All label probabilities and pattern parameters are fixed rational
numbers; any new draw's modulus may depend on earlier draws.

The source's sampler uses `ceil(log2 m)` fresh bits per proposal below a
positive integer `m`, returning zero without randomness when `m=1`.
Acceptance is at least `1/2`; conditional on acceptance each value below
`m` is equally likely. Capping at `n^3` rejected proposals has failure at
most `2^(-n^3)`. Couple to the ideal unlimited sampler using the same fresh
bits until a cap is reached. Conditional failure at each new call has the
same bound, including adaptively selected moduli. Thus

```text
Pr[checksum bounded and ideal draws differ] <= 2^(H-n^3). (10)
```

No independence between failure events is assumed. Since `q,dmax<=n`,
put `z=ceil(max(a,sigma))=17`. Then
`H<=(g+z+d)n+3`. The exact conditions

```text
N^2>=2(g+z+d), N^3>=2(b+3)
```

imply (10) is at most `2^-b` for every `n>=N`. They hold at `N=2^61`
for `d=2,b=130`. The proposal length and rational Bernoulli denominators
have polynomial bit length, so the rejection cap adds polynomial overhead.
There is no prime-selection term in (10): the checksum accepts an arbitrary
integer modulus satisfying its size condition.
Do not replace this coupling argument by averaging under the law conditioned
on every sampler succeeding: an adaptive call count can bias earlier phases
under that conditional law. The finite controls include this failure mode.

The explicit call count also permits a smaller proposed rejection cap:
**`H+b+ceil(log2 n)` proposals per ordinary uniform call** gives total
coupling failure at most `2^(-b-ceil(log2 n))<=2^-b/n=o(1)`.
The logarithmic reserve preserves vanishing failure, as well as the fixed
finite allowance. This cap is an integer determined before the draws,
and is `O(n)` because `H<=(g+z+d)n+3`. It replaces the source's
`n^3` cap while preserving its coupling allowance and bounded execution.
It does not alter the ideal law or improve the leading algorithm exponent.

## 2. Bound deterministic complex-evaluation and rounding error

All counts, ownership decisions and multiplicities are exact relative to
the fixed approximate-score rule. They are not recomputed at higher
precision. A retained sampled match has `mu>=1`, and the entry cap gives
`L<=C^2`. Consequently its rational weight `L/(P mu)` is at most `C^2`.
The direct weight `1/(P mu)` is at most one. We may use the common bound
`C^2` on every nonaborted history.

Approximate each unit phase to absolute complex error
`delta=2^(-2n^2)`. The Fourier product has `n` factors of norm at most two
and one leading unit phase. Telescoping gives error at most

```text
(n+1)delta(2+delta)^n <= (n+1)2^(-2n^2+2n).
```

Multiply by the exact rational weight, then round each complex contribution
on the grid of step `2^(-n^2)`, with a fixed rule for both coordinates.
Even truncation in both coordinates has complex error at most twice that
step. Each pattern has at most `C` direct terms or sampled terms, because
`sigma<tau+20epsilon`. Division by the common power-of-two sample count
is exact and cannot increase the error. Exact addition aligns the dyadic
denominators; it adds no further rounding. For each target, uniformly over
all nonaborted histories, total arithmetic error is at most

```text
(n+1)2^[-2n^2+2n+3c+gq] + 2^[1-n^2+c+gq].             (11)
```

There is no target-count multiplier in this deterministic per-target bound.
Since `tau+20epsilon<1`, `c<=n`; also `q<=n` and
`ceil(log2(n+1))<=n`. The two exponents in (11) are therefore at most

```text
-2n^2+(g+6)n, and 1-n^2+(g+1)n.
```

For `n>=max(4,2(g+6))=210`, each term is at most `1/32`.
Thus (11) is at most the source's **`1/16` arithmetic allowance**.
The receipt retains the sharper point exponents without materializing
the very small probabilities.

### Linear precision instead of quadratic precision

More generally, let `p` be the unit-phase accuracy in bits and `r` the
fractional bits used to round a weighted contribution. The same argument
gives total error at most

```text
(n+1)2^[-p+2n+3c+gq] + 2^[1-r+c+gq].                  (11a)
```

Choose the proposed deterministic schedules

```text
p = 2n+3c+gq+ceil(log2(n+1))+5,
r = c+gq+6.                                            (11b)
```

Each term in (11a) is at most `1/32`, so the same `1/16` allowance holds
on every nonaborted history. Because `c<=n`, `q<=n`, and
`ceil(log2(n+1))<=n`, these require at most `(g+6)n+5` phase bits and
`(g+1)n+6` contribution bits: **linear rather than quadratic precision**.
At the retained parameters they are
`p=3.46836800 n+99q+O(log n)` and `r=0.48945600 n+99q+O(1)`.
The source uses phase accuracy `2^(-2n^2)` and `n^2` contribution bits.

This change concerns the final complex checksum evaluation. The fixed
score-approximation rule, ownership comparisons and exact weights are
preserved. Averaging and adding dyadics remain exact. Weighted phase
products can still have polynomial-size intermediate rational values;
no specific improvement to the complete polynomial operation count is
claimed. The tests include an exact pessimistic weighted product and a
coarse-precision control whose error exceeds `1/16`.

The phase-accuracy hypothesis has an explicit rational realization.
For requested error `2^-p`, use `p+8` terms of the two alternating series
in `pi/4=arctan(1/2)+arctan(1/3)`. Their combined pi error is at most
`2^(-2p-14)`. Reduce the phase residue first and form its rational angle
using this pi approximation. The angle error is at most `2^(-2p-13)`
and the approximate angle lies in `[0,8]`.
For `m=2(p+8)+256`, the degree `m-1` complex exponential Taylor
polynomial has error at most `8^m/m!<=2^(-m/2)<=2^(-p-136)`:
use `m!>=(m/2)^(m/2)` and `m>=256`.
Round its coordinates to `p+8` fractional bits. The angle, Taylor and
rounding errors sum to less than `2^-p`. All operations are rational with
polynomial bit complexity. Set `p=2n^2` for (11), or use (11b).
The repository does not implement this full evaluator; the exact finite
tests exercise the error inequalities rather than execute it at the cutoff.

## 3. Account for the table and pattern exponential overhead

The table scans used in digit tests and multiplicity retests have exponent
at most `ceil(log2 q)+1+a dmax`. Including enumeration of patterns, put

```text
W(n)=gq+ceil(log2 q)+1+a dmax.
```

The minimum-block partition has `q<=2sqrt(n)`,
`dmax<=max(ceil(sqrt(n)),2t(n))`, and

```text
t(n) <= D0+[log2 n+g+9]/(2epsilon)+2.
```

Consequently

```text
W(n) <= A sqrt(n)+B log2(n)+K,
A=2g+a, B=1+a/epsilon,
K=2+5a+2a D0+a(g+9)/epsilon.                            (12)
```

For `n>=4`, each of `sqrt(n)/n`, `log2(n)/n` and `1/n` decreases.
At a dyadic `N=2^k`, the sufficient rational check

```text
A ceil(sqrt(N))+B k+K <= epsilon N
```

therefore implies `W(n)<=epsilon n` for every integer `n>=N`.
The first dyadic power satisfying this particular envelope is **93**.
At `2^61`, the table-scan bound alone exceeds `epsilon n`; the smaller
probability threshold cannot be presented as a matching runtime threshold.

Combining `W(n)<=epsilon n` with the source's cap exponent bounds the
checksum work's exponential factor by `2^[(tau+21epsilon)n]`, with
polynomial factors still to be quantified. The fixed filter evaluator,
arithmetic, sorting and storage costs are not all assigned explicit
constants here. This is not a complete finite-size machine resource bound.

## 4. Compose failure and accuracy

At `N=2^61`, rerun the coverage/cap and sampling conditions with random
failure allowance `2^-130` each, and allow `T<=n^2`. All uniform conditions
still pass. Equation (10) supplies a third failure at most `2^-130`.
On the remaining event, ideal sampling error is at most `1/16` and (11)
adds at most `1/16`, giving error at most `1/8` for all targets.
The total random failure is at most

```text
3*2^-130 < 2^-128.
```

At `N=2^93` the same error statement holds and (12) also pays for table
and pattern overhead. No independence between the three failures is
required. The continuous filtering certificate and written precursor
proofs remain conditional inputs; these calculations are not independent
scientific review or a new complete Subset Sum running-time record.
The proposed linear rejection and precision schedules preserve these
same failure and accuracy allowances.
