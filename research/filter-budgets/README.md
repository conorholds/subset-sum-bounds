# Subset Sum: an initial filter and parameter exploration

This pilot found two composable improvements in a numerical budget model of
OpenAI's Subset Sum construction: retuning its representation parameters, then
adding four cosine terms to each filtering polynomial. The retained witnesses
pass outward-rounded interval checks for the continuous filtering moments and
the sufficient inequalities below. They are **conditional parameter witnesses**.
We have not independently proved every algorithmic lemma, implemented the full
Subset Sum algorithm, established a new running-time record, or measured practical
Subset Sum performance.

The source is [Subset Sum in Time 2^{0.49n}](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026),
pinned to `openai/math@1581d4cbdd69a0faf34b2089614fe2a814316a6a`.
Its headline bound is randomized worst-case time `O_c(2^{0.49n})` for `n`
integers with at most `n^c` bits each, with success probability at least `2/3`
on every input. The proof's guarded implementation bound is
`2^{0.489995n} poly(n+b)`; its preliminary cap is `2^{0.489985n}` and main work
parameter is `tau = 0.48996`. These are different budgets.
[SOURCE.json](SOURCE.json) records the pinned source sections, hashes, and
attribution. The original algorithm, architecture, and degree-six coefficients
are OpenAI's.

## Retained results

Every decimal in [witnesses.json](witnesses.json) is an exact rational, supplied
as a string. The table rounds interval upper endpoints upward. `alpha` is the
outer exponent budget tested against all the listed components; its algorithmic
interpretation depends on the parameterized proof argument below.

| Witness | Filter degree | Main work `tau` | Outer budget `alpha` | Upper bound on log2 moment 1 | Upper bound on log2 moment 2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `source` | 6 | 0.48996 | 0.49000 | 0.514974050 | 0.514851452 |
| `parameters_only` | 6 | 0.48980 | 0.48981 | 0.517524050 | 0.517572506 |
| `filter_only` | 10 | 0.48996 | 0.49000 | 0.513650736 | 0.513408256 |
| `combined` | 10 | 0.48975 | 0.48976 | 0.517815736 | 0.517819693 |

Changing only the filters lowers both moments at the original parameters. It
does not automatically lower the whole algorithm's exponent. Retuning the
original construction yields the `0.48981` budget; retuning with the new filters
yields `0.48976`. The latter uses 16 representation blocks instead of 17 and
changes the modulus and representation densities. The final gain therefore
survives the other modeled costs, rather than just improving a local objective.

For `combined`, both continuous log2 moments are below `chi = 0.5179` and
`eta - chi - 48 epsilon = 0.00009952 > 0`. Its filtered-list and visited-pair
exponents are below `tau` by more than `0.00001573` and `0.00000919`,
respectively. All retained witnesses have positive coverage, representation,
preliminary, and assembly margins. [CERTIFICATE.json](CERTIFICATE.json) records
every enclosure and binds the witnesses, checker, source manifest, and runtime.

The initial search examined degree-six and degree-ten filters, then alternated
three rounds of simplex refinement and parameter screening. Degree-six tuning
made little progress; degree ten supplied the useful additional gain. Historical
floating-point outputs are retained in [discovery.json](discovery.json), including
unhelpful trials. Subsequent 4,001-evaluation mutation probes starting from the
combined design improved its worst log2 moment by only about `0.00000150`
(degree 10) and `0.00000472` (degree 14). Those probes are exploratory and do not
establish another outer exponent. They suggest that easy local tuning is already
slowing down. No search here establishes an optimum, worldwide novelty, or an
advantage over human-designed numerical optimization.

## Replay and explore

Use a Python environment with [requirements.txt](requirements.txt). The recorded
run uses Python 3.14 and `mpmath 1.3.0`; NumPy is used for exploration and the
quadrature controls. From the repository root:

```sh
python3 -B research/filter-budgets/test_package.py
python3 -B research/filter-budgets/verify.py --output .context/subset-sum-certificate.json
python3 -B research/filter-budgets/search.py --start combined --parameters combined --degree 14 --steps 4000 --output .context/subset-sum-candidate.json
```

`verify.py` uses direct interval sums, not the FFT calculations used by the
explorer. Its default grid has 128 nodes and 50 decimal interval digits. A failed
inequality raises an error; the output labels its limited proof scope.
`test_package.py` checks published scalar nodes, grid refinement, the filter
ablation, deterministic bounded search, and rejection of invalid budgets and a
uniform filter that fails the moment bound.

`search.py` is a bounded mutation explorer for fresh attempts; the retained
witnesses can be checked independently of the historical optimizers. It starts
with an existing design, optionally pads it with new cosine coefficients, and
minimizes the larger log2 moment at fixed parameters. Its output remains
exploratory. To check a candidate, place its `design` and explicit `parameters`
in the witness schema with a unique `id`, then pass that file to
`verify.py --witnesses ...`. Reducing a moment alone does not justify changing
`alpha`; every budget and the new moments must pass again.

## The conditional budget model

The following is a rederivation of sufficient inequalities from the pinned
manuscript's `setup.tex`, `aliases.tex`, `checksum.tex`, and
`implementation.tex`. Checking this finite list does not prove that no
algorithmic requirement was omitted. The remaining dependency is the validity
of the manuscript's qualitative constructions and their extension from fixed
displayed constants to the parameters described here.

Let `eta` be the large modulus bit density, `tau` the main cap, `rho = m/n`
the representation-block density, `J` the number of disjoint `m/2` blocks in a
random half, `d` the half-deficiency threshold divided by `n`, `delta` the
local deficiency threshold divided by `m`, `gamma` the smaller modulus density,
`zeta` the bit-length guard (`b <= 2^{zeta n}`), and `p` the preliminary
cap density. The JSON names are `blocks`, `deficiency`, `local_deficiency`,
`guard`, and `preliminary`, respectively. Let `h` denote binary entropy.

The checker requires `0 < tau < eta < 2 tau`, `0 < rho < 1/2`, positive
`d, delta, gamma, zeta, epsilon`, and the following strict inequalities:

| Inequality | Reason in the source construction |
| --- | --- |
| `p > 1/2 - d` | Chebyshev's inequality leaves a positive gap for both distinct-sum preliminary lists and integer rounding. |
| `J rho < 1` | A random half contains `J` disjoint blocks of size `m/2` for sufficiently large `n`. |
| `J rho delta > d` | Superadditivity of deficiency and Markov's inequality give a constant positive probability of a useful block. |
| `1/2 - delta > gamma` | The useful block has more than `2^{gamma m}` distinct mask weights after the polynomial loss for choosing a cardinality. |
| `h(rho) + rho + zeta < eta` | A union bound preserves every sparse integer distinction modulo the large prime. |
| `rho + 1/4 < tau` | Filtered-list generation setup fits the main work cap. |
| `1/2 + (h(1/4) - 1/2 - gamma) rho < tau` | Expected filtered lists fit their cap with exponentially small cap probability. |
| `1 + (2h(1/4) - 1 - gamma) rho - eta + zeta < tau` | The exact-integer exclusion and random-prime collision bound keep visited pairs below their cap. |
| `chi < eta - 48 epsilon` | Both finite-circle moment bounds give exponentially decreasing checksum variance. |
| `alpha > p`, `alpha > tau + 20 epsilon`, `alpha > eta/2` | All preliminary, checksum, alias, and prime-testing costs fit the outer budget. |
| `alpha < 1/2` | The tested budget improves on the meet-in-the-middle exponent. |

In the preliminary proof choose a deficiency threshold strictly between
`1/2 - p` and `d`; variance `O(n)` then gives failure probability `O(1/n)`.
For coverage, replace the source's 17 blocks by `J`. The limiting useful-block
probability is at least `1 - d/(J rho delta) > 0`. A smaller positive constant
and the same polynomial number of independent trials suffice asymptotically.
The source's list and pair estimates directly give the two entropy expressions
in the table. Their strict gaps absorb floors and polynomially many trials.

For checksums put `Delta = eta - tau`, `Gamma = 2 tau - eta`,
`a = Delta/Gamma`. Keep the source's domination slopes `8.5` and `8.3`.
The checker also requires `8.3(4a+3) < 30`; the other rounding loss is
`8.5 * 3 = 25.5 < 30`. With both finite moments below `2^chi`, the source's
unconditional variance argument has exponent
`-tau - Delta + chi + 48 epsilon = -eta + chi + 48 epsilon < 0`.
No conditioning on a successful extraction or a preliminary continuation is
introduced in this calculation.

The source proves its entropy ordering for any positive normalized densities,
so replacing the fixed-degree polynomials by other fixed finite polynomials and
changing the positive power `kappa` preserves that part of the argument. Their
derivatives and positive denominator bounds remain effective constants.
Strict continuous moment slack implies the required finite-circle bounds for
all sufficiently large moduli by uniform convergence, with an effective fixed
cutoff. The checker certifies the continuous slack; it does not compute that
algorithmic cutoff.

If the inherited qualitative arguments and this parameterization are correct,
all remaining factors are polynomial or `2^{O(sqrt(n))}`. Strict assembly gaps
allow them to be absorbed below `2^{alpha n}` for each fixed polynomial input
bit bound. Failed guards and small `n` use the source's exact fallback. This is
an asymptotic interpretation, with no useful finite-size crossover established.

## Moment certificate

All integrals below are normalized circle integrals and logarithms have base 2.
For cosine polynomials `r,s`, set `F=2^r`, `G=2^s`, `Z_f=int F`,
`Z_g=int G`, and

```text
P(x) = int F(x+w) G(w) dw
Y(x) = int F(x+w) G(w) (r(x+w)+s(w)) dw
L = log2 P - log2(Z_f Z_g)
B = Y/P - log2(Z_f Z_g)
H = log2(int F^kappa) + log2(int G^kappa) - kappa log2(Z_f Z_g)
u = L + (kappa-1) B - H
v = (2-kappa) B - 2 L + H
I1 = int (2+2cos x) 2^{8.5(Delta-u(x))} dx
I2 = int (2+2cos x) 2^{8.3(a v(x)-u(x))} dx
```

These are the source appendix's score identities. [certify.py](certify.py)
encloses their continuous values with `mpmath.iv`; this is a computational
interval certificate with that arithmetic library in its trusted base, rather
than a proof-assistant formalization. Coefficients, budgets, widths, and slopes
enter as exact decimals. Endpoints are serialized outward using exact binary
fractions and directed decimal rounding.

For a periodic analytic function bounded by `B0` on a strip of half-width
`sigma`, its `N`-point trapezoidal error is at most
`2 B0 / (exp(N sigma)-1)`. This follows by bounding each Fourier coefficient
by `B0 exp(-|k| sigma)` and summing the nonzero multiples of `N`.
Writing `R_sigma = sum |r_j| cosh(j sigma)` and likewise `S_sigma`, the
inner bounds are `2^{R_sigma}` for `F`, `2^{kappa R_sigma}` for `F^kappa`,
`2^{R_sigma+S_sigma}` for `P`, and
`(R_sigma+S_sigma) 2^{R_sigma+S_sigma}` for `Y`. Corresponding bounds for
`G` are identical. Direct interval quadratures are widened by these errors
before divisions and logarithms. Positive denominator lower bounds are checked.

The outer integrands contain `log Q`, where `Q=P/(Z_f Z_g)`, so a zero-free
complex strip must be justified. Set `beta=0.12` and `sigma=0.3`. Define
`r_beta(x)=sum r_j cosh(j beta) cos(jx)` and compute upper bounds

```text
U_f >= int |f(w+i beta)|^2 dw = int 2^{2 r_beta(w)} dw / Z_f^2
U_g >= int |g(w+i beta)|^2 dw
theta = sqrt((U_f-1)(U_g-1)) < 1.
```

These integrals get their own inner error bounds. Parseval's identity makes
the shifted squared norms nondecreasing in absolute shift, and shifted means
remain one. Splitting the contour at `z=x+iy` as
`Q(z)=int f(x+w+iy/2) g(w-iy/2) dw` therefore gives
`|Q(z)-1| <= theta` throughout `|Im z| <= 2 beta`. Thus `L` is analytic on
a neighborhood of the closed strip and `-Re L <= -log2(1-theta)`.
Cauchy--Schwarz on the same split contour bounds

```text
|B| <= sqrt(U_f U_g)
       * (R_beta + S_beta + |log2(Z_f Z_g)|) / (1-theta).
```

Substituting these bounds into `Delta-u` and
`a v-u=-(1+2a)L+[a(2-kappa)-(kappa-1)]B+(1+a)H` bounds the absolute
values of the two analytic outer integrands; their prefactor is at most
`2+2cosh(2 beta)`. Apply the trapezoid error bound again with half-width
`2 beta`. This accounts for both quadratures and all nonlinear operations,
not just floating-point agreement at several grid sizes. The retained
outer error bounds are at most about `2.01e-9` in moment units.

## Candidate assessment

This is worth a small collaborative pilot: the first gains are inexpensive to
check, different contribution types compose, and every claimed improvement has
a fixed input artifact. The evidence supports local headroom. It does not yet
show a long sequence of substantial improvements; the simple follow-up searches
already move slowly. Useful next work is an independent audit of the
parameter-to-algorithm argument, testing qualitatively different filter families,
and finding a sharper representation or checksum inequality that changes the
frontier rather than only its last digits.
