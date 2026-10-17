# Asymmetric Fourier filters

Allowing sine terms and entering a different filter basin gives the
interval-checked conditional budget **`alpha = 0.48946580`**. The preceding
fractional-coverage candidate was `0.48972043`; the preceding
tunable-multiplier candidate was `0.48972642`.

The retained construction has degree twenty, including both cosine and sine
coefficients. Its asymmetry cannot be removed by a common rotation. The
source's normalized-power tests and separable product kernel are retained.
The fractional coverage argument is essential: `rho=0.05929405` exceeds
`1/17`, and the old integer-packing coverage test fails for every fitting `J`.

The result is a numerical proposal within the explicitly generalized model,
with a proposed extension of its analytic bounds. Independent mathematical
review of the model-to-algorithm implication remains pending. No selected
result is promoted here.

## Source and exact predecessor

The scientific source remains
[`openai/math@1581d4c`](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026).
`filtering.tex` starts with arbitrary strictly positive continuous densities
and derives the entropy ordering for their kernel. Its later explicit
numerical example uses cosine polynomials. The change here is to that
numerical family; all densities remain positive smooth functions fixed
independently of the input.

[predecessor.json](predecessor.json) retains the exact preceding
[fractional-coverage](../fractional-coverage/README.md) candidate and its full
source-file hash, canonical project main `f25a303`, and the manuscript pin.
The source's algorithm and filtering architecture are OpenAI's. The budget
reserves, derivative approach, and multiplier tuning build on PRs
[12](https://github.com/conorholds/subset-sum-bounds/pull/12),
[13](https://github.com/conorholds/subset-sum-bounds/pull/13), and the retained
[tunable-multiplier package](../chernoff-multipliers/README.md).

## 1. Construction and the changed analytic bound

Use exact rational coefficients in

```text
r(x)=sum_j [r_j cos(jx)+r_sin_j sin(jx)],
s(x)=sum_j [s_j cos(jx)+s_sin_j sin(jx)].
```

As before, normalize `F=2^r`, `G=2^s` to densities `f,g`, use
`K_x(w)=f(x+w)g(w)/Q(x)`, and score against the normalized `kappa` powers
of `f,g`. The score identities, relative-entropy ordering, table coverage,
list balancing, and checksum domination use positive probability densities;
they require no reflection symmetry. The finite-circle convergence and fixed
accuracy computation arguments extend to fixed real trigonometric polynomials
with sine terms, whose values and derivatives are uniformly bounded.

The new design schema explicitly includes all five fields
`r,s,r_sin,s_sin,kappa`. Its [checker](interval_check.py) retains the direct
inner/outer interval quadratures and multiplier-specific remainder bounds,
with these changes:

1. Evaluate the complete sine/cosine polynomials on the full circular grid.
2. Replace each strip norm by
   `sum_j sqrt(c_j^2+d_j^2) cosh(j sigma)`.
3. For shifted squared norms, double **both** coefficients by `cosh(j beta)`.

For the strip norm, write one harmonic as
`sqrt(c^2+d^2) cos(jz-phi)`. Its magnitude on `|Im z|<=sigma` is at most
`sqrt(c^2+d^2) cosh(j sigma)`. For the shifted norm,

```text
|2^r(x+i beta)|^2
 = 2^(sum_j 2 cosh(j beta)[r_j cos(jx)+r_sin_j sin(jx)]).
```

Thus the same interval trapezoid error applies to the normalizers and these
shifted norms. This also controls numerator errors involving `r+s`.

The zero-free disk argument does not assume even functions. Split a complex
difference `x+iy` between `f(x+w+iy/2)` and `g(w-iy/2)` by contour translation.
Both shifted means remain one. Cauchy--Schwarz applied after subtracting one
gives `|Q(x+iy)-1|<=sqrt((U_f-1)(U_g-1))` for `|y|<=2 beta`.
For real-valued densities, opposite Fourier coefficients are conjugates;
Parseval expresses each shifted squared norm as a sum with `cosh(2j beta)`.
It is therefore even and nondecreasing with the strip width. Neither step
requires reflection symmetry of `f` or `g`.

Once that disk has radius below one, the logarithm bound and the analytic
outer trapezoid remainder are unchanged. At the improved candidate the
certified radius is below `0.396777`. The first log moment is enclosed by
`[0.518359984927380674,0.518359987587243486]`, and the second by
`[0.518359932551748532,0.518359933839640089]`; both are below
`chi=0.51836094`.

The checker imports the explicit fractional budget model and binds its code.
It preserves the old checkers and receipts. Its separate model identifier is
`asymmetric-filters-v1`; an old cosine-only receipt cannot certify this model.

## 2. Search and retained controls

A local sine perturbation at an even filter can miss a different basin.
[density_probe.py](density_probe.py) instead explores positive marginal test
densities on 128 circular nodes. It tries six specified initializations:
the predecessor, a uniform fixed-point control, opposed single peaks, split
and single peaks, shifted Poisson kernels, and triangular supports.
It allows 1,200 accepted iterations per initialization, 36 bounded bisections
for the conditional kernel temperature, and 24 line-search trials per step.
The scalar objective combines the two moments with weight `0.4245`.

During this exploratory stage, a difference-dependent temperature selects
the conditional product law. With fixed marginal tests, its first-order
condition is bisected between `(1+a)/(1+2a)` and one. With the conditional
law fixed, a weighted marginal update decreases the scalar objective; each
step is checked by a line search. The uniform initialization remains a fixed
point despite its poor objective, illustrating the method's limits.

[density-discovery.json](density-discovery.json) retains all six summaries,
the selected positive densities and temperatures, code/input hashes, and
runtime versions. These finite-grid values are exploratory. For the final
construction, project the selected **log densities** onto twenty cosine and
sine harmonics, set a single fixed `kappa=1/0.96735` rounded to twelve
decimals, and restore the source's product kernel. The adaptive-temperature
law is not used by the retained algorithmic candidate.
[start.json](start.json) is that exact projection.

[spectral.py](spectral.py) differentiates the full real Fourier representation
through the normalized correlations. [search.py](search.py) optimizes its
two implicit `eta` boundaries on the fractional-coverage cost slice, including
the derivative of that slice. It uses 128 nodes, forty root bisections in
`[.5175,.5195]`, smoothing temperatures `.01,.001,.0001,.00001`, at most
100 iterations per stage, 24 backtracks per step, and 4,000 evaluations.
Coefficients stay within `.35` of the projected start, `kappa` in
`[.95,1.12]`, and multipliers in `[7,9.99]` and `[7,9.15]`.
The retained run uses 1,028 frontier evaluations, including rounded retention.

The interval certificate supplies acceptance; search gradients and finite
grids do not. Coefficients, `kappa`, and multipliers use exact twelve-decimal
strings. The preceding conservative parameter-rounding policy and all its
minimum positive reserves are preserved.

| Candidate | Conditional `alpha` | Purpose |
| --- | ---: | --- |
| `even_predecessor` | `0.48972043` | Embed the preceding design with zero sine coefficients. |
| `asymmetric_filter_only` | `0.48972043` | Improve moments while preserving every preceding budget. |
| `asymmetric_projection_retuned` | `0.48958698` | Retune the projected filters with the preceding multipliers. |
| `asymmetric_degree20_retuned` | `0.48946580` | Optimize all coefficients, `kappa`, multipliers, and the frontier. |

[candidates.json](candidates.json) and [CERTIFICATE.json](CERTIFICATE.json)
contain all exact inputs and generated interval bounds.

### Genuine asymmetry and retained reserves

Let the first two cosine/sine pairs of `r` be `(c1,d1),(c2,d2)`.
A necessary condition for any rotation to make `r` even is

```text
2 c2 c1 d1-d2(c1^2-d1^2)=0.
```

It is the imaginary part of the rotation-invariant expression
`(c2-i d2)(c1+i d1)^2`. At the retained candidate its exact rational value
is nonzero, approximately `-0.23824215`. Therefore the sine terms cannot
all be removed by changing the common origin.

The best tuple retains `epsilon=1e-8`, bit guard `1e-5`, deficiency reserve
`0.00002`, useful-block separation at least `0.0000016` per `n`, representation
reserve `0.00003`, list reserve `0.000016`, visited-pair reserve `0.000010`,
moment reserve `0.0001`, and assembly reserve `0.00001`. Its largest strictly
fitting integer packing count is 16; that old coverage inequality fails.
The imported fractional lemma supplies the replacement argument.

## Replay and limitations

From the repository root, with the existing requirements installed:

```sh
python3 -B research/asymmetric-filters/density_probe.py --author conorholds --output /tmp/asymmetric-density.json
python3 -B research/asymmetric-filters/search.py --author conorholds --output /tmp/asymmetric-joint.json
python3 -B research/asymmetric-filters/interval_check.py --author conorholds --check-receipt research/asymmetric-filters/CERTIFICATE.json
python3 -B research/asymmetric-filters/test_package.py
make check
git diff --check
```

The test router replays this explicit model as part of `make check`.
The legacy-only `make verify` command does not discover it automatically.
Hosted verification remains manual and maintainer-only.

The searches are bounded local experiments, not a global optimum over positive
densities. The exact candidates are reproducible independently of whether a
different numerical environment follows precisely the same search trajectory.
The generalized strip proof, fractional coverage application, and complete
parameter-to-algorithm implication still need independent review. Practical
crossover sizes and an executable full Subset Sum solver are not established.

## Next question

Can a continuously certified difference-dependent kernel improve on the
constant-power projection retained here?
