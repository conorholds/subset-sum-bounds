# Fixed-filter block-count frontier

Retuning the existing degree-ten filters with **17 blocks** gives the conditional
outer budget **`alpha = 0.4897358`**, improving the predecessor's `0.48976` by
`0.0000242`. A simpler rounded witness gives `0.48974`. The filter coefficients
and `kappa` are unchanged. All three new witnesses pass the existing continuous
interval certificate and all 14 budget constraints.

This investigation also proves that **18 or more blocks cannot beat `0.49` in
the current budget model, for any filter**. The proof uses only preliminary
coverage, block space, representation coverage, and the filtered-list cost.

The predecessor is
[`f25a303`](https://github.com/conorholds/subset-sum-bounds/tree/f25a303525df636c079e6a7c574527cccf56e611).
The investigation is coordinated in
[Discussion 6](https://github.com/conorholds/subset-sum-bounds/discussions/6#discussioncomment-18948726).
The full parameter-to-algorithm argument remains conditional and requires
independent mathematical review.

## Results and complete costs

| Witness | Blocks | Main cap `tau` | Outer budget `alpha` |
| --- | ---: | ---: | ---: |
| Predecessor `combined` | 16 | 0.48975 | 0.48976 |
| `block16_retuned` control | 16 | 0.4897470 | 0.4897570 |
| `block17_retuned` | 17 | 0.4897258 | **0.4897358** |
| `block17_rounded` | 17 | 0.48973 | 0.48974 |

The control and improvement use the same reserve policy. Refining the 16-block
branch alone gives a much smaller gain. Searching the 17-block branch separately
supplies the additional progress; the pilot's local grid refinement stays on
the block count selected by its coarse search.

[witnesses.json](witnesses.json) gives exact decimal strings for every coefficient
and parameter. The best witness has `eta = 0.5179812`, `rho = 0.0577528`,
`d = 0.0102942`, `delta = 0.010511536`, `gamma = 0.489458464`,
`p = 0.4897158`, and `chi = 0.5178812`.

Its certified log2 moment upper bounds are `0.517861635829880532` and
`0.517879496266824808`, both below `chi`. The variance reserve is
`eta - chi - 48 epsilon = 0.00009952`, equal to the predecessor's. Both use
`epsilon = 0.00000001` and the bit-length guard `zeta = 0.00001`.

| Component | Certified upper bound or exact exponent |
| --- | ---: |
| Filtered lists | < 0.489709587 |
| Visited pairs | < 0.489715570 |
| List-generation setup | 0.3077528 |
| Preliminary search | 0.4897158 |
| Main work plus `20 epsilon` | 0.4897260 |
| Prime testing | 0.2589906 |

The preliminary and representation gaps, `alpha-p`, and
`alpha-(tau+20 epsilon)` are unchanged. Coverage, filtered-list, and
visited-pair reserves are larger than the predecessor's.
The new block-space margin is `1 - J rho = 0.0182024`, compared with `0.07496`
before; the larger block count uses more of that space. Its coverage margin
is `J rho delta - d = 0.0000260008171136 > 0`, and its representation gap is
`1/2 - delta - gamma = 0.00003`.

[CERTIFICATE.json](CERTIFICATE.json) records every margin and binds the new
witnesses to the unchanged shared checker and
[source manifest](../filter-budgets/SOURCE.json). Polynomial and subexponential
assembly costs still depend on the
[conditional model](../filter-budgets/README.md#the-conditional-budget-model).

## A necessary block-count bound

Write `H = h(1/4)` and `c = 1-H > 0`. Any witness admitted by the current model
must satisfy:

```text
alpha > p > 1/2 - d
J rho < 1
J rho delta > d
gamma < 1/2 - delta
alpha > tau > 1/2 + (H - 1/2 - gamma) rho.
```

Consequently `d > 1/2-alpha` and

```text
alpha > 1/2 - c rho + delta rho
      > 1/2 - c rho + (1/2-alpha)/J.
```

Rearranging, then using `rho < 1/J`, gives the strict necessary inequality

```text
alpha > 1/2 - J c rho/(J+1) > 1/2 - c/(J+1).
```

The final expression increases with `J`. At `J=18` it exceeds
`0.49006726970837541388`, so every `J >= 18` is excluded when searching for
`alpha < 0.49`. The simpler exact integer inequality `3^75 < 2^119` proves
`H = 2-(3/4)log2(3) > 0.81`, which already gives the cutoff without computing
any logarithm. Also `3^3 > 2^4` proves `H < 1`.

[block_bound.py](block_bound.py) independently encloses the sharper number with
standard-library rational arithmetic. For `1 < x <= 2`, put `t=(x-1)/(x+1)`;
the positive series has bounds

```text
S_m = 2 sum_{k=0}^{m-1} t^(2k+1)/(2k+1)
S_m < ln(x) < S_m + 2 t^(2m+1)/((2m+1)(1-t^2)).
```

Twenty-four terms for `ln(3/2)` and `ln(2)` enclose `log2(3)` and hence `H`.
[BLOCK-BOUND.json](BLOCK-BOUND.json) retains the outward-rounded decimal bounds
and the proof/checker hashes. This is a bound on this budget parameterization;
other representation schemes may obey different inequalities.

## Eliminating the auxiliary budgets

The fixed search reserves are `a=0.00002` for deficiency, `C=0.000026` for
coverage, `R=0.00003` for representation slack, `ell=0.000016` for filtered
lists, and `s=0.00001` for visited pairs. Set

```text
d = 1/2 - tau + a
delta = (d+C)/(J rho)
gamma = 1/2 - delta - R
p = tau - 0.00001
alpha = tau + 0.00001
chi = eta - 0.0001.
```

With `D=1/2+a+C`, the two list/pair constraints become

```text
tau >= [1/2 + (H-1+R) rho + D/J + ell] / (1+1/J)
tau >= [1 + (2H-3/2+R) rho + D/J - eta + zeta + s] / (1+1/J).
```

The first right-hand side decreases with `rho`; the second increases.
Their intersection is

```text
rho = (eta-1/2-zeta-s+ell)/(H-1/2).
```

[frontier.py](frontier.py) tests this balance, clipped to `J rho <= 0.999`,
and screens the remaining constraints. The moment cap depends only on `eta`
for the fixed filter, so its bisections are shared between block-count branches.
The sampled `v` scores must be positive for that bisection; continuous acceptance
is supplied separately by the interval checker.

Each of `J=1,...,17` gets 401 coarse `eta` samples in `[0.510,0.550]` and two
201-point local refinements. Each moment root uses at most 40 bisections.
The retained campaign used 13,651 parameter screens and 219,718 moment
evaluations, taking about 1.65 seconds on the recorded machine. Environment,
limits, inputs, hashes, and every branch's best point are in
[discovery.json](discovery.json). The search is deterministic and bounded;
it does not establish global optimality over all densities, reserve choices,
moduli, or filter families.

## Replay and next question

From the repository root, with the existing requirements installed:

```sh
mkdir -p .context
python3 -B research/block-count-frontier/frontier.py --output .context/block-frontier-discovery.json --witnesses-output .context/block-frontier-witnesses.json
python3 -B research/block-count-frontier/block_bound.py --output .context/block-count-bound.json
make check
```

The search recreates all three retained witnesses. The
[package controls](test_package.py) check the exact bound, admissible-count
boundary, reserve preservation, bounded independent branch refinements, and
opposite density perturbations rejected by the shared budget checker. The full
project check also replays every numerical package. Hosted CI remains manually
triggered by the maintainer.

The next useful filter experiment is to reduce the second moment at the new
17-block parameters while preserving the first: its certified slack below
`chi` is only about `0.00000170`, versus about `0.00001956` for the first moment.
Any filter gain must be followed by the complete budget retuning and replay.
