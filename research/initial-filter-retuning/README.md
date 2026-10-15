# Initial finding: filter and parameter retuning

The first experiment found two improvements that compose in the
[conditional budget model](../filter-budgets/README.md#the-conditional-budget-model):
retuning the representation parameters, then extending both degree-six cosine
polynomials to degree ten. The resulting outer budget is `alpha = 0.48976`.
Independent review of the parameter-to-algorithm argument is pending.

| Change | Filter degree | Blocks `J` | Main cap `tau` | Outer budget `alpha` |
| --- | ---: | ---: | ---: | ---: |
| Source construction | 6 | 17 | 0.48996 | 0.49000 |
| Parameters only | 6 | 17 | 0.48980 | 0.48981 |
| Filters only | 10 | 17 | 0.48996 | 0.49000 |
| Filters and parameters | 10 | 16 | 0.48975 | **0.48976** |

The source row uses the manuscript's headline exponent; its internal guarded
bound is `2^(0.489995n) poly(n+b)`. Exact coefficients and parameter tuples for
all four rows are in [witnesses.json](../filter-budgets/witnesses.json), with
interval evidence in [CERTIFICATE.json](../filter-budgets/CERTIFICATE.json).
The mathematical source is pinned in [SOURCE.json](../filter-budgets/SOURCE.json).

## What changes

At the source parameters, the new filters lower the certified upper bounds on
the two log2 moments from `0.514974050` and `0.514851452` to `0.513650736` and
`0.513408256`. Every displayed bound here rounds upward. This local gain leaves
room to reduce the main cap, but the representation and list costs must also
fit. The combined witness changes the block count from 17 to 16 and retunes
the modulus and representation densities.

For the combined witness, `eta = 0.51800`, `rho = 0.057815`, `d = 0.01027`,
`delta = 0.01113`, `gamma = 0.48884`, `zeta = 0.00001`, `p = 0.48974`,
`epsilon = 0.00000001`, and `chi = 0.5179`. Its continuous log2 moments are
below `0.517815736` and `0.517819693`, both less than `chi`. Consequently the
modeled checksum-variance exponent is
`chi - eta + 48 epsilon = -0.00009952`.

The other modeled costs fit as follows. These are exponent coefficients;
polynomial and `2^(O(sqrt(n)))` factors still require the source argument and
the strict assembly slack.

| Component | Certified upper bound or exact exponent |
| --- | ---: |
| Filtered lists | < 0.489734261 |
| Visited pairs | < 0.489740805 |
| List-generation setup | 0.307815 |
| Preliminary search | 0.48974 |
| Main checksum work, including rounding allowance | 0.4897502 |
| Prime testing | 0.259 |

Coverage is also strict: `J rho = 0.92504 < 1`,
`J rho delta - d = 0.0000256952 > 0`, and
`1/2 - delta - gamma = 0.00003 > 0`. The shared certificate records all 14
budget margins, including sparse distinctions and preliminary-search coverage.
This is a conditional numerical candidate, not an independently established
Subset Sum running-time theorem or a practical speed measurement.

## Replay the discovery

[replay.py](replay.py) runs the original bounded coordinate search and alternating
simplex/parameter refinement, using the shared numerical filtering functions.
It begins with the source design; the retained combined design is used only
for the final comparison. Historical outputs, including unsuccessful degree-six
trials and later mutation probes, are in
[discovery.json](../filter-budgets/discovery.json).

The seed is `20261015`. Each of degrees six and ten gets at most 65 coordinate
sweeps, then three rounds of at most 3,500 simplex iterations. Each parameter
scan tests block counts 11 through 22, three initial block-fill ratios, and 81
modulus densities, followed by two bounded local grid refinements. Moment
screening uses a 128-node grid, a fixed `0.0001` reserve, and at most 32
bisection steps per screened point. These floating-point screens guide the
search; only the separate interval checker validates a retained witness.

From the repository root, after installing the existing requirements:

```sh
mkdir -p .context
python3 -B research/initial-filter-retuning/replay.py --output .context/initial-retuning-replay.json
make check
```

The replay records its seed, limits, environment, input hashes, intermediate
designs, evaluation counts, and whether its final degree-ten coefficients
round to the retained six-decimal design. Floating-point trajectories can
vary across environments. `--quick` runs a smaller smoke campaign; it does not
reproduce the retained result. Certificate replay checks the exact retained
witnesses independently of either search trajectory.

An observed full replay with Python 3.14.3 and NumPy 2.5.1 took 22.5 seconds,
made 1,922 coordinate-objective evaluations and 56,746 simplex-objective
evaluations, and reproduced the retained degree-ten design after rounding.
Parameter-screen evaluations are additional. This timing describes discovery
on one machine, not Subset Sum performance.

## Questions for review

- Does the generalized budget model account for every use of the source's
  fixed parameters, particularly the change to 16 blocks and its coverage
  probability? See the [proof-audit Discussion](https://github.com/conorholds/subset-sum-bounds/discussions/5).
- Do the unconditional checksum-variance argument, finite-circle cutoff, and
  final operation count remain valid with the changed filters and parameters?
  The interval checker certifies continuous moments and numerical inequalities;
  it does not settle those implications.
- Which constraint limits the next improvement: the filter moments, the
  representation costs, or the block-count boundary? Coordinate parameter work
  in [Discussion 6](https://github.com/conorholds/subset-sum-bounds/discussions/6)
  and changed filter methods in [Discussion 7](https://github.com/conorholds/subset-sum-bounds/discussions/7).
