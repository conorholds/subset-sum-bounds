# Analytic-gradient minimax filters

Changing the filter-search method gives a degree-fourteen witness with
**conditional outer budget `alpha = 0.48973396`**. This improves the proposed
17-block predecessor's `0.4897358` by `0.00000184`; canonical main's selected
candidate remains `0.48976`. The new coefficients, continuous moments, and all
14 budget constraints pass the unchanged shared interval checker.

This is a numerical contribution inside the existing conditional model, not
an independently proved Subset Sum running-time record. Neither predecessor
promotion nor the full parameter-to-algorithm proof is supplied here.

## Starting point and changed mechanism

The checkout starts at canonical main
[`f25a303`](https://github.com/conorholds/subset-sum-bounds/tree/f25a303525df636c079e6a7c574527cccf56e611).
The comparison is the exact `block17_retuned` witness from
[PR 12](https://github.com/conorholds/subset-sum-bounds/pull/12), pinned at
[`1767944`](https://github.com/conorholds/subset-sum-bounds/tree/176794411d1facdc299e1f4dae99201f851e2d4f);
[baseline.json](baseline.json) retains that witness and its provenance. The
mathematical source and shared checker are still pinned by
[SOURCE.json](../filter-budgets/SOURCE.json).

This investigation is coordinated in
[Discussion 7](https://github.com/conorholds/subset-sum-bounds/discussions/7#discussioncomment-18950342).
It changes the filter coefficients and optimization method, complementing the
predecessor's fixed-filter block-count analysis. It does not repeat that search.

The existing explorer accepts random coefficient mutations that reduce the
larger moment. Here, [optimize.py](optimize.py) differentiates both moments
through normalization and convolution, then minimizes a smooth approximation
to their maximum by deterministic BFGS descent with backtracking. Four decreasing
smoothing temperatures avoid requiring each trial to improve the nonsmooth
objective. The best actual maximum seen across all trials is retained, and
the exact nine-decimal design is evaluated again after rounding.

## Results and controls

Both methods start from the same degree-ten design, extended by zero
coefficients for degree fourteen, at the predecessor's fixed
`eta = 0.5179812`, `tau = 0.4897258`.

| Filter search | Degree | Evaluations | Worse sampled log2 moment |
| --- | ---: | ---: | ---: |
| Starting filter | 10 | 1 | 0.517879496 |
| Existing mutation control | 10 | 4,001 | 0.517873172 |
| Analytic-gradient descent | 10 | 450 | 0.517868878 |
| Existing mutation control | 14 | 4,001 | 0.517865279 |
| Analytic-gradient descent | 14 | 353 | 0.517839486 |

The displayed values round upward. The degree-fourteen descent lowers the
sampled maximum by about `0.00004001`, compared with `0.00001422` for the
fixed-seed mutation control. A gradient evaluation computes both values and
derivatives, so evaluation counts are not equal units of computational cost.
This is a comparison of two bounded campaigns, not evidence that one optimizer
wins for every start, seed, family, or machine.

The degree-ten result slightly worsens the first moment while improving the
maximum. The degree-fourteen filter improves both moments at the fixed
parameters. Filter-only ablations keep every predecessor budget unchanged;
retuned ablations restore all other costs and positive reserves.

| Witness | Outer budget `alpha` | Certified upper log2 moment 1 | Certified upper log2 moment 2 |
| --- | ---: | ---: | ---: |
| `predecessor_block17` | 0.4897358 | 0.517861636 | 0.517879497 |
| `degree10_filter_only` | 0.4897358 | 0.517868849 | 0.517868879 |
| `degree10_retuned` | 0.48973525 | 0.517880154 | 0.517880861 |
| `degree14_filter_only` | 0.4897358 | 0.517839337 | 0.517839487 |
| `degree14_retuned` | **0.48973396** | 0.517880647 | 0.517883184 |

[witnesses.json](witnesses.json) contains exact decimal coefficients and every
budget parameter. [CERTIFICATE.json](CERTIFICATE.json) supplies continuous
interval enclosures, source/checker hashes, and all margins for all five rows.
The certificate uses the original 128-point grid, 50 interval digits, and
unchanged strip bounds; the degree increase requires no checker relaxation.

The best witness has `eta = 0.51798422`, `tau = 0.48972396`,
`rho = 0.05776256`, `J = 17`, `d = 0.01029604`,
`delta = 0.0105116329`, `gamma = 0.4894583671`, and `chi = 0.51788422`.
It preserves `epsilon = 0.00000001`, `zeta = 0.00001`, and the variance reserve
`eta-chi-48 epsilon = 0.00009952`. The preliminary gap, representation gap,
and preliminary/checksum assembly gaps are unchanged.

| Component | Certified upper bound or exact exponent |
| --- | ---: |
| Filtered lists | < 0.489707854 |
| Visited pairs | < 0.489713855 |
| List-generation setup | 0.30776256 |
| Preliminary search | 0.48971396 |
| Main work plus `20 epsilon` | 0.48972416 |
| Prime testing | 0.25899211 |

The coverage margin exceeds `0.000026`; filtered-list and visited-pair margins
exceed `0.000016` and `0.000010`, respectively. Their rounded margins are a little
smaller than the predecessor's, while preserving its prescribed minimum reserve
policy. Block-space slack is `0.01803648 > 0`, down from `0.0182024`.

## Differentiating the moments

Use the score definitions in the
[shared derivation](../filter-budgets/README.md#moment-certificate). Write
`C(A,B)(x) = int A(x+w) B(w) dw`, `F=2^r`, `G=2^s`,
`P=C(F,G)`, `Y=C(Fr,G)+C(F,Gs)`, and `z=log2(Z_f Z_g)`.
For an `r` coefficient with basis function `c_j=cos(jx)`,

```text
dF = ln(2) F c_j
d(Fr) = F c_j (1 + ln(2) r)
dP = C(dF,G)
dY = C(d(Fr),G) + C(dF,Gs)
dz = int F c_j / Z_f
dH = kappa int F^kappa c_j / int F^kappa - kappa dz.
```

The `s` derivatives are symmetric. Holding coefficients fixed,

```text
dH/dkappa = int F^kappa r / int F^kappa
         + int G^kappa s / int G^kappa - z.
dL = dP/(ln(2) P) - dz
dB = dY/P - Y dP/P^2 - dz
du = dL + (kappa-1)dB + B d(kappa) - dH
dv = (2-kappa)dB - 2dL - B d(kappa) + dH.
```

For `L_i=log2(I_i)`, differentiating a log integral gives the normalized
integrand-weighted average of the derivative of its base-two exponent.
Thus the two gradients average `-8.5 du` and `8.3(a dv-du)`.
The implementation batches the correlation derivatives using NumPy FFTs.
Central differences of the independent shared scorer check every coefficient
and `kappa` at both 128 and 256 points; a uniform-density control has known
moments and zero first derivatives.

At temperature `T`, the surrogate is
`T log(exp(L_1/T)+exp(L_2/T))`, with natural logarithm here. It lies between
`max(L_1,L_2)` and `max(L_1,L_2)+T log(2)`; its gradient is the softmax-weighted
combination of the two gradients. The temperatures are `0.001`, `0.0001`,
`0.00001`, and `0.000001`. Each stage permits 200 steps and 24 backtracking
trials per step, with an overall cap of 12,000 evaluations per degree.
Coefficients stay within `0.35` of their start, and `kappa` stays in `[0.95,1.12]`.
Only degrees ten and fourteen were searched.

## Retuning and replay

[retune.py](retune.py) uses the predecessor's eliminated-budget formulas, with
fixed deficiency, coverage, representation, list, and pair reserves
`0.00002`, `0.000026`, `0.00003`, `0.000016`, and `0.00001`.
It balances list and pair bounds at 17 blocks, then bisects the local moment
boundary over `eta` in `[0.5179,0.51805]` for 40 steps. Decimal rounding restores
the strict margins before independent certification. This is one bounded
parameter slice, not a global optimum.

From the repository root, with the existing requirements installed:

```sh
python3 -B research/minimax-filters/optimize.py --degree 10 --output .context/minimax-degree10.json
python3 -B research/minimax-filters/optimize.py --degree 14 --output .context/minimax-degree14.json
python3 -B research/minimax-filters/compare.py --output .context/minimax-mutation-control.json
python3 -B research/minimax-filters/retune.py --output .context/minimax-retuning.json --witnesses-output .context/minimax-witnesses.json
python3 -B research/filter-budgets/verify.py --witnesses research/minimax-filters/witnesses.json --output .context/minimax-certificate.json
make check
git diff --check
```

Create `.context/` first if absent. The retuner reads the retained discovery
files; to retune a fresh optimizer result, write that result to the corresponding
`research/minimax-filters/discovery-degree10.json` or `discovery-degree14.json`
path instead. Fresh results can differ in their last digits across environments;
the exact retained witnesses are the certificate inputs.

The two descent campaigns used 803 evaluations, the mutation controls 8,002,
and the retuning replays 86. Discovery receipts record the configured limits,
actual counts, environment, coefficients, and input/code hashes. The search
evidence remains floating point; the continuous certificate is separate.
Package tests cover derivative agreement, invalid inputs, exact evaluation
caps, deterministic bounded descent, retained reserves, and rejected caps.

## Next question

At the retuned degree-fourteen witness, moment 2 has about `0.00000104` slack
below `chi`, while moment 1 has about `0.00000357`. Can the exposed Jacobian
optimize the complete retuned frontier directly, allowing this unequal slack,
rather than minimizing the worst moment at the fixed predecessor parameters?
