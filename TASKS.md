# Initial research tasks

These are starting questions, not promises of improvement. Read the
[selected result](certificates/selected-result.json) and check
[Issues](https://github.com/conorholds/subset-sum-bounds/issues),
[Discussions](https://github.com/conorholds/subset-sum-bounds/discussions), and
[PRs](https://github.com/conorholds/subset-sum-bounds/pulls) before claiming work.
Claims coordinate effort; they do not exclude independent approaches.

## 1. Audit the generalized budget argument

**Priority:** highest. **Good fit:** mathematical reading and adversarial review.

The current verifier checks moments and 14 budget margins. Its connection to
an actual algorithm uses a rederivation of the manuscript's fixed constants.
Independently check one complete interface: preliminary deficiency versus
coverage, alias-list and pair costs, or checksum variance versus the moments.

Start with the conditional model in
[the derivation](research/filter-budgets/README.md) and the exact sections in
[SOURCE.json](research/filter-budgets/SOURCE.json).

```sh
make verify
```

Choose one implication, map every parameter to its source use, and try to find a
missing constraint or conditioning step. An acceptable result is a written
proof of that implication with its remaining dependencies, or a precise gap and
counterexample. Replaying the same arithmetic alone does not resolve this task.
Stop after resolving that interface or isolating a specific blocker.

## 2. Map the joint parameter frontier

**Good fit:** constrained optimization and asymptotic cost analysis.

The initial parameter tuning used a coarse screen and manual conservative
rounding. It is not an optimality result. Hold the source or combined filter
fixed, and search `eta`, `tau`, block density/count, and the deficiency/auxiliary
modulus budgets together. Include the 16-versus-17-block boundary and check
which constraint actually binds. Keep fixed positive rounding and variance
reserves; do not earn an apparent gain by silently dropping them.

```sh
python research/filter-budgets/verify.py
```

Use `budget_margins` in
[certify.py](research/filter-budgets/certify.py) as the acceptance check, not
as a differentiable search objective. Provide a bounded search with exact
retained witnesses and a complete certificate. A rigorous ceiling for a named
parameter slice or a sensitivity map is also useful. Stop if the claimed gain
disappears when every cost and strict margin is restored.

## 3. Change the filter search

**Good fit:** numerical methods, harmonic analysis, or new constructions.

The initial degree-ten filter improves both original moments. Another 4,001
mutation evaluations at degrees 10 and 14 made only tiny gains. Try a genuinely
different family, joint objective, initialization, or structural argument.
Explain why the chosen change could improve the complete frontier.

```sh
make search
```

Start from the combined design in
[witnesses.json](research/filter-budgets/witnesses.json). For the existing cosine
schema, keep coefficients exact and obtain continuous interval bounds with the
shared verifier. Other positive densities need a derivation and their own
certificate path; sampled grid agreement alone is insufficient. Retain a local
moment improvement even if parameter retuning is a separate handoff. Record the
search budget and stop if the proposed family cannot be certified or does not
survive the other costs.

## Smaller ways to help

An independent derivation of the zero-free strip bound, a checker failure case,
a clearer source-to-constraint explanation, and an attribution correction are
all welcome. Name the exact question and expected evidence. Avoid creating
routine metadata or documentation changes merely to produce a PR.
