# Subset Sum Bounds

[![Verify](https://github.com/conorholds/subset-sum-bounds/actions/workflows/verify.yml/badge.svg)](https://github.com/conorholds/subset-sum-bounds/actions/workflows/verify.yml)

An open research project to improve the exponent in worst-case Subset Sum,
one reproducible contribution at a time. Humans and coding agents are welcome.

Given a list of integers and a target, Subset Sum asks whether some subset adds
up to that target. We are exploring the randomized construction in
[OpenAI's family 138 manuscript](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026).
The shared numerical objective is a smaller exponent `alpha` in `2^(alpha n)`,
where `n` counts the integers, while paying every cost in the construction.

**Current scope: interval-checked numerical witnesses in a conditional budget
model.** The full algorithmic argument and its parameterization still need
independent review. This is an initial research candidate, not an established
new Subset Sum running-time record or a practical solver benchmark.

## Start contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md), then pick one bounded question from
[TASKS.md](TASKS.md). If you want to donate an agent session, use
[autoresearch.md](autoresearch.md). Coordinate in
[Discussions](https://github.com/conorholds/subset-sum-bounds/discussions) and
check [open work](https://github.com/conorholds/subset-sum-bounds/issues) and
[pull requests](https://github.com/conorholds/subset-sum-bounds/pulls) first.

| You could contribute | What would help |
| --- | --- |
| A parameter improvement | An exact witness with a lower complete budget, checked against every constraint. |
| A different filter | Coefficients, a reproducible search, and interval bounds; local moment gains can be combined with later retuning. |
| A mathematical audit | A checked implication, a missing hypothesis, or a concrete counterexample in the budget derivation. |
| An independent verification | A separate numerical implementation or proof of a specific certificate step. |
| A useful negative result | A precise family, search budget, failure reason, and a next question. |

You do not need a new best exponent to make a useful contribution. Preserve
predecessor credit and describe exactly what your result establishes.

## Initial checkpoint

| Witness | Filter degree | Outer exponent budget `alpha` |
| --- | ---: | ---: |
| Source construction | 6 | 0.49000 |
| Retuned parameters with source filters | 6 | 0.48981 |
| Retuned parameters with new filters | 10 | **0.48976** |

These budgets are conditional on the manuscript's qualitative algorithmic
arguments and our generalization of its fixed constants. The source paper also
has an internal guarded bound of `2^(0.489995n) poly(n+b)`; the table's source
row uses its headline exponent. Smaller budgets describe asymptotic exponents,
not measured speedups on small instances.

[Selected result](certificates/selected-result.json) ·
[Current status and open obligations](docs/current-status.md) ·
[Derivation and experiments](research/filter-budgets/README.md) ·
[Exact witnesses](research/filter-budgets/witnesses.json) ·
[Certificate](research/filter-budgets/CERTIFICATE.json)

The selected result is the single starting point for comparisons. An unchecked
PR or an old README does not replace it. Our initial searches found two gains
that compose, but simple follow-up searches slowed down; sustained headroom is
an open research question.

## Reproduce locally

Use Python 3.12 or later; the initial local replay used Python 3.14.

```sh
git clone https://github.com/conorholds/subset-sum-bounds.git
cd subset-sum-bounds
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make check
```

`make test` runs source controls and rejection tests. `make verify` replays all
`research/*/witnesses.json` packages sequentially, checks their committed
receipts, and validates the selected record. Fresh receipts go in `.context/`.
The initial four-witness interval replay takes roughly 13 seconds on the
maintainer's machine; timing varies by environment.

For a bounded new experiment:

```sh
make search
# Or choose an explicit search budget and write an exploratory candidate:
python research/filter-budgets/search.py --start combined --parameters combined --degree 14 --steps 4000 --output .context/candidate.json
```

Search output is exploratory. A contribution needs a separate certificate and
an explanation of how its improvement affects the complete budget.

## License

[Apache 2.0](LICENSE). See [NOTICE](NOTICE) for source attribution.
