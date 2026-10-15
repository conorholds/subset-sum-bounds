# Current status

The initial selected candidate is `combined`, with outer budget
`alpha = 0.48976`. Its status is **interval-checked and conditional**;
independent mathematical review is pending. Resolve its files and command
through [selected-result.json](../certificates/selected-result.json).

The seed work reproduced the source's printed scalar quadrature controls,
retuned its original parameters, and found a degree-ten filter that permits a
further small retuning gain. All four retained ablations pass interval moment
bounds and the modeled budget constraints. The evidence is in
[the research package](../research/filter-budgets/README.md).

What the current checks establish:

- Continuous filtering moment enclosures, including inner and outer quadrature
  errors and the complex zero-free strip used for logarithms.
- Strict margins for the rederived sufficient-budget inequalities.
- Exact binding of witnesses, checker files, dependency versions, and the pinned
  mathematical source manifest.

What remains open:

- An independent audit that the parameterized inequalities suffice for every
  part of the source algorithm, including probability laws, all-size cutoffs,
  precision, and the complete operation count.
- Independent validation beyond the current `mpmath.iv` arithmetic and analytic
  certificate implementation; a formal proof is not supplied.
- Evidence of sustained structural headroom beyond the initial tuning gains.
  Practical crossover sizes and a full executable Subset Sum solver are not known.

The source manuscript is an input to this project, not a theorem independently
established by our certificate replay. Scientific priority and worldwide
optimality have not been established. A newly discovered gap should be recorded
even if it invalidates the initial candidate; fixing the shared understanding
is a useful contribution.

New numerical packages are proposals until reviewed. A maintainer can update
the selected record after replay and independent review of the changed argument,
while preserving the prior evidence and contributors' credit. CI success alone
does not promote a candidate. See [CONTRIBUTING.md](../CONTRIBUTING.md).
