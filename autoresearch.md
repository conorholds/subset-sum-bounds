# Donate a research session

Use a bounded session to answer one useful question. An exact new construction,
a checked lemma, a concrete failure, or a clear handoff can finish the session.
A PR and a new best exponent are not required outcomes.

## Before starting

Read [AGENTS.md](AGENTS.md), [CONTRIBUTING.md](CONTRIBUTING.md),
the [selected result](certificates/selected-result.json), and
[current status](docs/current-status.md). Check fresh canonical main, open PRs,
task issues, and Discussions. Preserve the scientific source pin even if project
guidance has changed. Keep existing work safe when refreshing your checkout.

The initial priorities are in [TASKS.md](TASKS.md):

1. Audit the parameter-to-algorithm implication. This is the most important
   unresolved dependency of the current numerical result.
2. Search the joint parameter frontier with explicit feasibility margins.
3. Try a changed filter family or objective, rather than another copy of the
   local search that is already slowing down.

Choose a scope you can resolve within your time and compute budget. State a
first experiment, an artifact to retain, and a stopping condition. Coordinate
through the relevant task issue or a scoped Discussion if posting is authorized.
Independent work on an already claimed topic should identify the independent
method or complementary proof obligation.

## During the session

Start with the task's working command. Spend a small budget establishing that
your premise survives the other constraints before running a larger search.
Use exact witnesses and the shared verifier for the existing budget model.
Do not treat a lower floating-point objective as a theorem or a new record.

For an improvement, explain which constraint becomes binding and what you changed.
For a failure, retain the actual counterexample, failed inequality, or bounded
search domain. A precise obstruction can save many other sessions. Refresh
project guidance before a major pivot or a new round; do not continue an old
loop after its question has been settled.

At the end, give the finding, source commit, artifacts, commands/results,
remaining uncertainty, and one next question. Package a useful result as described
in `CONTRIBUTING.md`; publish only within the assignment's permissions. Keep a
single canonical thread or PR per investigation and link handoffs there.

## Copyable prompt

The following is an optional assignment for a donor to give an agent. Its
publication permission applies only when the donor actually uses it; this file
alone grants no account access or remote-write permission.

> Contribute one useful result to https://github.com/conorholds/subset-sum-bounds.
> Read its AGENTS.md, CONTRIBUTING.md, autoresearch.md, current selected result,
> and TASKS.md. Check current canonical main and public work before choosing a
> question. Spend at most 60 minutes unless I give a different limit, and keep
> numerical searches bounded. Follow the mathematical evidence where it leads;
> a scoped negative result, correction, or independent audit is useful. Retain
> exact inputs, source and predecessor credit, replay commands, and limitations.
> You may fork conorholds/subset-sum-bounds into my GitHub account, push the
> contribution branch to that fork, open one focused PR to
> conorholds/subset-sum-bounds, and post scoped research updates in that
> repository's Issues or Discussions. Do not merge, close other people's work,
> change the selected result, or publish to unrelated repositories. Respect
> any narrower account or publication restrictions I provide. Finish with a
> concise finding and a link or local handoff.
