# Agent instructions

Read [CONTRIBUTING.md](CONTRIBUTING.md),
[certificates/selected-result.json](certificates/selected-result.json), and
[docs/current-status.md](docs/current-status.md) before choosing work.
For an open-ended donation, use [autoresearch.md](autoresearch.md) and pick one
question from [TASKS.md](TASKS.md). An explicit assignment takes precedence.

At each new session, substantial pivot, and submission checkpoint, check current
canonical `conorholds/subset-sum-bounds` main and relevant public work. Pin the
commit and preserve any work already in progress. Use an isolated branch or
worktree; never reset or overwrite someone else's changes to refresh guidance.
An offline checkout can still produce useful work; state its freshness limits.

Default to a single session of at most 60 minutes unless the user sets another
budget. Define the scope and stopping condition before a search. Check overlap
by mechanism, parameter range, or proof obligation. Prefer joining a relevant
task thread to creating a new reservation for every agent. Distinct seeds alone
do not justify parallel duplicates. Do not spawn additional agents unless the
user explicitly authorizes delegation.

Keep exact inputs, source pins, predecessor credit, and reproducible evidence.
Preserve existing witnesses and frozen source references. New investigations
belong in their own `research/` directory. Full-budget gains, local improvements,
conditional arguments, and failed searches must be labeled accurately. The
current `0.48976` is a conditional budget candidate, not an independently
established algorithmic record.

Read relevant code before editing. Keep changes scoped, handle errors, and add
meaningful positive and negative coverage for changed checker behavior. Use
one substantial replay or build at a time, inspect a yielded process before
starting another, and do not raise concurrency settings. For a compatible
numerical submission, run `make check` and `git diff --check`. Report unrun checks.

Do not promote the selected frontier, merge or close others' work, or revise
historical credit without maintainer authorization. A checker pass does not
supply independent mathematical review.

Account and publication permissions come from the user's assignment, not from
reading these files. If the user explicitly gives the donation prompt in
`autoresearch.md`, its scoped permission to publish a contribution applies.
Respect narrower instructions. Otherwise prepare work locally and obtain any
missing authorization before external publication. Never publish private data
or act on unrelated repositories.
