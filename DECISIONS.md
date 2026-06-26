# Design decisions

Short rationales for the choices that shape the deliverable. Each is the kind of
question a downstream integrator is likely to ask.

## 1. Deviation form: normalize by `sum |alpha|`, not `sum alpha`

The mean-preserving multiplier is
`M_x(t) = 1 + q_x * sum_i alpha_i (TF_i(t)/TF_i_mean - 1) / sum_i |alpha_i|`.
Dividing by the signed `sum_i alpha_i` (the obvious construction) is unsafe:
activator and repressor contributions cancel and drive the denominator toward
zero, and net-repressor genes make it negative — 17 of 2,406 genes go negative
over the cell cycle. The absolute-value denominator removes every blow-up and
negative-rate event, preserves the time average exactly for any mix of signs, and
reduces to the obvious form when all couplings are positive.

## 2. `M_x = 1` for genes with no usable regulators

A substrate with no fitted regulators, or all-zero `alpha`, returns `1.0`
(`sum |alpha_i| = 0`), leaving its baseline rate unmodulated. Enforced directly in
`multiplier` so partial input never errors.

## 3. Dependency-light, serialization optional

The multiplier needs only `CSV` + `DataFrames`. `save_handoff` round-trips through
the same CSVs; `export_handoff_json` is a hand-written, write-only JSON emitter, so
the package gains no Arrow/JSON dependency. Anything heavier is the consumer's call.

## 4. TF mRNA as a proxy for TF activity

The couplings regress substrate mRNA on regulator *mRNA*. This is a deliberate,
stated approximation: many cell-cycle regulators are timed post-translationally,
which the proxy cannot resolve (see `known_limitations.md`). The multiplier
exploits the periodic mRNA program, not faithful TF-protein dynamics.

## 5. L1 regression on a curated candidate set

One L1 linear program per substrate (sparse, signed couplings) over the regulators
in the Teufel candidate network. The value over naive co-expression is signed,
sparse, transferable couplings — not a better gold-standard ranking, which on the
curation-aligned rectangle is dominated by candidate-set membership.

## 6. `tau = 20 min` Teufel-only fit is the shipped default; joint is opt-in

The shipped couplings use a 20-minute TF→target delay on the Teufel course. The
joint multi-dataset refit recovers more topology and transfers better across
conditions but at the cost of noisier signs, so it is offered as an opt-in tool
rather than the default. A leakage-free nested-CV setting (`w2 = 1, tau = 15`) is
recommended for *new* fits.
