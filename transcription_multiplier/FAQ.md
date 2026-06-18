# FAQ

## Which file do I load?

For the multiplier you need two things: the fitted edges and the TF means.
`load_handoff()` returns both from `data/tf_network_fitted.csv` and
`data/tf_means.csv`. If you prefer a dense matrix, the handoff folder also ships
`alpha_matrix.csv` (n_substrates x n_tfs). `alpha_edges.csv` is the long form
including the L1-zeroed candidate edges, useful if you want to know which edges the
fit rejected.

## What does `alpha = 0` mean?

The L1 penalty drove that candidate edge to zero — the regulator carries no
detectable signal for that substrate under the fit. Zero edges are omitted from
`tf_network_fitted.csv` and present in `alpha_edges.csv`.

## What happens for a gene with no regulators?

`multiplier` returns `1.0` (the baseline rate is unmodulated). Same for a gene
whose couplings are all zero. So you can call it for every gene without special-casing.

## Teufel-only vs joint couplings — which should I use?

The shipped `tf_network_fitted.csv` (Teufel-only, tau = 20) is the default
within-condition multiplier; its edge signs are cleaner, which is what the
multiplier acts on. The joint multi-dataset couplings
(`joint_multidataset_alpha.csv`) recover more topology and transfer better across
conditions, at the cost of noisier signs — an opt-in cross-condition tool.

## How do I save couplings back out programmatically?

`save_handoff(dir, edges, means)` writes the dicts back to the same CSV schema
(round-trips with `load_handoff`). `export_handoff_json(path, edges, means)` emits
one portable JSON document for non-Julia consumers, with no extra dependency.

## Is `multiplier` safe to call from multiple threads?

Yes. It reads only its arguments and holds no shared mutable state, so concurrent
calls on a shared `(edges, means)` are race-free. `runtests.jl` asserts this with a
`Threads.@threads` sweep over all genes.

## How do I reproduce the headline numbers?

`julia transcription_multiplier/fetch_datasets.jl --gold-only` fetches two small
public files, then `runtests.jl` runs the joint AUROC + paired-DeLong test
(0.804 -> 0.899, +0.071 [0.052, 0.091], p = 1.4e-12). Without the fetch the test
skips and the rest of the suite still runs.
