# Transcription factor coupling handoff

Fitted coefficients for `S_i(t) = sum_j alpha_ij * TF_j(t - tau)` on the
Teufel 2019 *S. cerevisiae* WT unstressed time course. `tau = 20` min.
Both `S` and `TF` are in molecules per cell (RPM rescaled by `60000 / 1e6`).
The fit is doc-literal: no intercept, no degradation term. See section 10 of the
generating notebook for the intercept trade-off.

`tau = 20` min is the canonical handoff used here. A leakage-free nested
cross-validation prefers `(w2 = 1, tau = 15)`; that setting is recommended for any
new fit, while these `tau = 20` couplings remain the reference. These are the
Teufel-only couplings; a joint multi-dataset variant (fit across several cell-cycle
RNA-seq time courses) recovers more topology and transfers better across conditions,
at the cost of noisier activator/repressor signs. To turn these couplings into a
mean-preserving cell-cycle rate multiplier, see the `TranscriptionMultiplier` package in this repository; its
per-TF means are computed on the interpolated grid so the multiplier preserves each
gene's average rate to machine precision.

## Files

- `metadata.json`: fit parameters, source data, counts, timestamp.
- `alpha_matrix.csv`: dense `n_substrates x n_tfs` matrix, header row gives TF order.
- `alpha_edges.csv`: long form `substrate, tf, alpha`. All candidate edges, including `alpha = 0`.
- `tf_network_candidate.csv`: structural connectivity from Teufel, filtered to measured TFs and substrates. Use to refit alpha with a different method.
- `tf_network_fitted.csv`: only edges with `abs(alpha) > 0` after L1 sparsification. Use if you trust the fit.

## Conventions

- `alpha = 0` means "fit drove it to zero", not "no edge in the network". Absent edges do not appear in any file.
- Sign of alpha is free: positive activates, negative represses.
- Linear (not log) space throughout.

## Load in Julia

```julia
using DataFrames, CSV, JSON3

meta  = JSON3.read(read("metadata.json", String))
alpha = CSV.read("alpha_matrix.csv", DataFrame)
edges = CSV.read("tf_network_fitted.csv", DataFrame)

# alpha[alpha.substrate .== "CLN2", :SWI4]   # about 1.38
```

## Provenance

Generated from `transcription_fit.ipynb` at the timestamp in `metadata.json`.
Source: Teufel et al. 2019, *Scientific Reports* 9:3343,
[doi:10.1038/s41598-019-39850-7](https://doi.org/10.1038/s41598-019-39850-7).
