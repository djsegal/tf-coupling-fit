# transcription_multiplier

Turn the fitted TF couplings in `../output/tf_handoff/` into a cell-cycle dependent
transcription rate for a dynamic (e.g. whole-cell) model, without changing each
gene's average rate. Julia-native, matching the stack of the fit itself.

For a gene `x` with literature rate constant `k_x` and fitted couplings `alpha_xi`
to its regulators `i`:

```
k_x(t) = k_x * M_x(t)
M_x(t) = 1 + q_x * ( sum_i alpha_i * (TF_i(t)/TF_i_mean - 1) ) / sum_i |alpha_i|
```

`TF_i(t)` is the current level of regulator `i`, `TF_i_mean` its cell-cycle average,
and `q_x in [0,1]` an optional per-gene cell-cycle weight. This "deviation form":

- **preserves the mean** (`<M_x>_t = 1`) for any mix of activators and repressors
  and any `q_x`, so the long-run rate is unchanged;
- **never blows up** (`sum |alpha_i|` cannot vanish); the older signed denominator
  drove 17 of 2,406 genes negative, to as low as -8.5;
- is **identical** to dividing by `sum alpha_i` for the 2,322 genes with no
  repressor; only the 84 repressor-carrying genes differ.

## Use it (Julia)

```julia
include("multiplier.jl")
edges, means = load_handoff()                 # reads data/
levels = Dict("SWI4" => 220.0, "STE12" => 80.0)  # current TF amounts (any consistent unit)
M = multiplier("CLN2", levels, edges, means; q = 1.0)
k_eff = k_base * M
```

Units cancel in `TF_i(t)/TF_i_mean`, so feed TF levels in whatever unit the model
stores. **Use `q_x` rather than `q = 1`** (see validation below).

## Reproduce the fit and run the tests (Julia)

```bash
julia --project=.. transcription_multiplier/runtests.jl
```

`refit.jl` re-implements the fitting linear program (JuMP/HiGHS) and reproduces the
committed `alpha` exactly (max |delta alpha| = 0). `runtests.jl` checks that, plus
the multiplier properties (unit mean, non-negativity, all-activator identity,
`q_x` amplitude scaling, and machine-precision mean preservation across all genes).
All 28 checks passing.

## Data dictionary (`data/`)

| File | Contents |
| --- | --- |
| `tf_means.csv` | per-TF cell-cycle mean (the required multiplier input), computed on the interpolated 22-point grid so the multiplier is mean-preserving to machine precision |
| `tf_network_fitted.csv` | fitted non-zero couplings `(substrate, tf, alpha)` |
| `qx_scores.csv` | per-gene `q_x` from Cyclebase 3.0 (rank, p-value, peaktime) + three [0,1] transforms; the rank-linear transform (`q_rank`) is the default |
| `tf_normalization_table.csv` | per-gene `sum alpha`, `sum |alpha|`, flags, min multiplier |
| `multiplier_examples.csv` | signed vs deviation multiplier over the cell cycle |
| `multisource_support.csv` | each fitted edge vs six independent networks + sign agreement |
| `augmented_network_candidate.csv` | optional v2 network (Teufel + YEASTRACT-strict + literature) |
| `augmented_alpha_fitted.csv` | `alpha` refit on the augmented network |
| `union_network_alpha_fitted.csv` | `alpha` refit on the union of all sources (an alternative coupling set) |
| `ground_truth_cellcycle_edges.csv` | 62 curated, cited textbook cell-cycle edges |

## What the multiplier is (and is not)

The multiplier is a **steady-state, mean-preserving descriptor** of how a gene's
transcription rate leans across the cell cycle. It is a relative-rate prior, not a
dynamic predictor: in a first-order ODE driven by `k_x(t)`, it does **not** beat a
flat rate at predicting held-out mRNA timepoints (the amplitude cap is structural,
and per-gene phase prediction is at chance genome-wide). Use it to bias rates in the
right direction with the average preserved, not to forecast mRNA dynamics.

## Validation summary

- **Structure:** 86% of fitted edges are corroborated by at least one of six
  independent sources (ChIP binding, two YEASTRACT vintages, a curated benchmark,
  SNAP, textbook edges); 49% by two or more. Scored as a ranked classifier against
  the Inferelator gold standard, the network reaches AUROC 0.815 (precision 0.48, a
  29x lift over the edge prior; recall 0.68). This is structural corroboration
  within a curation-aligned candidate set, not independent validation.
- **Signs:** agree with database labels 71% of the time overall, 94% against
  curated textbook edges. Coupling *magnitude* is fit-derived and does not track
  database presence (corr ~ 0.01). Activator signs are well supported; repressor
  edges are individually lower-confidence, so treat the repressor contribution as a
  sign (a direction), not a magnitude estimate.
- **Mean preservation:** the per-TF means in `tf_means.csv` are computed on the same
  NaN-interpolated 22-point grid the multiplier is evaluated on, so `<M_x>_t = 1`
  holds to machine precision for every gene (earlier means were averaged over each
  TF's non-NaN timepoints, which left up to a ~13% offset on NaN-affected genes).
- **Cell-cycle weighting:** at `q_x = 1` the modulation is modest and skill is
  confined to genuinely periodic genes; weighting by `q_x` improves 93.8% of genes
  and recovers flat-baseline performance, so always apply it. The default `q_x` is
  the rank-linear transform (`q_rank`); `q_sigmoid` is a documented alternative.

## Recommendations

- **Default coupling set: Teufel-only** (`tf_network_fitted.csv`). It carries the
  cleanest signs (94% on curated edges), which is what the sign-aware multiplier
  depends on, so it is the default for within-condition use.
- **Joint multi-dataset couplings** (fit across Teufel plus other cell-cycle RNA-seq
  time courses) recover more topology and transfer better across conditions in a
  leave-one-dataset-out test; the independent signal here is that the joint fit beats
  the single-dataset fit, not the absolute skill. The cost is noisier signs, so
  prefer the joint set for cross-condition transfer and topology, and keep Teufel-only
  for the sign-aware multiplier. The within-source union variant
  (`union_network_alpha_fitted.csv`) is also provided.
- **Time delay `tau`:** the released couplings use `tau = 20` min, which is the
  canonical handoff. A leakage-free nested cross-validation prefers `(w2 = 1, tau = 15)`;
  use that setting for any **new** fit, but the shipped `tau = 20` couplings remain
  the reference.

## License and citation

MIT (see `../LICENSE`). See `CITATION.cff` for how to cite this work. Please also
cite the source data (Teufel et al. 2019, Sci Rep 9:3343) and link back to
`https://github.com/djsegal/tf-coupling-fit`. `CHECKSUMS.sha256` lists the SHA-256
of every released `data/` CSV.
