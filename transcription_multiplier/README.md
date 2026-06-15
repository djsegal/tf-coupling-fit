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
| `qx_scores.csv` | per-gene `q_x` from Cyclebase 3.0 (rank, p-value, peaktime) + three [0,1] transforms |
| `tf_normalization_table.csv` | per-gene `sum alpha`, `sum |alpha|`, flags, min multiplier |
| `multiplier_examples.csv` | signed vs deviation multiplier over the cell cycle |
| `multisource_support.csv` | each fitted edge vs six independent networks + sign agreement |
| `augmented_network_candidate.csv` | optional v2 network (Teufel + YEASTRACT-strict + literature) |
| `augmented_alpha_fitted.csv` | `alpha` refit on the augmented network |
| `union_network_alpha_fitted.csv` | `alpha` refit on the union of all sources (marginally better out-of-sample) |
| `ground_truth_cellcycle_edges.csv` | 62 curated, cited textbook cell-cycle edges |

## Validation summary

- **Structure:** 86% of fitted edges are corroborated by at least one of six
  independent sources (ChIP binding, two YEASTRACT vintages, a curated benchmark,
  SNAP, textbook edges); 49% by two or more. Scored as a ranked classifier against
  the Inferelator gold standard, the network reaches AUROC 0.815 (precision 0.48, a
  29x lift over the edge prior; recall 0.68).
- **Signs:** agree with database labels 71% of the time overall, 94% against
  curated textbook edges. Coupling *magnitude* is fit-derived and does not track
  database presence (corr ~ 0.01).
- **Out of sample:** at `q_x = 1` the multiplier does not beat a flat rate on
  held-out timepoints; skill is confined to genuinely periodic genes. Weighting by
  `q_x` improves 93.8% of genes and recovers flat-baseline performance, so always
  apply it. Treat the cell-cycle modulation as modest and structurally grounded.
- **Network choice:** the union of all sources is marginally better out-of-sample
  than Teufel alone (`union_network_alpha_fitted.csv`); consensus-only and
  ChIP-only are worse. Teufel-based remains the validated default.

## License and citation

MIT (see `../LICENSE`). Please cite the source data — Teufel et al. 2019, Sci Rep
9:3343 — and link back to `https://github.com/djsegal/tf-coupling-fit`.
