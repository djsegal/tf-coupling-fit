# tf-coupling-fit

Per-edge transcription factor coupling strength inference for *S. cerevisiae*,
intended as the transcription-module input for a whole-cell model.

**Project context.** A team of ~20 students is rebuilding a whole-cell model
of budding yeast. The transcription module needs to know, for every gene with
known regulators, how strongly each TF couples to that gene's transcription
rate. This repo provides those coupling values (~3,900 non-zero edges across
~2,400 substrate genes and ~150 TFs), plus the analysis notebook that
generated them.

## How to use this repo

**If you just want the fitted alpha values**, go straight to `output/tf_handoff/`.
The files are committed to the repo so you don't need to run anything. The
README inside that folder explains the conventions and units.

**If you want to inspect or rerun the analysis**, open
`notebook/transcription_fit.ipynb`. End-to-end runtime is about a minute on a
modern laptop.

**If you don't want to run the [notebook](http://github.com/djsegal/tf-coupling-fit/blob/main/notebook/transcription_fit.ipynb) yourself**, the reference figures in
`figures_reference/` show what each plot should look like.

## Quick start (re-running the notebook)

You need Julia 1.8 or later.

```bash
git clone https://github.com/djsegal/tf-coupling-fit.git
cd tf-coupling-fit
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Then open `notebook/transcription_fit.ipynb` in Jupyter (with the Julia kernel)
and Run All. The notebook reads from `data/` and overwrites the contents of
`output/tf_handoff/` plus `output/tf_handoff.xlsx`. The .zip archive of the
handoff is in `.gitignore` since it can be re-created from the folder contents.

## What's in this repo

```
tf-coupling-fit/
|-- README.md                   <- you are here
|-- LICENSE                     <- MIT
|-- Project.toml                <- Julia dependencies (pinned)
|-- .gitignore
|-- notebook/
|   `-- transcription_fit.ipynb    <- the analysis, end to end
|-- data/
|   |-- WT_unstressed_readspermillionreads.csv   <- Teufel 2019 RNA-seq
|   `-- TF_von_Teufel.csv                        <- Teufel 2019 TF network
|-- figures_reference/
|   |-- 01_nan_distribution.png        <- preview of the notebook plots,
|   |-- 02_tf_coverage.png             <-   so you can see what to expect
|   |-- 03_canaries_max_nans_4.png     <-   without running anything
|   |-- 04_cln2_strict_vs_relaxed.png
|   |-- 05_lcurve.png
|   `-- 06_fit_quality.png
|-- output/
|   |-- tf_handoff.xlsx                <- everything below, as one workbook
|   `-- tf_handoff/                    <- the handoff folder
|       |-- README.md                   <-   conventions, units, load snippet
|       |-- metadata.json               <-   fit params, source data, counts
|       |-- alpha_matrix.csv            <-   dense n_substrates x n_tfs matrix
|       |-- alpha_edges.csv             <-   long form, all candidate edges
|       |-- tf_network_candidate.csv    <-   structural network from Teufel
|       `-- tf_network_fitted.csv       <-   only edges with abs(alpha) > 0
`-- transcription_multiplier/       <- turning the couplings into a cell-cycle rate
    |-- README.md                   <-   how to use it + validation summary
    |-- multiplier.jl               <-   drop-in mean-preserving multiplier
    |-- refit.jl                    <-   fit-reproduction + CV harness (JuMP/HiGHS)
    |-- runtests.jl                 <-   tests (reproduces the committed fit exactly)
    `-- data/                       <-   TF means, q_x scores, validation tables
```

## How to use the output in the whole-cell model

The fitted coupling values are in `output/tf_handoff/`. Pick the file that
matches how your downstream code wants the data:

| File                       | What to do with it                                                                                                                                 |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `alpha_matrix.csv`         | Load as a dense `n_substrates x n_tfs` matrix; multiply with the current TF mRNA vector (delayed by tau) to get per-substrate transcription rates. |
| `tf_network_fitted.csv`    | Long-form edge list of fitted couplings, if you prefer sparse representation.                                                                      |
| `tf_network_candidate.csv` | Structural network from Teufel, if you want to refit alpha with a different method.                                                                |
| `alpha_edges.csv`          | Long-form including the L1-zeroed edges. Useful if you need to know which candidate edges the fit rejected.                                        |
| `metadata.json`            | Fit parameters, timestamp, counts. Read it before trusting anything.                                                                               |
| `README.md`                | Conventions, units, what `alpha = 0` means.                                                                                                        |

`output/tf_handoff.xlsx` is the same six files bundled into one Excel
workbook, one sheet per file. Convenient for emailing.

### Driving a cell-cycle dependent rate

If you want to turn these couplings into a time-varying transcription rate for a
dynamic model, see `transcription_multiplier/`. It provides a mean-preserving
multiplier `k_x(t) = k_x * (1 + q_x * sum_i alpha_i (TF_i(t)/TF_i_mean - 1) / sum_i |alpha_i|)`,
the per-TF cell-cycle means it needs (computed on the interpolated grid, so the
multiplier preserves each gene's average rate to machine precision), per-gene
periodicity weights `q_x`, a drop-in Julia implementation, a fit-reproduction
harness with tests (reproduces the committed fit exactly), and a six-source
structural corroboration of the network. The multiplier is a steady-state,
mean-preserving relative-rate descriptor, not a dynamic mRNA predictor. See its
README for the formula, data dictionary, and validation summary.

Example Julia code (from inside `output/tf_handoff/`):

```julia
using DataFrames, CSV, JSON3

meta  = JSON3.read(read("metadata.json", String))
alpha = CSV.read("alpha_matrix.csv", DataFrame)
edges = CSV.read("tf_network_fitted.csv", DataFrame)

# Get CLN2's regulators and coupling strengths:
edges[edges.substrate .== "CLN2", :]
#  substrate    tf      alpha
#  CLN2         STE12   0.532...
#  CLN2         SWI4    1.375...
```

## The model and the math

For each substrate gene `i`:

$$S_i(t) = \sum_{j \in R(i)} \alpha_{ij} \cdot TF_j(t - \tau)$$

We fit alpha by minimizing `L1(residual) + w2 * L1(alpha)` via a linear
program (HiGHS solver, via JuMP). One LP per substrate, all ~2,400
substrates in a few seconds.

The notebook walks through how we picked the regularization weight, why we
default to a NaN tolerance of 4 (with linear interpolation), and why we
intentionally omit the intercept term. See sections 8, 9, and 10 of the
notebook.

The handoff couplings in `output/tf_handoff/` use the time delay `tau = 20` min
and are the canonical, Teufel-only fit. Two variants are also available: a joint
multi-dataset fit (across several cell-cycle RNA-seq time courses) that recovers more
topology and transfers better across conditions, at the cost of noisier signs; and a
leakage-free nested-CV setting `(w2 = 1, tau = 15)` recommended for new fits. The
shipped `tau = 20` couplings remain the reference.

## What this is **not**

- **Not a publication-ready fit.** Cross-validation, bootstrap CIs, and a
  random-network baseline are all on the v2+ list. See section 14 of the
  notebook.
- **Not a complete TF list.** Only TFs in the Teufel 2019 supplementary
  table 4 are considered. Known regulators absent from that table (e.g.
  MBP1 -> CLN2) won't be captured.
- **Not protein-level.** We use TF *mRNA* as a proxy for active TF protein.
  Post-transcriptional regulation, phosphorylation, and localization are
  out of scope here.

## Data attribution

RNA-seq and TF network: Teufel L. *et al.* (2019). "A transcriptome-wide
analysis deciphers distinct roles of G1 cyclins in temporal organization
of the yeast cell cycle." *Scientific Reports* 9: 3343.
[doi:10.1038/s41598-019-39850-7](https://doi.org/10.1038/s41598-019-39850-7).

## Citation

If you use these alpha values in your work, please credit the source paper
(Teufel et al. 2019) and link back to this repo so others can find the
fitting methodology.
