# TranscriptionMultiplier.jl

Turn the fitted TF couplings in `output/tf_handoff/` into a cell-cycle dependent
transcription rate for a dynamic (e.g. whole-cell) model, without changing each
gene's average rate. Julia-native, matching the stack of the fit itself.

This repository **is** the `TranscriptionMultiplier` Julia package (`src/`,
`test/`, `docs/`), and it also ships the fitted transcription-factor coupling
handoff for direct use: the ready-to-consume α values are in
`output/tf_handoff/`, and the end-to-end fit that produced them is in
`notebook/transcription_fit.ipynb` (≈1 min, reads `data/`). If you only need the
couplings, take `output/tf_handoff/`; for the cell-cycle rate multiplier, use the
package as below.

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
julia --project=. runtests.jl
```

`refit.jl` re-implements the fitting linear program (JuMP/HiGHS) and reproduces the
committed `alpha` exactly (max |delta alpha| = 0). `runtests.jl` checks that, plus
the multiplier properties (unit mean, non-negativity, all-activator identity,
`q_x` amplitude scaling, and machine-precision mean preservation across all genes).
With no third-party data present, 28 checks pass and the joint-fit check is skipped;
with the gold standard fetched (below), all 35 pass.

## Reproduce the headline cross-dataset result (Julia)

The single-dataset coupling set (`tf_network_fitted.csv`) is the within-condition
default. Fitting jointly across five cell-cycle expression courses
(Teufel + Spellman + Pramila + Orlando + Kelliher) recovers more of the gold-standard
topology. Scored against the signed Inferelator gold standard on the **common
rectangle** (the 82 TF x 606 target block all three coupling sets and the gold
cover, 864 gold-positive edges), area under the ROC rises from

> **Teufel-only AUROC 0.804 -> joint AUROC 0.899**, paired DeLong test
> **+0.071, 95% CI [0.052, 0.091], p = 1.4e-12**.

Regenerate it end to end:

```bash
# 1. fetch the gold standard + ORF name map (two small public files)
julia fetch_datasets.jl --gold-only

# 2. score the SHIPPED coupling sets and run the paired DeLong test
julia --project=. joint_score.jl
#    -> prints 0.804, 0.899, +0.071 [0.052, 0.091], p = 1.4e-12
```

This *scored* path needs no expression data: the joint couplings
(`joint_multidataset_alpha.csv`) and the leave-Kelliher-out couplings
(`loo_kelliher_alpha.csv`) are shipped in `data/`, so an independent reader
reproduces the headline from released files plus the two-file fetch.

To regenerate the joint couplings themselves from raw expression (the deeper
"from scratch" path), fetch all five datasets and re-solve the shared-alpha L1 LP:

```bash
julia fetch_datasets.jl          # all datasets + gold
julia --project=. joint_score.jl --refit
#    re-solves the joint LP (JuMP/HiGHS, ~6 s) and scores those couplings;
#    same 0.804 -> 0.899, +0.071.
```

`joint_fit.jl` builds the joint multi-dataset design and shared-alpha L1 LP (the
same residual-L1 + `w2`*L1(alpha) program as `refit.jl`, with every dataset's rows
stacked). `joint_score.jl` loads the signed gold, scores trapezoid AUROC on the
common rectangle, and runs the fast Sun & Xu DeLong test. The from-raw re-solve
reproduces the shipped joint couplings to < 5e-7 per edge.

The **single-dataset AUROC 0.815** (Teufel-only on its own gold rectangle, the
structural-corroboration number) is likewise re-scorable from released files plus
the fetched gold standard.

### What reproduces what

| Number | Reproduce from | Needs fetch? |
| --- | --- | --- |
| `alpha` (single-dataset couplings), max delta = 0 | `refit.jl` + `data/` (shipped) | no |
| single-dataset AUROC 0.815 | shipped `tf_network_fitted.csv` + gold | gold + SGD only |
| joint AUROC 0.899 / DeLong +0.071 (scored) | shipped joint couplings + gold | gold + SGD only |
| joint couplings themselves, then 0.899 / +0.071 | `joint_fit.jl` re-solve from raw | all five datasets + gold |

### Shipped vs fetched

**Shipped** (in `data/`, checksummed in `CHECKSUMS.sha256`): the single-dataset
couplings, the joint couplings (`joint_multidataset_alpha.csv`), the
leave-Kelliher-out couplings (`loo_kelliher_alpha.csv`), and all derived tables.

**Fetched, not redistributed** (downloaded by `fetch_datasets.jl` into
`data/external/`, each from its own public source and verified by SHA-256):

| File | Source | Citation |
| --- | --- | --- |
| `inferelator_gold_standard.tsv` | Flatiron Institute `inferelator` repo, `data/yeast/gold_standard.tsv` | Tchourine, Vogel & Bonneau 2018, Cell Rep 23:376 |
| `SGD_features.tab` | SGD archive (`sgd-archive.yeastgenome.org`) | Saccharomyces Genome Database |
| `spellman_combined.txt` | Spellman cell-cycle combined table (archived) | Spellman et al. 1998, MBC 9:3273 (PMID 9843569) |
| `pramila.pcl` | SGD archive `GSE4987_setA_family.pcl` | Pramila et al. 2006, Genes Dev 20:2266 (GEO GSE4987) |
| `orlando_setA.pcl` | SGD archive `GSE8799_setA_family.pcl` | Orlando et al. 2008, Nature 453:944 (GEO GSE8799) |
| `GSE80474_Scerevisiae_normalized.txt` | GEO supplementary (`.gz`, decompressed) | Kelliher et al. 2016, PLoS Genet 12:e1006453 (GEO GSE80474) |

If the archived Spellman table is unreachable, `fetch_datasets.jl` prints a manual
fallback (the Bioconductor `yeastCC` package, or the SGD-archive alpha-factor PCL);
the *scored* headline does not depend on it, since the joint couplings are shipped.

## Data dictionary (`data/`)

| File | Contents |
| --- | --- |
| `tf_means.csv` | per-TF cell-cycle mean (the required multiplier input), computed on the interpolated 22-point grid so the multiplier is mean-preserving to machine precision |
| `tf_network_fitted.csv` | fitted non-zero couplings `(substrate, tf, alpha)`; the single-dataset (Teufel-only) default |
| `joint_multidataset_alpha.csv` | joint couplings fit across all five cell-cycle courses (the AUROC 0.899 set); `(substrate, tf, alpha)` |
| `loo_kelliher_alpha.csv` | leave-Kelliher-out joint couplings (leakage-free transfer reference); `(substrate, tf, alpha)` |
| `qx_scores.csv` | per-gene `q_x` from Cyclebase 3.0 (rank, p-value, peaktime) + three [0,1] transforms; the rank-linear transform (`q_rank`) is the default |
| `tf_normalization_table.csv` | per-gene `sum alpha`, `sum |alpha|`, flags, min multiplier |
| `multiplier_examples.csv` | signed vs deviation multiplier over the cell cycle |
| `multisource_support.csv` | each fitted edge vs six independent networks + sign agreement |
| `augmented_network_candidate.csv` | optional v2 network (Teufel + YEASTRACT-strict + literature) |
| `augmented_alpha_fitted.csv` | `alpha` refit on the augmented network |
| `union_network_alpha_fitted.csv` | `alpha` refit on the union of all sources (an alternative coupling set) |
| `ground_truth_cellcycle_edges.csv` | 62 curated, cited textbook cell-cycle edges |
| `multiplier_confidence.csv` | per-gene 0-1 trust score (blend of periodicity, bootstrap robustness, source support, cross-dataset transfer, repressor fragility) + USE/CAUTION/DO-NOT-USE tier and reason; a heuristic flag for downstream use, not a calibrated probability |

### Per-column reference

Every column of every shipped `data/` CSV is described below. The `(substrate, tf,
alpha)` triple has a fixed meaning throughout: `substrate` is the regulated
(target) gene, `tf` is the regulating transcription factor, and `alpha` is the
fitted signed coupling (positive = activator, negative = repressor; zero where an
edge is a candidate but sparsified out). Counts below were verified against the
shipped files on the date of this release.

**Naming carryovers (documented, not yet renamed).** Three column names are kept
for backward compatibility and will be standardized in a future MAJOR release (see
`VERSIONING.md`). Until then:

- `q_x` and `q_rank` (in `qx_scores.csv`) are **identical for every row** (verified:
  0 of 2,406 rows differ). `q_x` is the active per-gene weight the multiplier uses;
  it is the rank-linear transform, i.e. equal to `q_rank`. `q_x` is also duplicated
  into `multiplier_confidence.csv` with the same value.
- `target` (in `ground_truth_cellcycle_edges.csv`) means the same thing as
  `substrate` everywhere else: the regulated gene. Only this one file uses `target`.
- `source` (singular, in `augmented_network_candidate.csv` and
  `augmented_alpha_fitted.csv`) and `sources` (plural, in
  `union_network_alpha_fitted.csv`) carry the same kind of provenance information.
  The singular column holds one origin label per edge; the plural column holds a
  `|`-delimited list of origin labels for a unioned edge.

#### `tf_means.csv` (110 rows, one per measured TF) — the required multiplier input

| Column | Meaning |
| --- | --- |
| `tf` | transcription factor (gene name) |
| `tf_mean_rpm` | cell-cycle mean expression of the TF in reads per million (RPM), averaged on the NaN-interpolated 22-point grid so the multiplier is mean-preserving to machine precision |
| `tf_mean_molecules_per_cell` | the same mean expressed as molecules per cell, equal to `tf_mean_rpm * 60000 / 1e6` (i.e. `tf_mean_rpm * 0.06`); a convenience unit only, since the multiplier uses `TF(t)/TF_mean` and units cancel |
| `n_timepoints_available` | number of cell-cycle timepoints with a measured (non-NaN) value for this TF before interpolation |
| `n_timepoints_total` | total number of timepoints on the grid (22 for every TF); `available < total` indicates NaN-interpolated points |

#### `tf_network_fitted.csv` (3,929 rows), `joint_multidataset_alpha.csv` (7,890), `loo_kelliher_alpha.csv` (7,503)

| Column | Meaning |
| --- | --- |
| `substrate` | regulated gene |
| `tf` | regulating transcription factor |
| `alpha` | fitted signed coupling. In `tf_network_fitted.csv` only non-zero edges are listed; the joint and leave-Kelliher-out sets list candidate edges and may carry `alpha = 0.000000` for sparsified edges |

#### `qx_scores.csv` (2,406 rows, one per substrate gene)

| Column | Meaning |
| --- | --- |
| `substrate` | gene |
| `cyclebase_rank` | periodicity rank from Cyclebase 3.0 (lower = more periodic) |
| `periodicity_pvalue` | Cyclebase periodicity p-value (smaller = more periodic) |
| `peaktime_pct` | peak-expression phase as a percentage of the cell cycle (0-99) |
| `q_x` | per-gene cell-cycle weight in [0,1] actually used by the multiplier; identical to `q_rank` (see naming carryovers above) |
| `q_rank` | rank-linear [0,1] transform of periodicity; the default and equal to `q_x` |
| `q_sigmoid` | sigmoid [0,1] transform of periodicity; a documented alternative weight |
| `q_pcomplement` | p-value-complement [0,1] transform of periodicity; the third alternative weight |
| `multiplier_amplitude` | peak-to-trough amplitude of the gene's multiplier `M_x(t)` at `q_x = 1`, a descriptive size of the modulation |

Note: 15 genes absent from Cyclebase have empty `q_x`, `q_rank`, `q_sigmoid`,
`q_pcomplement`, `cyclebase_rank`, `periodicity_pvalue`, and `peaktime_pct`.

#### `tf_normalization_table.csv` (2,406 rows, one per substrate gene)

| Column | Meaning |
| --- | --- |
| `substrate` | gene |
| `n_tf` | number of TFs (fitted edges) regulating this gene |
| `n_activators` | number of those edges with `alpha > 0` |
| `n_repressors` | number of those edges with `alpha < 0` |
| `sum_alpha` | signed sum of the gene's couplings (`sum_i alpha_i`) |
| `sum_abs_alpha` | sum of absolute couplings (`sum_i |alpha_i|`), the deviation-form denominator |
| `baseline_ratio_signed_over_abs` | `sum_alpha / sum_abs_alpha`; equals 1.0 when the gene has no repressor |
| `flag_has_repressor` | 1 if the gene carries at least one repressor edge, else 0 |
| `flag_net_negative` | 1 if `sum_alpha < 0` (repressor-dominated), else 0 |
| `flag_near_zero` | 1 if `sum_alpha` is near zero (signed denominator would be unstable), else 0 |
| `min_mult_signed` | minimum over the cell cycle of the older signed-denominator multiplier (can be < 0 for repressor genes) |
| `min_mult_absdev` | minimum over the cell cycle of the deviation-form multiplier (the shipped form; stays non-negative) |
| `ever_negative_signed` | 1 if `min_mult_signed < 0` (the signed form drove this gene negative), else 0 |

#### `multiplier_examples.csv` (352 rows: 16 example genes x 22 timepoints)

| Column | Meaning |
| --- | --- |
| `substrate` | example gene |
| `time_min` | cell-cycle time in minutes (0, 5, ..., 150) |
| `mult_signed` | the older signed-denominator multiplier `M_x(t)` at this time; **may be negative** (53 of 352 example rows are), which is exactly the failure mode the deviation form fixes |
| `mult_absdev` | the shipped deviation-form multiplier `M_x(t)` at this time; stays non-negative |

#### `multisource_support.csv` (3,929 rows, one per fitted edge)

| Column | Meaning |
| --- | --- |
| `tf` | regulating transcription factor |
| `substrate` | regulated gene |
| `alpha` | fitted signed coupling for the edge |
| `abs_alpha` | `|alpha|` |
| `n_independent_sources` | number of the six independent networks (below) that contain this edge; equals the sum of the six non-`teufel` `in_*` flags (verified for all rows) |
| `in_teufel` | 1 if the edge is in the Teufel network (the fit's own source); always 1 here |
| `in_yeastract_ito` | 1 if the edge is in the YEASTRACT (Ito) network |
| `in_yeastract2019` | 1 if the edge is in the YEASTRACT 2019 network |
| `in_macisaac_chip` | 1 if the edge is in the MacIsaac ChIP-binding network |
| `in_inferelator_gold` | 1 if the edge is in the Inferelator gold standard |
| `in_snap` | 1 if the edge is in the SNAP network (0 for every row in this shipped file) |
| `in_curated` | 1 if the edge is in the curated benchmark |
| `sign_sources_n` | number of independent sources that provide a sign for this edge (0-3) |
| `sign_sources_agree` | number of those sources whose sign agrees with the fitted `alpha` sign (0-2) |

The "six independent networks" are the non-`teufel` `in_*` columns
(`in_yeastract_ito`, `in_yeastract2019`, `in_macisaac_chip`, `in_inferelator_gold`,
`in_snap`, `in_curated`); `in_teufel` is the seventh and is the fit's own source.

#### `augmented_network_candidate.csv` (7,186 rows) — optional v2 candidate edges

| Column | Meaning |
| --- | --- |
| `tf` | regulating transcription factor |
| `substrate` | regulated gene |
| `source` | origin of the candidate edge: one of `teufel`, `yeastract_strict`, `literature` |
| `yeastract_sign` | sign reported by YEASTRACT for the edge: `+`, `-`, `both`, or `NA` (no YEASTRACT sign) |
| `in_teufel` | 1 if the edge is also in the Teufel network, else 0 |

#### `augmented_alpha_fitted.csv` (4,802 rows) — alpha refit on the augmented network

| Column | Meaning |
| --- | --- |
| `substrate` | regulated gene |
| `tf` | regulating transcription factor |
| `alpha` | coupling refit on the augmented network |
| `source` | origin of the edge: `teufel` or `yeastract_strict` (singular column; see naming carryovers) |

#### `union_network_alpha_fitted.csv` (5,999 rows) — alpha refit on the union of all sources

| Column | Meaning |
| --- | --- |
| `substrate` | regulated gene |
| `tf` | regulating transcription factor |
| `alpha` | coupling refit on the union network |
| `sources` | `|`-delimited list of every origin that contributed the edge (e.g. `teufel|macisaac_chip|curated`); plural column (see naming carryovers) |

#### `ground_truth_cellcycle_edges.csv` (62 rows) — curated, cited textbook edges

| Column | Meaning |
| --- | --- |
| `tf` | regulating transcription factor |
| `target` | regulated gene (same role as `substrate` elsewhere; see naming carryovers) |
| `sign` | `activator` or `repressor` |
| `phase` | cell-cycle phase of the regulation: `G1`, `G1/S`, `G2/M`, or `M/G1` |
| `citation` | `;`-delimited literature reference(s) supporting the edge (e.g. `Iyer2001;Simon2001`) |

#### `multiplier_confidence.csv` (2,406 rows, one per substrate gene)

| Column | Meaning |
| --- | --- |
| `substrate` | gene |
| `q_x` | per-gene cell-cycle weight, duplicated from `qx_scores.csv` (same value) |
| `boot_robust_frac` | fraction of bootstrap resamples in which the gene's edges are robust (0-1) |
| `mean_n_sources` | mean number of independent sources supporting the gene's edges (0-4 observed) |
| `transfer_r_median` | median cross-dataset transfer correlation for the gene (can be negative; -0.78 to 0.95 observed) |
| `repressor_score` | repressor-fragility component (0-1); higher is more confident |
| `carries_repressor` | `true`/`false`: whether the gene has any repressor edge |
| `n_repressor_edges` | number of repressor edges for the gene |
| `n_edges` | total number of fitted edges for the gene |
| `n_robust_edges` | number of edges judged robust (e.g. under bootstrap) |
| `periodicity_pvalue` | Cyclebase periodicity p-value (as in `qx_scores.csv`) |
| `score` | overall 0-1 trust score, a heuristic blend of the components above (not a calibrated probability) |
| `tier` | `USE`, `CAUTION`, or `DO-NOT-USE`, derived from `score` and the reasons |
| `reason` | `;`-delimited reason codes when not a clean USE (e.g. `non-periodic`, `over-parameterized`, `sign-unstable`, `low-overall-score`); empty for clean genes |
| `apply_multiplier` | boolean drop-in gate: `true` (1,058 genes) means apply `M_x(t)`; `false` (1,348) means fall back to a flat rate (`M = 1`). `true` iff `tier == USE` and the gene is periodic. The CAUTION tier shows no held-out transfer advantage over flat, so the practical gate is binary |
| `recommended_fallback` | `none` for apply genes, `flat_M1` for DO-NOT-USE, `damped` for CAUTION (damping toward 1 is an optional engineering choice, not data-supported) |

> The exact weighting formula that combines these components into `score`, and the
> thresholds that map `score` to `tier`, live in the generating code, not in this
> table. Treat `score`/`tier` as a heuristic downstream flag.

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
  TF's non-NaN timepoints, which left up to a ~11% offset on NaN-affected genes).
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

MIT (see `LICENSE`). See `CITATION.cff` for how to cite this work. Please also
cite the source data (Teufel et al. 2019, Sci Rep 9:3343) and link back to
`https://github.com/djsegal/TranscriptionMultiplier.jl`. `CHECKSUMS.sha256` lists the SHA-256
of every released `data/` CSV.
