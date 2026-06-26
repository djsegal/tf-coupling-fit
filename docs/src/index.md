```@meta
CurrentModule = TranscriptionMultiplier
```

# TranscriptionMultiplier.jl

A mean-preserving cell-cycle transcription-rate multiplier for dynamic models of
the yeast cell cycle, plus the L1 linear-program fit that produces its
coefficients.

The gated transcription rate of a substrate gene `x` is

```math
k_x(t) = k_{x,\text{base}} \cdot M_x(t), \qquad
M_x(t) = 1 + q_x \frac{\sum_i \alpha_i \left(\frac{\mathrm{TF}_i(t)}{\overline{\mathrm{TF}_i}} - 1\right)}{\sum_i |\alpha_i|}.
```

The deviation form is **mean-preserving**: ``\langle M_x \rangle_t = 1`` for any
sign pattern of the ``\alpha_i`` and any ``q_x``, it never blows up, and units
cancel in ``\mathrm{TF}_i(t)/\overline{\mathrm{TF}_i}``.

```@docs
TranscriptionMultiplier
```

## Quick start

```julia
using TranscriptionMultiplier

edges, means = load_handoff()                       # fitted (tf, alpha) + TF means
M = multiplier("CLN2", Dict("SWI4" => 120.0, "STE12" => 50.0), edges, means; q = 1.0)
```

There is also a low-level vector method and a small CLI
(`bin/multiplier_cli.jl`) that prints ``M_x`` and the gated rate for a gene.

## Mean preservation

The headline invariant, as a runnable example:

```@example
using TranscriptionMultiplier, Statistics
alphas = [2.0, -1.0, 0.5]
levels = [1.0 3.0 2.0; 4.0 1.0 1.0; 2.0 2.0 5.0]   # 3 TFs x 3 timepoints
mus    = vec(mean(levels; dims = 2))
Ms     = [multiplier(alphas, levels[:, t] ./ mus; q = 0.7) for t in 1:3]
mean(Ms)                                            # == 1 to machine precision
```

## The multiplier

```@docs
multiplier
```

## Loading the fitted handoff

```@docs
load_handoff
DATA_DIR
```

## Fitting (L1 linear program)

```@docs
load_rna_seq
load_tf_network
build_regulators
fit_substrate
fit_all
linear_interp
```

## Constants

```@docs
MOLECULES_PER_CELL
RPM_TO_MOLECULES
DEFAULT_TAU
DEFAULT_W2
DEFAULT_MAX_NANS
```

## Index

```@index
```
