"""
    TranscriptionMultiplier

Module wrapper for the yeast cell-cycle transcription-rate multiplier package.

This is an *additive* wrapper: it does not move or rewrite any of the existing
top-level scripts (`multiplier.jl`, `refit.jl`, `joint_fit.jl`, `joint_score.jl`,
`runtests.jl`, `regen_*.jl`). Those remain usable standalone exactly as before.
The module simply `include`s the canonical implementations and re-exports the
core public API so that `using TranscriptionMultiplier` works.

Included (with their original `@__DIR__`, so relative `data/` paths still resolve):

  - `../multiplier.jl` : `load_handoff`, `multiplier(substrate, tf_levels, edges, means; ...)`
  - `../refit.jl`      : `load_rna_seq`, `load_tf_network`, `build_regulators`,
                         `fit_substrate`, `fit_all`, `linear_interp`,
                         `multiplier(alphas::Vector, ratios::Vector; q=...)`

Both `multiplier.jl` and `refit.jl` add methods to a single `multiplier`
generic function (a Dict-based convenience method and a low-level vector method);
both are re-exported.
"""
module TranscriptionMultiplier

# The canonical implementations live as scripts one directory up. Including them
# by their real paths keeps each file's `@__DIR__` pointing at the package root,
# so `const DATA_DIR = joinpath(@__DIR__, "data")` etc. continue to resolve.
include(joinpath(@__DIR__, "..", "multiplier.jl"))
include(joinpath(@__DIR__, "..", "refit.jl"))

# ---- public API -----------------------------------------------------------
# mean-preserving multiplier (both Dict and Vector methods)
export multiplier
# load / refit helpers
export load_handoff, load_rna_seq, load_tf_network, build_regulators,
       fit_substrate, fit_all, linear_interp
# default data directory used by load_handoff
export DATA_DIR

# Useful constants from the fit pipeline (re-exported for downstream callers).
export MOLECULES_PER_CELL, RPM_TO_MOLECULES, DEFAULT_TAU, DEFAULT_W2, DEFAULT_MAX_NANS

end # module TranscriptionMultiplier
