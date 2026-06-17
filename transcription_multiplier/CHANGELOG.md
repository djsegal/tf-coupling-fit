# Changelog

All notable changes to the `TranscriptionMultiplier` package are documented here.
This project loosely follows [Keep a Changelog](https://keepachangelog.com/) and
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **Module-ization.** `src/TranscriptionMultiplier.jl` wraps the existing
  scripts as a proper Julia module, so `using TranscriptionMultiplier` exposes
  the core API: the mean-preserving `multiplier` (both the low-level
  `multiplier(alphas::Vector, ratios::Vector; q)` method and the convenience
  `multiplier(substrate, tf_levels, edges, means; q, clamp)` method), and the
  load/refit helpers `load_handoff`, `load_rna_seq`, `load_tf_network`,
  `build_regulators`, `fit_substrate`, `fit_all`, `linear_interp`. This is
  purely additive: the module `include`s the original `multiplier.jl` and
  `refit.jl` at their existing paths (preserving each file's `@__DIR__` so the
  relative `data/` lookups still resolve). The standalone scripts and
  `runtests.jl` continue to work unchanged.
- **Package metadata.** A package-level `Project.toml` (`name`, `uuid`,
  `[deps]`, `[compat]`) so the package can be activated with
  `julia --project=.` from `transcription_multiplier/`.
- **Property/invariant tests.** `proptests.jl` exercises the multiplier's
  stated invariants on randomized inputs: (a) `<M>_t == 1` to machine precision
  for any alpha signs and any `q`; (b) `M >= 0` over the grid (clamped method);
  (c) monotonicity of deviation amplitude `|M-1|` in `q`; (d) edge cases
  (all-zero alpha, single edge, repressor-only). Independent of `runtests.jl`.
- **Performance spot check.** `perf_check.jl` runs `@code_warntype` and
  `@allocated`/`@time` on the multiplier hot path. Result: both methods are
  type-stable (`Body::Float64`); the vector method is non-allocating (0
  bytes/call) and the Dict method allocates 16 bytes/call.

### Unchanged

- No existing script was moved, renamed, or behaviorally altered. The notebook
  entry points, the committed handoff, and `runtests.jl` (35 tests) all still
  run and pass exactly as before.
