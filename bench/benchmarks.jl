#!/usr/bin/env julia
# EE1: BenchmarkTools suite for the transcription-multiplier hot paths.
#
# Run from the repo root:
#
#   julia --project=bench bench/benchmarks.jl          # print results
#   julia --project=bench bench/benchmarks.jl --write  # also (re)write BENCHMARKS.md
#
# ADDITIVE: this lives in its own environment (bench/Project.toml) so the
# package's own Project.toml / test stack stay untouched. It loads the package
# by adding the parent directory's project to LOAD_PATH at runtime.
#
# Times the three hot paths:
#   1. vector-form  multiplier(::Vector, ::Vector; q)        (the Mseries inner loop)
#   2. Dict-form    multiplier(substrate, levels, edges, means; q)
#   3. single-substrate L1 fit  fit_substrate(sub_series, tf_matrix, w2)

using BenchmarkTools
using Printf

# --- load the package from the sibling project (../) ------------------------
const PKG_DIR = normpath(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, PKG_DIR)
try
    @eval using TranscriptionMultiplier
    @eval const M             = TranscriptionMultiplier.multiplier
    @eval const FIT_SUBSTRATE = TranscriptionMultiplier.fit_substrate
catch err
    @warn "Could not load TranscriptionMultiplier as a module; including scripts directly." err
    include(joinpath(PKG_DIR, "multiplier.jl"))
    include(joinpath(PKG_DIR, "refit.jl"))
    @eval const M             = multiplier
    @eval const FIT_SUBSTRATE = fit_substrate
end

# --------------------------------------------------------------- inputs -----
# (1) vector-form multiplier: a representative regulator set + ratios.
const ALPHAS = [2.0, -1.0, 0.5, 1.25]
const RATIOS = [1.1, 0.9, 1.3, 0.8]

# (2) Dict-form multiplier: CLN2-like substrate with three regulators.
const EDGES  = Dict("CLN2" => [("SWI4", 2.0), ("MBP1", -1.0), ("SWI6", 0.5)])
const MEANS  = Dict("SWI4" => 100.0, "MBP1" => 50.0, "SWI6" => 75.0)
const LEVELS = Dict("SWI4" => 120.0, "MBP1" => 40.0,  "SWI6" => 80.0)

# (3) single-substrate L1 fit: a small but realistic k x T problem.
#     18 effective timepoints (22-point grid minus the tau lag) x 4 regulators.
const T_FIT  = 18
const K_FIT  = 4
const TF_MAT = let
    # deterministic, well-conditioned synthetic TF trajectories
    [sin(0.4 * j + 0.7 * t) + 1.5 for j in 1:K_FIT, t in 1:T_FIT]
end
const SUB_SERIES = let
    true_alpha = [1.0, -0.5, 0.25, 0.75]
    vec(true_alpha' * TF_MAT) .+ 0.01 .* cos.(1:T_FIT)
end
const W2 = 5.0

# --------------------------------------------------------------- suite ------
const SUITE = BenchmarkGroup()
SUITE["multiplier_vector"] = @benchmarkable M($ALPHAS, $RATIOS; q = 1.0)
SUITE["multiplier_dict"]   = @benchmarkable M("CLN2", $LEVELS, $EDGES, $MEANS; q = 1.0)
SUITE["fit_substrate"]     = @benchmarkable FIT_SUBSTRATE($SUB_SERIES, $TF_MAT, $W2)

println("Tuning + warming up benchmark suite ...")
tune!(SUITE)
results = run(SUITE; verbose = true)

# --------------------------------------------------------------- report -----
const ORDER = ["multiplier_vector", "multiplier_dict", "fit_substrate"]
const LABEL = Dict(
    "multiplier_vector" => "vector multiplier (`multiplier(::Vector, ::Vector; q)`)",
    "multiplier_dict"   => "Dict multiplier (`multiplier(substrate, levels, edges, means; q)`)",
    "fit_substrate"     => "single-substrate L1 fit (`fit_substrate`)",
)

fmt_time(ns) = ns < 1e3  ? @sprintf("%.1f ns", ns) :
               ns < 1e6  ? @sprintf("%.3f us", ns / 1e3) :
               ns < 1e9  ? @sprintf("%.3f ms", ns / 1e6) :
                           @sprintf("%.3f s",  ns / 1e9)
fmt_bytes(b) = b < 1024 ? @sprintf("%d B", b) :
               b < 1024^2 ? @sprintf("%.2f KiB", b / 1024) :
                            @sprintf("%.2f MiB", b / 1024^2)

println("\n=== baseline benchmark results (median) ===")
rows = String[]
for key in ORDER
    t  = results[key]
    md = median(t)
    mn = minimum(t)
    @printf("  %-22s median %-12s  min %-12s  %-10s  %d allocs\n",
            key, fmt_time(time(md)), fmt_time(time(mn)),
            fmt_bytes(memory(md)), allocs(md))
    push!(rows, "| $(LABEL[key]) | $(fmt_time(time(md))) | $(fmt_time(time(mn))) | $(fmt_bytes(memory(md))) | $(allocs(md)) |")
end

# --------------------------------------------------------------- write ------
if "--write" in ARGS
    out = joinpath(PKG_DIR, "BENCHMARKS.md")
    ver = try string(VERSION) catch; "?" end
    body = """
    # Benchmarks

    Baseline median timings for the transcription-multiplier hot paths, captured
    with [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl). These
    are the numbers a future change is compared against, so regressions show up.

    Reproduce (standalone, from the package directory):

    ```
    julia --project=bench transcription_multiplier/bench/benchmarks.jl --write
    ```

    The suite times three hot paths:

    1. the vector-form `multiplier` call (the inner loop of the M-series / scoring),
    2. the Dict-form `multiplier` public call, and
    3. a single-substrate L1 fit via `fit_substrate` (JuMP/HiGHS).

    ## Baseline (median over the BenchmarkTools sample)

    Captured on Julia $ver.

    | hot path | median | min | memory | allocations |
    |---|---|---|---|---|
    $(join(rows, "\n"))

    Notes:

    - `multiplier_vector` / `multiplier_dict` are closed-form and allocation-light;
      they are the calls executed millions of times when scoring trajectories.
    - `fit_substrate` is dominated by JuMP model construction + the HiGHS solve, so
      its time/allocations are orders of magnitude larger and naturally noisier.
    - Times depend on the machine; treat the committed numbers as a *relative*
      regression guard rather than an absolute spec.
    """
    write(out, body)
    println("\nWrote baseline to ", out)
end
