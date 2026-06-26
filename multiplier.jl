#=
Drop-in (Julia): cell-cycle transcription-rate multiplier for a dynamic model.

    k_x(t) = k_x_base * Mₓ(t)
    Mₓ(t)  = 1 + qₓ * ( Σᵢ αᵢ (TFᵢ(t)/TFᵢ_mean − 1) ) / Σᵢ |αᵢ|

Mean-preserving deviation form (see report/report.pdf): ⟨M⟩ₜ = 1 for any α signs
and any qₓ, never blows up, keeps repressor signs, identical to the signed form
when all αᵢ > 0. Units cancel in TFᵢ(t)/TFᵢ_mean.

The TFᵢ_mean values in data/tf_means.csv are the discrete means of the same
NaN-interpolated 22-point trajectory the multiplier is evaluated on, so ⟨M⟩ₜ = 1
holds to machine precision for every gene (not just genes with complete TF data).
=#
using CSV, DataFrames

"""
    DATA_DIR

Absolute path to the package's bundled `data/` directory (fitted network, TF
means, q_x scores). Default argument to [`load_handoff`](@ref).
"""
const DATA_DIR = joinpath(@__DIR__, "data")

"""
    load_handoff(dir::AbstractString=DATA_DIR) -> (edges, means)

Load the fitted handoff produced by the fit pipeline:

  - `edges::Dict{String,Vector{Tuple{String,Float64}}}` maps each substrate gene to
    its list of `(tf, alpha)` regulator coefficients (`data/tf_network_fitted.csv`).
  - `means::Dict{String,Float64}` maps each TF to its mean RPM level over the
    NaN-interpolated cell-cycle trajectory (`data/tf_means.csv`).

These are exactly the inputs the [`multiplier`](@ref) Dict method expects.

Throws an informative `ArgumentError` if `dir` does not exist or is missing
either CSV (instead of a cryptic file-not-found from the CSV reader).

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> edges, means = load_handoff();

julia> edges isa Dict{String,Vector{Tuple{String,Float64}}}
true

julia> haskey(edges, "CLN2") && haskey(means, first(first(edges["CLN2"])))
true
```
"""
function load_handoff(dir::AbstractString=DATA_DIR)
    isdir(dir) ||
        throw(ArgumentError("load_handoff: data directory does not exist: $dir"))
    for f in ("tf_network_fitted.csv", "tf_means.csv")
        isfile(joinpath(dir, f)) ||
            throw(ArgumentError("load_handoff: missing required file '$f' in $dir"))
    end
    e = CSV.read(joinpath(dir, "tf_network_fitted.csv"), DataFrame)
    edges = Dict{String,Vector{Tuple{String,Float64}}}()
    for row in eachrow(e)
        push!(get!(edges, row.substrate, Tuple{String,Float64}[]),
              (row.tf, Float64(row.alpha)))
    end
    m = CSV.read(joinpath(dir, "tf_means.csv"), DataFrame)
    means = Dict{String,Float64}()
    for row in eachrow(m)
        ismissing(row.tf_mean_rpm) || (means[row.tf] = Float64(row.tf_mean_rpm))
    end
    return edges, means
end

"""
    multiplier(substrate, tf_levels, edges, means; q=1.0, clamp=true) -> Float64

Mean-preserving cell-cycle transcription-rate multiplier `Mₓ(t)` for one
`substrate` gene, computed from a Dict of current TF amounts:

    Mₓ(t) = 1 + q * ( Σᵢ αᵢ (TFᵢ(t)/TFᵢ_mean − 1) ) / Σᵢ |αᵢ|

Arguments:

  - `substrate`            : gene name (key into `edges`).
  - `tf_levels::Dict`      : current TF amounts, `Dict{String,<:Real}`.
  - `edges::Dict`          : `substrate => [(tf, alpha), ...]` (see [`load_handoff`](@ref)).
  - `means::Dict`          : `tf => mean_level` over the cell cycle.

Keywords: `q` scales the cycling amplitude (the report uses `q ∈ [0,1]`, with
`q=0` giving a flat `M≡1`); `clamp` applies a `max(0, M)` floor.

A substrate with no regulators (or all-zero `alpha`) returns `1.0`. TFs absent
from `tf_levels`/`means` (or with zero mean) are skipped, so partial input is
tolerated rather than erroring.

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> edges = Dict("CLN2" => [("SWI4", 2.0), ("MBP1", -1.0)]);

julia> means = Dict("SWI4" => 100.0, "MBP1" => 50.0);

julia> multiplier("CLN2", Dict("SWI4" => 100.0, "MBP1" => 50.0), edges, means)  # at the mean
1.0

julia> multiplier("UNKNOWN_GENE", Dict("SWI4" => 999.0), edges, means)  # no regulators
1.0
```
"""
function multiplier(substrate, tf_levels, edges, means; q::Float64=1.0, clamp::Bool=true)
    haskey(edges, substrate) || return 1.0
    regs = edges[substrate]
    sabs = sum(abs(a) for (_, a) in regs)
    sabs == 0 && return 1.0
    dev = 0.0
    for (tf, a) in regs
        (haskey(means, tf) && means[tf] != 0 && haskey(tf_levels, tf)) || continue
        dev += a * (tf_levels[tf] / means[tf] - 1.0)
    end
    M = 1.0 + q * dev / sabs
    return clamp ? max(0.0, M) : M
end

"""
    save_handoff(dir, edges, means) -> dir

Inverse of [`load_handoff`](@ref): write `edges` and `means` back out to
`tf_network_fitted.csv` and `tf_means.csv` under `dir` (created if needed), in
the exact schema `load_handoff` reads. Round-trips: after
`save_handoff(dir, e, m)`, `load_handoff(dir)` returns dicts equal to `(e, m)`.
Rows are sorted (substrate, then TF order as given; means by TF) for a
deterministic, diff-friendly file. Uses only the package's existing `CSV`
dependency — no extra serialization libraries.
"""
function save_handoff(dir::AbstractString, edges::AbstractDict, means::AbstractDict)
    isdir(dir) || mkpath(dir)
    esub = String[]; etf = String[]; ealpha = Float64[]
    for sub in sort(collect(keys(edges))), (tf, a) in edges[sub]
        push!(esub, sub); push!(etf, tf); push!(ealpha, Float64(a))
    end
    CSV.write(joinpath(dir, "tf_network_fitted.csv"),
              DataFrame(substrate = esub, tf = etf, alpha = ealpha))
    mtf = String[]; mval = Float64[]
    for tf in sort(collect(keys(means)))
        push!(mtf, tf); push!(mval, Float64(means[tf]))
    end
    CSV.write(joinpath(dir, "tf_means.csv"),
              DataFrame(tf = mtf, tf_mean_rpm = mval))
    return dir
end

"""
    export_handoff_json(path, edges, means) -> path

Write the handoff as one portable JSON document for non-Julia consumers (e.g. the
whole-cell-model team's Python stack):

    {"edges": {substrate: [[tf, alpha], ...], ...}, "means": {tf: mean, ...}}

Hand-written and write-only, so the package gains no JSON dependency (the board's
"optional, not a heavy hard dep" rule); read it back with any JSON library. Keys
are plain gene/TF symbols, but strings are escaped defensively. `alpha`/`mean`
are finite by construction, so the numeric output is always valid JSON.
"""
function export_handoff_json(path::AbstractString, edges::AbstractDict, means::AbstractDict)
    esc(s) = replace(string(s), "\\" => "\\\\", "\"" => "\\\"")
    io = IOBuffer()
    print(io, "{\n  \"edges\": {")
    subs = sort(collect(keys(edges)))
    for (i, sub) in enumerate(subs)
        print(io, i == 1 ? "\n" : ",\n", "    \"", esc(sub), "\": [")
        for (j, (tf, a)) in enumerate(edges[sub])
            print(io, j == 1 ? "" : ", ", "[\"", esc(tf), "\", ", Float64(a), "]")
        end
        print(io, "]")
    end
    print(io, "\n  },\n  \"means\": {")
    tfs = sort(collect(keys(means)))
    for (i, tf) in enumerate(tfs)
        print(io, i == 1 ? "\n" : ",\n", "    \"", esc(tf), "\": ", Float64(means[tf]))
    end
    print(io, "\n  }\n}\n")
    write(path, take!(io))
    return path
end
