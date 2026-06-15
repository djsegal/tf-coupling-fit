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

const DATA_DIR = joinpath(@__DIR__, "data")

"Load (edges::Dict{String,Vector{Tuple{String,Float64}}}, means::Dict{String,Float64})."
function load_handoff(dir::AbstractString=DATA_DIR)
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

`tf_levels::Dict{String,<:Real}` current TF amounts. `q ∈ [0,1]` scales cycling
amplitude; `clamp` applies a max(0, M) floor.
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
