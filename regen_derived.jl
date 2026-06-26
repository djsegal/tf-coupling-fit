#!/usr/bin/env julia
# Regenerate the public derived tables that depend on tf_means.csv so the package
# stays self-consistent after the TF means are recomputed on the interpolated grid:
#   data/multiplier_examples.csv     (mult_signed vs mult_absdev over the cell cycle)
#   data/tf_normalization_table.csv  (per-substrate alpha sums, flags, min multipliers)
# Reuses refit.jl::load_rna_seq for the exact interpolated 22-point TF trajectories.
#   julia --project=. regen_derived.jl [--validate]
using CSV, DataFrames, Statistics
include(joinpath(@__DIR__, "refit.jl"))

const REPO = @__DIR__
const EXPR = joinpath(REPO, "data", "WT_unstressed_readspermillionreads.csv")
const DATA = joinpath(@__DIR__, "data")

validate = "--validate" in ARGS
means_file = validate ? "HEAD:data/tf_means.csv" : joinpath(DATA, "tf_means.csv")

# Interpolated 22-point TF trajectories (RPM->molecules; ratios are scale-free).
time_axis, genes, expr = load_rna_seq(EXPR)
gidx = Dict(g => i for (i, g) in enumerate(genes))
T = length(time_axis)

# Edges
edf = CSV.read(joinpath(DATA, "tf_network_fitted.csv"), DataFrame)
edges = Dict{String,Vector{Tuple{String,Float64}}}()
for row in eachrow(edf)
    push!(get!(edges, String(row.substrate), Tuple{String,Float64}[]),
          (String(row.tf), Float64(row.alpha)))
end

# Means (from file, or from HEAD when validating the logic against committed tables)
mdf = validate ? CSV.read(IOBuffer(read(`git -C $REPO show $means_file`, String)), DataFrame) :
                 CSV.read(means_file, DataFrame)
means = Dict{String,Float64}()
for row in eachrow(mdf)
    ismissing(row.tf_mean_rpm) || (means[String(row.tf)] = Float64(row.tf_mean_rpm))
end

# Ratio r_i(t) = TF_i(t)/TF_i_mean on the interpolated grid (in RPM units to match means).
ratio(tf, k) = (expr[gidx[tf], k] / RPM_TO_MOLECULES) / means[tf]

# Form (I) signed and form (II) deviation multipliers for a substrate at timepoint k.
function mults(regs, k)
    sa  = sum(a for (_, a) in regs)
    sab = sum(abs(a) for (_, a) in regs)
    msigned = sa == 0 ? NaN : sum(a * ratio(tf, k) for (tf, a) in regs) / sa
    mabs    = sab == 0 ? 1.0 : 1.0 + sum(a * (ratio(tf, k) - 1.0) for (tf, a) in regs) / sab
    return msigned, mabs
end

# ---- multiplier_examples.csv (preserve substrate set + ordering from committed file) ----
ex_old = CSV.read(joinpath(DATA, "multiplier_examples.csv"), DataFrame)
ex_subs = unique(String.(ex_old.substrate))
ex_rows = NamedTuple[]
for sub in ex_subs
    regs = edges[sub]
    for k in 1:T
        ms, ma = mults(regs, k)
        push!(ex_rows, (substrate = sub, time_min = Int(time_axis[k]),
                        mult_signed = round(ms, digits = 5),
                        mult_absdev = round(ma, digits = 5)))
    end
end
ex_new = DataFrame(ex_rows)

# ---- tf_normalization_table.csv (all substrates, alphabetical as committed) ----
nt_old = CSV.read(joinpath(DATA, "tf_normalization_table.csv"), DataFrame)
nt_rows = NamedTuple[]
for sub in String.(nt_old.substrate)
    regs = edges[sub]
    nact = count(a > 0 for (_, a) in regs)
    nrep = count(a < 0 for (_, a) in regs)
    sa  = sum(a for (_, a) in regs)
    sab = sum(abs(a) for (_, a) in regs)
    baseline = sab == 0 ? 1.0 : sa / sab
    msig = Inf; mabs = Inf
    for k in 1:T
        ms, ma = mults(regs, k)
        msig = min(msig, ms); mabs = min(mabs, ma)
    end
    push!(nt_rows, (substrate = sub, n_tf = length(regs),
        n_activators = nact, n_repressors = nrep,
        sum_alpha = round(sa, digits = 6), sum_abs_alpha = round(sab, digits = 6),
        baseline_ratio_signed_over_abs = round(baseline, digits = 6),
        flag_has_repressor = nrep > 0 ? 1 : 0,
        flag_net_negative = sa < 0 ? 1 : 0,
        flag_near_zero = abs(sa) < 0.05 ? 1 : 0,
        min_mult_signed = round(msig, digits = 4),
        min_mult_absdev = round(mabs, digits = 4),
        ever_negative_signed = msig < 0 ? 1 : 0))
end
nt_new = DataFrame(nt_rows)

if validate
    # Compare against committed files (computed with HEAD means) to prove the logic.
    function maxdiff(a::DataFrame, b::DataFrame, cols)
        m = 0.0
        for c in cols, i in 1:nrow(a)
            av = a[i, c]; bv = b[i, c]
            (ismissing(av) || ismissing(bv)) && continue
            (av isa Real && isnan(av)) && continue
            m = max(m, abs(Float64(av) - Float64(bv)))
        end
        m
    end
    println("VALIDATE multiplier_examples: max|Δ| = ",
            maxdiff(ex_new, ex_old, [:mult_signed, :mult_absdev]))
    println("VALIDATE norm table sums/min: max|Δ| = ",
            maxdiff(nt_new, nt_old, [:sum_alpha, :sum_abs_alpha,
                :baseline_ratio_signed_over_abs, :min_mult_signed, :min_mult_absdev]))
    flagcols = [:n_tf, :n_activators, :n_repressors, :flag_has_repressor,
                :flag_net_negative, :flag_near_zero, :ever_negative_signed]
    println("VALIDATE norm table flags:    max|Δ| = ", maxdiff(nt_new, nt_old, flagcols))
else
    CSV.write(joinpath(DATA, "multiplier_examples.csv"), ex_new)
    CSV.write(joinpath(DATA, "tf_normalization_table.csv"), nt_new)
    println("Rewrote multiplier_examples.csv ($(nrow(ex_new)) rows) and ",
            "tf_normalization_table.csv ($(nrow(nt_new)) rows).")
end
