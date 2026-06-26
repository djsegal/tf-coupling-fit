#=
refit.jl: Julia reproduction of the tf-coupling-fit linear program (JuMP/HiGHS),
plus leave-one-timepoint-out cross-validation. Same solver stack as the generating
notebook, so it reproduces the committed alpha exactly. Used by runtests.jl.

    include("refit.jl")
    ta, genes, expr = load_rna_seq(joinpath(@__DIR__,"data","WT_unstressed_readspermillionreads.csv"))
    edges = load_tf_network(joinpath(@__DIR__,"data","TF_von_Teufel.csv"))
    reg   = build_regulators(edges, Set(genes))
    res   = fit_all(ta, Dict(g=>i for (i,g) in enumerate(genes)), expr, reg)
=#
using CSV, DataFrames, JuMP, HiGHS, Statistics

"Assumed mRNA molecules per cell at the genome-wide mean level (RPM→molecules scale)."
const MOLECULES_PER_CELL = 60_000
"RPM→molecules-per-cell conversion factor (`MOLECULES_PER_CELL / 1e6` = 0.06)."
const RPM_TO_MOLECULES   = MOLECULES_PER_CELL / 1_000_000   # 0.06
"Default TF→substrate regulatory lag (minutes) used by [`fit_all`](@ref)."
const DEFAULT_TAU      = 20.0
"Default L1 sparsity penalty on the fitted coefficients ([`fit_substrate`](@ref)/[`fit_all`](@ref))."
const DEFAULT_W2       = 5.0
"Default max number of missing timepoints a gene may have before it is dropped in [`load_rna_seq`](@ref)."
const DEFAULT_MAX_NANS = 4

"""
    linear_interp(x, xp, fp) -> Float64

Piecewise-linear interpolation of the points `(xp, fp)` at `x`, clamped to the
endpoints outside `[xp[1], xp[end]]`. `xp` must be sorted ascending. Used to
evaluate TF trajectories (and lag them by `tau`) on the cell-cycle grid.

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> xp = [0.0, 10.0, 20.0]; fp = [1.0, 3.0, 2.0];

julia> linear_interp(5, xp, fp)    # halfway between 1 and 3
2.0

julia> linear_interp(-5, xp, fp)   # below range -> first value
1.0

julia> linear_interp(25, xp, fp)   # above range -> last value
2.0
```
"""
function linear_interp(x::Real, xp::AbstractVector, fp::AbstractVector)
    x = Float64(x)
    if x <= xp[1]
        return Float64(fp[1])
    elseif x >= xp[end]
        return Float64(fp[end])
    end
    i = searchsortedlast(xp, x)
    x0, x1 = Float64(xp[i]), Float64(xp[i+1])
    f0, f1 = Float64(fp[i]), Float64(fp[i+1])
    return f0 + (f1 - f0) * (x - x0) / (x1 - x0)
end

"""
    load_rna_seq(path; max_nans=DEFAULT_MAX_NANS) -> (time_axis, gene_names, expression)

Read a wide RNA-seq CSV (a `name` column plus one numeric column per timepoint).
Numeric-named columns are detected and sorted into `time_axis`; genes with more
than `max_nans` missing timepoints are dropped, remaining gaps are filled by
[`linear_interp`](@ref), and levels are rescaled from RPM to molecules-per-cell
(`RPM_TO_MOLECULES`).

Returns `(time_axis::Vector{Float64}, gene_names::Vector{String},
expression::Matrix{Float64})` where `expression` is `genes × timepoints`.

Throws an `ArgumentError` if `path` does not exist.

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> ta, genes, expr = load_rna_seq(joinpath(DATA_DIR, "..", "..", "data",
                                                "WT_unstressed_readspermillionreads.csv"));

julia> size(expr, 1) == length(genes) && size(expr, 2) == length(ta)
true

julia> all(isfinite, expr)
true
```
"""
function load_rna_seq(path::String; max_nans::Int = DEFAULT_MAX_NANS)
    isfile(path) || throw(ArgumentError("load_rna_seq: file not found: $path"))
    df = CSV.read(path, DataFrame)
    time_cols = String[]; time_axis_raw = Float64[]
    for c in names(df)
        v = tryparse(Float64, string(c))
        if v !== nothing
            push!(time_cols, string(c)); push!(time_axis_raw, v)
        end
    end
    ord = sortperm(time_axis_raw)
    time_cols = time_cols[ord]; time_axis = time_axis_raw[ord]
    is_missing(v) = ismissing(v) || (v isa Number && isnan(v))
    keep = trues(nrow(df))
    for i in 1:nrow(df)
        nan_count = 0
        for c in time_cols
            is_missing(df[i, c]) && (nan_count += 1)
        end
        nan_count > max_nans && (keep[i] = false)
    end
    df = df[keep, :]
    df[!, :name] = String[ismissing(x) ? "" : String(strip(string(x))) for x in df.name]
    df = filter(row -> !isempty(row.name) && lowercase(row.name) != "nan", df)
    gene_names = Vector{String}(df.name)
    T = length(time_cols); G = nrow(df)
    expression = Matrix{Float64}(undef, G, T)
    for j in 1:T
        col = df[!, time_cols[j]]
        for i in 1:G
            v = col[i]; expression[i, j] = is_missing(v) ? NaN : Float64(v)
        end
    end
    keep_after = trues(G)
    for i in 1:G
        row = view(expression, i, :); nan_mask = isnan.(row)
        any(nan_mask) || continue
        good = .!nan_mask
        if !any(good); keep_after[i] = false; continue; end
        x_good = time_axis[good]; y_good = Float64.(expression[i, good])
        for k in findall(nan_mask)
            expression[i, k] = linear_interp(time_axis[k], x_good, y_good)
        end
    end
    expression = expression[keep_after, :]; gene_names = gene_names[keep_after]
    expression .*= RPM_TO_MOLECULES
    return time_axis, gene_names, expression
end

"""
    load_tf_network(csv_path) -> Vector{Tuple{String,String}}

Read a TF→substrate edge list (TF in column 1, substrate in column 3), skipping
blank or purely numeric cells. Returns the `(tf, substrate)` pairs as found
(use [`build_regulators`](@ref) to restrict to measured genes and deduplicate).

Throws an `ArgumentError` if `csv_path` does not exist.
"""
function load_tf_network(csv_path::String)
    isfile(csv_path) || throw(ArgumentError("load_tf_network: file not found: $csv_path"))
    df = CSV.read(csv_path, DataFrame)
    edges = Tuple{String,String}[]
    for row in eachrow(df)
        tf_val = row[1]; sub_val = row[3]
        (ismissing(tf_val) || ismissing(sub_val)) && continue
        tf = strip(string(tf_val)); sub = strip(string(sub_val))
        (isempty(tf) || isempty(sub)) && continue
        tryparse(Float64, tf)  !== nothing && continue
        tryparse(Float64, sub) !== nothing && continue
        push!(edges, (String(tf), String(sub)))
    end
    return edges
end

"""
    build_regulators(edges, measured) -> Dict{String,Vector{String}}

From a list of `(tf, substrate)` `edges`, keep only edges where BOTH endpoints
are in the `measured` set, then map each substrate to its sorted, deduplicated
list of regulator TFs.

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> edges = [("A","X"), ("B","X"), ("A","X"), ("C","Y")];

julia> reg = build_regulators(edges, Set(["A","B","X"]));

julia> reg["X"]               # deduped + sorted
2-element Vector{String}:
 "A"
 "B"

julia> haskey(reg, "Y")       # Y not measured -> dropped
false
```
"""
function build_regulators(edges::Vector{Tuple{String,String}}, measured::AbstractSet{String})
    reg = Dict{String,Vector{String}}()
    for (tf, sub) in edges
        if tf in measured && sub in measured
            push!(get!(reg, sub, String[]), tf)
        end
    end
    for k in keys(reg); reg[k] = sort(unique(reg[k])); end
    return reg
end

"""
    fit_substrate(sub_series, tf_matrix, w2) -> (alpha, resid, ok)

Fit one substrate's regulator coefficients by an L1 linear program (JuMP/HiGHS):
minimize `Σₜ |residualₜ| + w2 * Σⱼ |αⱼ|` subject to
`Σⱼ tf_matrix[j,t] * αⱼ + residualₜ == sub_seriesₜ`.

Arguments:

  - `sub_series::AbstractVector{Float64}` : substrate level at each of `T` timepoints.
  - `tf_matrix::AbstractMatrix{Float64}`  : `k × T` matrix of regulator TF levels
    (row `j` = TF `j`, column `t` = timepoint `t`).
  - `w2::Float64`                         : L1 sparsity penalty on the coefficients.

Returns `(alpha::Vector{Float64}, resid::Float64, ok::Bool)`. If the solver does
not reach an optimal status, returns `(zeros(k), Inf, false)`.

Throws an `ArgumentError` if `tf_matrix` has a different number of columns than
`length(sub_series)` (instead of a cryptic constraint-construction error).

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> tf = [1.0 2.0 3.0];           # 1 TF, 3 timepoints

julia> sub = [2.0, 4.0, 6.0];        # substrate = 2 * TF, exactly fittable

julia> alpha, resid, ok = fit_substrate(sub, tf, 0.0);

julia> ok && isapprox(alpha[1], 2.0; atol = 1e-6) && resid < 1e-6
true
```
"""
function fit_substrate(sub_series::AbstractVector{Float64}, tf_matrix::AbstractMatrix{Float64}, w2::Float64)
    k, T = size(tf_matrix)
    T == length(sub_series) ||
        throw(ArgumentError("fit_substrate: tf_matrix has $T columns (timepoints) " *
                            "but sub_series has length $(length(sub_series)); they must match"))
    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, a_plus[1:k]  >= 0); @variable(m, a_minus[1:k] >= 0)
    @variable(m, f_plus[1:T]  >= 0); @variable(m, f_minus[1:T] >= 0)
    @constraint(m, [t = 1:T],
        sum(tf_matrix[j, t] * (a_plus[j] - a_minus[j]) for j in 1:k) + f_plus[t] - f_minus[t] == sub_series[t])
    @objective(m, Min, sum(f_plus[t] + f_minus[t] for t in 1:T) + w2 * sum(a_plus[j] + a_minus[j] for j in 1:k))
    optimize!(m)
    if termination_status(m) != OPTIMAL
        return zeros(k), Inf, false
    end
    alpha = [value(a_plus[j]) - value(a_minus[j]) for j in 1:k]
    resid = sum(value(f_plus[t]) + value(f_minus[t]) for t in 1:T)
    return alpha, resid, true
end

"""
    fit_all(time_axis, gene_index, expression, regulators; tau=DEFAULT_TAU, w2=DEFAULT_W2)

Fit every substrate in `regulators` by calling [`fit_substrate`](@ref) on its
regulator TF trajectories, lagged by `tau` and interpolated onto the fit grid.

  - `gene_index::Dict{String,Int}` : gene name → row index into `expression`.
  - `expression`                   : `genes × timepoints` matrix (see [`load_rna_seq`](@ref)).
  - `regulators`                   : substrate → regulator-TF list (see [`build_regulators`](@ref)).

Returns `Dict{String,Vector{Tuple{String,Float64}}}` mapping each substrate to
its `(tf, alpha)` coefficients - the same shape returned by [`load_handoff`](@ref).
"""
function fit_all(time_axis, gene_index, expression, regulators; tau=DEFAULT_TAU, w2=DEFAULT_W2)
    fit_mask = time_axis .- tau .>= time_axis[1]
    fit_times = time_axis[fit_mask]; T_eff = length(fit_times)
    out = Dict{String,Vector{Tuple{String,Float64}}}()
    for sub in sort(collect(keys(regulators)))
        regs = regulators[sub]; k = length(regs); k == 0 && continue
        sub_series = expression[gene_index[sub], fit_mask]
        tf_matrix = Matrix{Float64}(undef, k, T_eff)
        for (j, tf) in enumerate(regs)
            tf_full = view(expression, gene_index[tf], :)
            for ti in 1:T_eff
                tf_matrix[j, ti] = linear_interp(fit_times[ti] - tau, time_axis, tf_full)
            end
        end
        alpha, _, _ = fit_substrate(sub_series, tf_matrix, w2)
        out[sub] = collect(zip(regs, alpha))
    end
    return out
end

# TF means (data/tf_means.csv) are the discrete means of this same NaN-interpolated
# 22-point grid (load_rna_seq), so the deviation-form multiplier below is
# mean-preserving (⟨M⟩ₜ = 1) to machine precision for every gene.
"""
    multiplier(alphas::Vector{Float64}, ratios::Vector{Float64}; q=1.0) -> Float64

Low-level mean-preserving deviation-form multiplier for one substrate:

    M = 1 + q * ( Σᵢ αᵢ (ratioᵢ − 1) ) / Σᵢ |αᵢ|

where `ratiosᵢ = TFᵢ(t) / TFᵢ_mean`. Returns `1.0` when `alphas` is empty or
sums to zero magnitude (no regulation). If each TF's ratio averages to 1 across
the trajectory, then `⟨M⟩ₜ = 1` to machine precision (mean preservation).

Throws an `ArgumentError` if `alphas` and `ratios` have different lengths
(instead of a cryptic `BoundsError` from the inner loop).

# Example
```jldoctest
julia> using TranscriptionMultiplier

julia> multiplier([2.0, -1.0], [1.0, 1.0])   # all TFs at their mean -> 1
1.0

julia> multiplier([2.0, -1.0], [1.5, 1.0])   # activator above mean
1.3333333333333333

julia> multiplier(Float64[], Float64[])      # no regulators -> 1
1.0
```

Mean preservation, `⟨M⟩ₜ = 1`, holds for any sign pattern and any `q`:

```jldoctest
julia> using TranscriptionMultiplier, Statistics

julia> alphas = [2.0, -1.0, 0.5];

julia> levels = [1.0 3.0 2.0; 4.0 1.0 1.0; 2.0 2.0 5.0];  # 3 TFs x 3 timepoints

julia> mus = vec(mean(levels; dims = 2));

julia> Ms = [multiplier(alphas, levels[:, t] ./ mus; q = 0.7) for t in 1:3];

julia> isapprox(mean(Ms), 1.0; atol = 1e-12)
true
```
"""
function multiplier(alphas::Vector{Float64}, ratios::Vector{Float64}; q=1.0)
    length(alphas) == length(ratios) ||
        throw(ArgumentError("multiplier: alphas and ratios must have equal length " *
                            "(got $(length(alphas)) and $(length(ratios)))"))
    sabs = sum(abs, alphas)
    sabs == 0 && return 1.0
    dev = sum(alphas[i] * (ratios[i] - 1.0) for i in eachindex(alphas); init = 0.0)
    return 1.0 + q * dev / sabs
end
