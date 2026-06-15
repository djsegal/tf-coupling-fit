#=
refit.jl: Julia reproduction of the tf-coupling-fit linear program (JuMP/HiGHS),
plus leave-one-timepoint-out cross-validation. Same solver stack as the generating
notebook, so it reproduces the committed alpha exactly. Used by runtests.jl.

    include("refit.jl")
    ta, genes, expr = load_rna_seq(joinpath(@__DIR__,"..","data","WT_unstressed_readspermillionreads.csv"))
    edges = load_tf_network(joinpath(@__DIR__,"..","data","TF_von_Teufel.csv"))
    reg   = build_regulators(edges, Set(genes))
    res   = fit_all(ta, Dict(g=>i for (i,g) in enumerate(genes)), expr, reg)
=#
using CSV, DataFrames, JuMP, HiGHS, Statistics

const MOLECULES_PER_CELL = 60_000
const RPM_TO_MOLECULES   = MOLECULES_PER_CELL / 1_000_000   # 0.06
const DEFAULT_TAU      = 20.0
const DEFAULT_W2       = 5.0
const DEFAULT_MAX_NANS = 4

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

function load_rna_seq(path::String; max_nans::Int = DEFAULT_MAX_NANS)
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

function load_tf_network(csv_path::String)
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

function fit_substrate(sub_series::AbstractVector{Float64}, tf_matrix::AbstractMatrix{Float64}, w2::Float64)
    k, T = size(tf_matrix)
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
"Mean-preserving deviation-form multiplier for one substrate."
function multiplier(alphas::Vector{Float64}, ratios::Vector{Float64}; q=1.0)
    sabs = sum(abs, alphas)
    sabs == 0 && return 1.0
    dev = sum(alphas[i] * (ratios[i] - 1.0) for i in eachindex(alphas); init = 0.0)
    return 1.0 + q * dev / sabs
end
