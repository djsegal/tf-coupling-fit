#=
joint_fit.jl: builds the joint, multi-dataset, shared-alpha L1 fit of the TF->gene
couplings, in pure Julia (JuMP/HiGHS), and provides the dataset parsers and
gene-name machinery it needs. This is the fit behind the headline cross-dataset
result: scored against the Inferelator gold standard on the common rectangle, the
joint couplings raise AUROC from 0.804 (single-dataset, Teufel-only) to 0.899.

THE MODEL
---------
For each substrate S with its Teufel-network regulators {TF_j}, and each training
dataset d (z-scored within itself), build the delayed-TF design on d's own time
grid:

    Shat_d(t) = sum_j alpha_j * TF_{j,d}(t - tau)

The shared coupling vector alpha is fit by stacking every dataset's rows into one
L1 linear program -- the same residual-L1 + w2*L1(alpha) program as the
single-dataset fit in refit.jl, but with the rows of all datasets stacked so a
single alpha must explain every dataset at once. Datasets where the substrate is
unmeasurable, or where no regulator is measurable, contribute no rows; an
unmeasurable regulator inside an otherwise-used dataset gets an all-zero design
row so the alpha index stays aligned.

INPUTS (see fetch_datasets.jl and the README reproducibility section)
  WT_unstressed_readspermillionreads.csv  the Teufel WT RNA-seq (in data)
  TF_von_Teufel.csv                        the candidate TF->gene network (in data)
  SGD_features.tab                         ORF<->standard name map (fetched)
  the four extra cell-cycle courses        Spellman, Pramila, Orlando, Kelliher (fetched)

OUTPUT
  a Dict substrate -> [(tf, alpha), ...] for every substrate with a nonzero fit,
  writable to the (substrate, tf, alpha) CSV format used throughout the package.

Reuses only the package's existing deps (JuMP, HiGHS, CSV, DataFrames, Statistics,
Printf). No new dependency, no Python.
=#
module JointFit

using JuMP, HiGHS, Statistics, Printf

export build_name_maps, key_factory, Dataset,
       parse_teufel, parse_spellman_alpha, parse_pramila, parse_orlando, parse_kelliher,
       load_teufel_regs, linear_interp, zscore_dict,
       joint_fit_substrate, joint_fit, write_alpha_csv

const DEFAULT_TAU = 20.0
const DEFAULT_W2  = 5.0

# --------------------------------------------------------------------- helpers
_is_num(s) = tryparse(Float64, strip(String(s))) !== nothing
upcasestrip(s) = uppercase(strip(String(s)))

"Minimal CSV line splitter: plain commas and double-quoted fields."
function split_csv(line::AbstractString)
    fields = String[]; buf = IOBuffer(); inq = false
    for c in line
        if c == '"'
            inq = !inq
        elseif c == ',' && !inq
            push!(fields, String(take!(buf)))
        else
            print(buf, c)
        end
    end
    push!(fields, String(take!(buf)))
    return fields
end

"Clamp-at-edges linear interpolation (matches refit.jl linear_interp)."
function linear_interp(x::Float64, xp::Vector{Float64}, fp::Vector{Float64})
    x <= xp[1]   && return fp[1]
    x >= xp[end] && return fp[end]
    i = searchsortedlast(xp, x)
    i >= length(xp) && return fp[end]
    x0, x1 = xp[i], xp[i+1]; f0, f1 = fp[i], fp[i+1]
    return f0 + (f1 - f0) * (x - x0) / (x1 - x0)
end

"Population std over finite entries; 0 if none. Used for keep-strongest-duplicate."
function _nanstd(v::AbstractVector{<:Real})
    g = Float64[x for x in v if isfinite(x)]
    isempty(g) && return 0.0
    return std(g; corrected = false)
end

# --------------------------------------------------------------------- name maps
"""orf<->std uppercased, from the WT RNA-seq table + SGD_features.tab. First
mapping wins (setdefault semantics)."""
function build_name_maps(wt_path::String, sgd_path::String)
    orf2std = Dict{String,String}(); std2orf = Dict{String,String}()
    add = function(orf, std)
        orf = upcasestrip(orf); std = upcasestrip(std)
        (isempty(orf) || isempty(std)) && return
        haskey(orf2std, orf) || (orf2std[orf] = std)
        haskey(std2orf, std) || (std2orf[std] = orf)
    end
    open(wt_path) do f
        header = split_csv(readline(f))
        idx = Dict(h => i for (i, h) in enumerate(header))
        oi = get(idx, "orf", 0); ni = get(idx, "name", 0)
        for line in eachline(f)
            isempty(line) && continue
            p = split_csv(line)
            orf = oi > 0 && oi <= length(p) ? p[oi] : ""
            nm  = ni > 0 && ni <= length(p) ? p[ni] : ""
            add(orf, nm)
        end
    end
    open(sgd_path) do f
        for line in eachline(f)
            p = split(line, '\t')
            if length(p) > 5 && p[2] == "ORF" && !isempty(p[4]) && !isempty(p[5])
                add(p[4], p[5])
            end
        end
    end
    return orf2std, std2orf
end

"""std<->orf alias key resolver against a z-scored series dict: literal name,
then std->orf, then orf->std. Returns the matching key or nothing."""
function key_factory(zser::Dict{String,Vector{Float64}},
                     std2orf::Dict{String,String}, orf2std::Dict{String,String})
    keys_set = Set(keys(zser))
    return function(name::String)
        n = uppercase(name)
        n in keys_set && return n
        a = get(std2orf, n, "")
        (a != "" && a in keys_set) && return a
        b = get(orf2std, n, "")
        (b != "" && b in keys_set) && return b
        return nothing
    end
end

# --------------------------------------------------------------------- z-scoring
"""Interpolate residual NaNs across time, then z-score each gene series
(population std). Genes with <2 finite points or zero variance are dropped."""
function zscore_dict(series::Dict{String,Vector{Float64}}, times::Vector{Float64})
    out = Dict{String,Vector{Float64}}()
    for (g, arr) in series
        a = copy(arr); mask = isnan.(a)
        if any(mask)
            good = .!mask
            sum(good) < 2 && continue
            xg = times[good]; yg = a[good]
            for k in findall(mask)
                a[k] = linear_interp(times[k], xg, yg)
            end
        end
        sd = std(a; corrected = false)
        (sd == 0 || !isfinite(sd)) && continue
        out[g] = (a .- mean(a)) ./ sd
    end
    return out
end

# --------------------------------------------------------------------- parsers
"""Teufel WT RNA-seq (data/WT_unstressed_readspermillionreads.csv) ->
(times, std-gene -> z-scored series). Key by standard `name` (uppercased),
falling back to orf2std[orf] or orf. Drop if >4 NaNs; keep strongest-varying
duplicate."""
function parse_teufel(wt_path::String, orf2std::Dict{String,String})
    rows = open(wt_path) do f
        [split_csv(line) for line in eachline(f)]
    end
    header = rows[1]
    tcols = [i for i in 1:length(header) if _is_num(header[i])]
    times0 = [parse(Float64, strip(String(header[i]))) for i in tcols]
    order = sortperm(times0; alg = MergeSort)
    times = times0[order]; tcols = tcols[order]
    hidx = Dict(strip(String(h)) => i for (i, h) in enumerate(header))
    ni = get(hidx, "name", 0); oi = get(hidx, "orf", 0)
    series = Dict{String,Vector{Float64}}()
    for r in rows[2:end]
        isempty(r) && continue
        nm  = ni > 0 && ni <= length(r) ? upcasestrip(r[ni]) : ""
        orf = oi > 0 && oi <= length(r) ? upcasestrip(r[oi]) : ""
        std = !isempty(nm) ? nm : get(orf2std, orf, orf)
        isempty(std) && continue
        vals = Float64[]
        for c in tcols
            s = c <= length(r) ? strip(String(r[c])) : ""
            push!(vals, (s == "" || s == "NaN" || s == "NA" || s == "nan") ? NaN : parse(Float64, s))
        end
        count(isnan, vals) > 4 && continue
        if haskey(series, std) && _nanstd(vals) <= _nanstd(series[std])
            continue
        end
        series[std] = vals
    end
    return times, zscore_dict(series, times)
end

"""Spellman cdc28-13 alpha-factor combined table: columns alpha0..alpha119 at
7-min spacing, rows systematic ORFs. Drop if >3 NaNs."""
function parse_spellman_alpha(path::String, orf2std::Dict{String,String})
    rows = open(path) do f; [split(line, '\t') for line in eachline(f)]; end
    header = rows[1]
    tcols = Int[]; tmins = Int[]
    for (i, h) in enumerate(header)
        m = match(r"^alpha(\d+)$", strip(String(h)))
        if m !== nothing
            push!(tcols, i); push!(tmins, parse(Int, m.captures[1]))
        end
    end
    order = sortperm(tmins); tcols = tcols[order]; tmins = tmins[order]
    times = Float64.(tmins)
    series = Dict{String,Vector{Float64}}()
    for r in rows[2:end]
        (isempty(r) || isempty(strip(String(r[1])))) && continue
        orf = upcasestrip(String(r[1])); std = get(orf2std, orf, orf)
        vals = Float64[]
        for c in tcols
            s = c <= length(r) ? strip(String(r[c])) : ""
            push!(vals, (s == "" || s == "NaN" || s == "NA") ? NaN : parse(Float64, s))
        end
        count(isnan, vals) > 3 && continue
        series[std] = vals
    end
    return times, zscore_dict(series, times)
end

"""Pramila alpha-factor PCL (SGD-archive GSE4987_setA_family.pcl). Header columns
read `... time point N min ...`; two replicate series per minute are averaged.
Rows: YORF + NAME (both systematic). Drop if >4 NaNs."""
function parse_pramila(path::String, orf2std::Dict{String,String})
    rows = open(path) do f; [split(line, '\t') for line in eachline(f)]; end
    header = rows[1]
    tcols = Int[]; tmins = Int[]
    for (i, h) in enumerate(header)
        m = match(r"time point\s+(\d+)\s+min", String(h))
        if m !== nothing
            push!(tcols, i); push!(tmins, parse(Int, m.captures[1]))
        end
    end
    uniq = sort(unique(tmins)); times = Float64.(uniq)
    body = [r for r in rows[3:end] if !isempty(r) && !isempty(strip(String(r[1])))]
    series = Dict{String,Vector{Float64}}()
    for r in body
        yorf = upcasestrip(r[1]); name2 = length(r) >= 2 ? upcasestrip(r[2]) : ""
        std = haskey(orf2std, yorf) ? orf2std[yorf] :
              (haskey(orf2std, name2) ? orf2std[name2] : yorf)
        vals = Float64[]
        for mn in uniq
            reps = Float64[]
            for (c, t) in zip(tcols, tmins)
                if t == mn && c <= length(r)
                    s = strip(String(r[c]))
                    if s != "" && s != "NaN" && s != "NA"
                        v = tryparse(Float64, s); v !== nothing && push!(reps, v)
                    end
                end
            end
            push!(vals, isempty(reps) ? NaN : mean(reps))
        end
        count(isnan, vals) > 4 && continue
        series[std] = vals
    end
    return times, zscore_dict(series, times)
end

"""Orlando 2008 set-A WildType course (SGD-archive GSE8799_setA_family.pcl).
Columns Cerevisiae_WildType_<min>min_rep{1,2}; pool reps per minute. Rows
YORF + NAME. Prefer a non-Y standard NAME; else orf2std[YORF] else NAME-or-YORF.
Drop if >4 NaNs; keep strongest-varying duplicate."""
function parse_orlando(path::String, orf2std::Dict{String,String})
    rows = open(path) do f; [split(line, '\t') for line in eachline(f)]; end
    header = rows[1]
    tcols = Int[]; tmins = Int[]
    for (i, h) in enumerate(header)
        m = match(r"Cerevisiae_WildType_(\d+)min_rep\d+", String(h))
        if m !== nothing
            push!(tcols, i); push!(tmins, parse(Int, m.captures[1]))
        end
    end
    uniq = sort(unique(tmins)); times = Float64.(uniq)
    body = [r for r in rows[2:end]
            if !isempty(r) && !isempty(strip(String(r[1]))) && strip(String(r[1])) != "EWEIGHT"]
    series = Dict{String,Vector{Float64}}()
    for r in body
        yorf = upcasestrip(r[1]); nm = length(r) >= 2 ? upcasestrip(r[2]) : ""
        std = (!isempty(nm) && !startswith(nm, "Y")) ? nm : ""
        if isempty(std)
            std = get(orf2std, yorf, isempty(nm) ? yorf : nm)
        end
        vals = Float64[]
        for mn in uniq
            reps = Float64[]
            for (c, t) in zip(tcols, tmins)
                if t == mn && c <= length(r)
                    s = strip(String(r[c]))
                    if s != "" && s != "NaN" && s != "NA"
                        v = tryparse(Float64, s); v !== nothing && push!(reps, v)
                    end
                end
            end
            push!(vals, isempty(reps) ? NaN : mean(reps))
        end
        count(isnan, vals) > 4 && continue
        if haskey(series, std) && _nanstd(vals) <= _nanstd(series[std])
            continue
        end
        series[std] = vals
    end
    return times, zscore_dict(series, times)
end

"""Kelliher GSE80474 RNA-seq (GSE80474_Scerevisiae_normalized.txt): header
`time_points` + minute columns, row labels are STANDARD gene names. Drop if >4
NaNs; keep strongest-varying duplicate."""
function parse_kelliher(path::String, orf2std::Dict{String,String})
    rows = open(path) do f; [split(line, '\t') for line in eachline(f)]; end
    header = rows[1]
    tcols = [i for i in 2:length(header) if _is_num(header[i])]
    times0 = [parse(Float64, strip(String(header[i]))) for i in tcols]
    order = sortperm(times0); times = times0[order]; tcols = tcols[order]
    series = Dict{String,Vector{Float64}}()
    for r in rows[2:end]
        (isempty(r) || isempty(strip(String(r[1])))) && continue
        std = upcasestrip(String(r[1]))
        vals = Float64[]
        for c in tcols
            s = c <= length(r) ? strip(String(r[c])) : ""
            push!(vals, (s == "" || s == "NaN" || s == "NA" || s == "nan") ? NaN : parse(Float64, s))
        end
        count(isnan, vals) > 4 && continue
        if haskey(series, std) && _nanstd(vals) <= _nanstd(series[std])
            continue
        end
        series[std] = vals
    end
    return times, zscore_dict(series, times)
end

# --------------------------------------------------------------------- regulators
"""Candidate regulator set from the Teufel network CSV: substrate (col 3) ->
sorted unique uppercased TFs (col 1). Skips rows whose tf/sub parse as floats.
No measurability filter -- every network edge is a candidate; the LP zeroes
unmeasurable regulators per-dataset."""
function load_teufel_regs(network_path::String)
    regs = Dict{String,Vector{String}}()
    open(network_path) do f
        readline(f)
        for line in eachline(f)
            isempty(line) && continue
            p = split_csv(line)
            length(p) < 3 && continue
            tf = strip(String(p[1])); sub = strip(String(p[3]))
            (isempty(tf) || isempty(sub)) && continue
            tryparse(Float64, tf)  !== nothing && continue
            tryparse(Float64, sub) !== nothing && continue
            lst = get!(regs, uppercase(sub), String[])
            uppercase(tf) in lst || push!(lst, uppercase(tf))
        end
    end
    for k in keys(regs); regs[k] = sort(regs[k]); end
    return regs
end

# --------------------------------------------------------------------- dataset
struct Dataset
    name::String
    times::Vector{Float64}
    zser::Dict{String,Vector{Float64}}
    key_fn::Function
end

# --------------------------------------------------------------------- design
"""(k x T) delayed-TF design for `regs` on `times`/`zser`, plus a `measurable`
mask. Unmeasurable regulators keep an all-zero row."""
function lagged_design(regs::Vector{String}, key_fn, times::Vector{Float64},
                       zser::Dict{String,Vector{Float64}}, tau::Float64)
    T = length(times); k = length(regs)
    M = zeros(k, T); measurable = falses(k)
    for (j, tf) in enumerate(regs)
        kk = key_fn(tf); kk === nothing && continue
        full = zser[kk]
        for ti in 1:T
            M[j, ti] = linear_interp(times[ti] - tau, times, full)
        end
        measurable[j] = true
    end
    return M, measurable
end

# --------------------------------------------------------------------- joint LP
"""Shared-alpha L1 fit across datasets for one substrate. Stacks the delayed-TF
rows + z-scored target of every dataset where the substrate and >=1 regulator are
measurable. Returns (alpha::Vector, ok::Bool, n_obs, n_datasets_used)."""
function joint_fit_substrate(regs::Vector{String}, datasets::Vector{Dataset},
                             sub::String; tau::Float64 = DEFAULT_TAU, w2::Float64 = DEFAULT_W2)
    k = length(regs)
    blocks = Tuple{Matrix{Float64},Vector{Float64}}[]
    for d in datasets
        sk = d.key_fn(sub); sk === nothing && continue
        M, meas = lagged_design(regs, d.key_fn, d.times, d.zser, tau)
        any(meas) || continue
        push!(blocks, (M, d.zser[sk]))
    end
    isempty(blocks) && return (Float64[], false, 0, 0)
    Ttot = sum(size(b[1], 2) for b in blocks)

    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, a_plus[1:k]  >= 0)
    @variable(m, a_minus[1:k] >= 0)
    @variable(m, r_plus[1:Ttot]  >= 0)
    @variable(m, r_minus[1:Ttot] >= 0)
    row = 0
    for (M, target) in blocks
        Td = size(M, 2)
        for t in 1:Td
            row += 1
            @constraint(m,
                sum(M[j, t] * (a_plus[j] - a_minus[j]) for j in 1:k) +
                r_plus[row] - r_minus[row] == target[t])
        end
    end
    @objective(m, Min,
        sum(r_plus[t] + r_minus[t] for t in 1:Ttot) +
        w2 * sum(a_plus[j] + a_minus[j] for j in 1:k))
    optimize!(m)
    if termination_status(m) != OPTIMAL
        return (Float64[], false, Ttot, length(blocks))
    end
    alpha = [value(a_plus[j]) - value(a_minus[j]) for j in 1:k]
    return (alpha, true, Ttot, length(blocks))
end

"""Joint fit over all Teufel substrates. Returns Dict substrate ->
[(tf, alpha), ...] for substrates with a successful, nonzero fit."""
function joint_fit(datasets::Vector{Dataset}, regs_of::Dict{String,Vector{String}};
                   tau::Float64 = DEFAULT_TAU, w2::Float64 = DEFAULT_W2)
    out = Dict{String,Vector{Tuple{String,Float64}}}()
    for sub in sort(collect(keys(regs_of)))
        regs = regs_of[sub]; isempty(regs) && continue
        alpha, ok, _, _ = joint_fit_substrate(regs, datasets, sub; tau = tau, w2 = w2)
        if ok && !isempty(alpha) && any(alpha .!= 0)
            out[sub] = collect(zip(regs, alpha))
        end
    end
    return out
end

"""Write a coupling dict to a (substrate, tf, alpha) CSV, 6dp, sorted by substrate
then regulator order. CRLF line endings + %.6f match the released alpha CSVs
byte-for-byte."""
function write_alpha_csv(path::String, joint_alpha::Dict{String,Vector{Tuple{String,Float64}}})
    open(path, "w") do f
        print(f, "substrate,tf,alpha\r\n")
        for sub in sort(collect(keys(joint_alpha)))
            for (tf, a) in joint_alpha[sub]
                av = a == 0.0 ? 0.0 : a   # normalize -0.0 -> 0.0
                print(f, sub, ",", tf, ",", @sprintf("%.6f", av), "\r\n")
            end
        end
    end
end

end # module
