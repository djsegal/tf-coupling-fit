#!/usr/bin/env julia
#=
joint_score.jl: reproduce the headline cross-dataset result.

Scores the single-dataset (Teufel-only) couplings and the joint multi-dataset
couplings against the signed Inferelator gold standard on the COMMON RECTANGLE
(the TF x target block all three networks -- Teufel, joint, leave-Kelliher-out --
and the gold cover), and runs the paired DeLong test of joint vs Teufel.

Prints, and asserts within tolerance:
  Teufel AUROC 0.804   joint AUROC 0.899   (trapezoid, on the common rectangle)
  DeLong joint - Teufel = +0.071  95% CI [0.052, 0.091]  p = 1.4e-12

Two modes (default = scored, which needs no dataset fetch):

  julia --project=.. transcription_multiplier/joint_score.jl
      Scores the SHIPPED coupling CSVs (data/tf_network_fitted.csv +
      data/joint_multidataset_alpha.csv + data/loo_kelliher_alpha.csv). Fast.
      This is the path an independent reader runs to regenerate the headline.

  julia --project=.. transcription_multiplier/joint_score.jl --refit
      Additionally re-solves the joint L1 LP from the fetched raw datasets
      (joint_fit.jl), writes the couplings to /tmp, and scores those instead, so
      the full raw-data-to-headline chain is exercised. Requires fetch_datasets.jl
      to have populated data/external/ first.

Gold standard provenance: the Inferelator yeast gold standard
(flatironinstitute/inferelator, data/yeast/gold_standard.tsv). Fetched by
fetch_datasets.jl to data/external/inferelator_gold_standard.tsv.

Pure Julia. The DeLong variance uses the fast Sun & Xu (2014) algorithm; the
normal tail and trapezoid AUROC are computed in-file (no new dependency).
=#

include(joinpath(@__DIR__, "joint_fit.jl"))
# Import only what the scorer needs, by name, so this file can be `include`d into
# a scope that already defines its own `linear_interp` (e.g. runtests.jl) without
# an export collision.
using .JointFit: build_name_maps, key_factory, Dataset,
                 parse_teufel, parse_spellman_alpha, parse_pramila, parse_orlando,
                 parse_kelliher, load_teufel_regs, joint_fit, write_alpha_csv
using Statistics, Printf

const HERE = @__DIR__
const REPO = normpath(joinpath(HERE, ".."))
const DATA = joinpath(HERE, "data")
const EXT  = joinpath(DATA, "external")

const WT       = joinpath(REPO, "data", "WT_unstressed_readspermillionreads.csv")
const NETWORK  = joinpath(REPO, "data", "TF_von_Teufel.csv")
const SGD      = joinpath(EXT, "SGD_features.tab")
const GOLD     = joinpath(EXT, "inferelator_gold_standard.tsv")

const TEUFEL_CSV = joinpath(DATA, "tf_network_fitted.csv")
const JOINT_CSV  = joinpath(DATA, "joint_multidataset_alpha.csv")
const LOO_CSV    = joinpath(DATA, "loo_kelliher_alpha.csv")

const TAU = 20.0
const W2  = 5.0

# headline anchors and tolerances
const ANCHOR_TEUFEL = 0.804
const ANCHOR_JOINT  = 0.899
const ANCHOR_DIFF   = 0.071
const ANCHOR_CI     = (0.052, 0.091)

# --------------------------------------------------------------------- gold
"""Load the signed Inferelator gold standard. Header row = TF ORFs (columns);
first cell of each data row = target ORF; nonzero signed entries are edges.
Returns (gold::Dict{(tf_std,tgt_std)=>sign}, all_tf_std::Set, all_tgt_std::Set)."""
function load_gold(path::String, orf2std::Dict{String,String})
    rows = [split(line, '\t') for line in eachline(path)]
    clean(c) = uppercase(strip(strip(String(c)), ['"', ' ']))
    tf_orfs = [clean(c) for c in rows[1]]
    gold = Dict{Tuple{String,String},Int}()
    tf_cov = Set{String}(); tgt_cov = Set{String}()
    for r in rows[2:end]
        (isempty(r) || isempty(strip(String(r[1])))) && continue
        tgt_std = get(orf2std, clean(r[1]), clean(r[1]))
        push!(tgt_cov, tgt_std)
        for c in 2:length(r)
            tfi = c - 1
            tfi > length(tf_orfs) && break
            cell = strip(strip(String(r[c]), ['"', ' ']))
            (cell == "" || cell == "0") && continue
            v = tryparse(Float64, cell); v === nothing && continue
            s = Int(round(v)); s == 0 && continue
            tf_std = get(orf2std, tf_orfs[tfi], tf_orfs[tfi])
            push!(tf_cov, tf_std)
            gold[(tf_std, tgt_std)] = s
        end
    end
    all_tf  = Set(get(orf2std, o, o) for o in tf_orfs if !isempty(o))
    all_tgt = Set{String}()
    for r in rows[2:end]
        (isempty(r) || isempty(strip(String(r[1])))) && continue
        push!(all_tgt, get(orf2std, clean(r[1]), clean(r[1])))
    end
    return gold, all_tf, all_tgt
end

"""Load (substrate, tf, alpha) couplings uppercased; drop alpha==0 unless asked."""
function load_edges(path::String; include_zeros::Bool = false)
    out = Tuple{String,String,Float64}[]
    open(path) do f
        header = JointFit.split_csv(readline(f))
        idx = Dict(h => i for (i, h) in enumerate(header))
        si = idx["substrate"]; ti = idx["tf"]; ai = idx["alpha"]
        for line in eachline(f)
            isempty(line) && continue
            p = JointFit.split_csv(line)
            a = parse(Float64, strip(p[ai]))
            (!include_zeros && a == 0.0) && continue
            push!(out, (JointFit.upcasestrip(p[si]), JointFit.upcasestrip(p[ti]), a))
        end
    end
    return out
end

# --------------------------------------------------------------------- AUROC
"""Trapezoidal integral of y wrt x (np.trapz)."""
function trapz(y::Vector{Float64}, x::Vector{Float64})
    s = 0.0
    @inbounds for i in 2:length(x)
        s += (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
    end
    return s
end

"""Trapezoid AUROC: stable descending sort by score (ties broken by index),
prepend the origin, integrate TPR over FPR."""
function auroc_trap(scores::Vector{Float64}, labels::Vector{Int})
    n = length(scores)
    order = sortperm(1:n; by = i -> (-scores[i], i))
    y = labels[order]; P = sum(y); N = n - P
    (P == 0 || N == 0) && return NaN
    tp = cumsum(y); fp = cumsum(1 .- y)
    tpr = vcat(0.0, tp ./ P); fpr = vcat(0.0, fp ./ N)
    return trapz(Float64.(tpr), Float64.(fpr))
end

# --------------------------------------------------------------------- DeLong
"Midranks with average ranks for ties (Sun & Xu)."
function midrank(x::Vector{Float64})
    order = sortperm(x; alg = MergeSort)
    xs = x[order]; r = Vector{Float64}(undef, length(x)); i = 1
    while i <= length(xs)
        j = i
        while j + 1 <= length(xs) && xs[j+1] == xs[i]; j += 1; end
        for k in i:j; r[order[k]] = 0.5 * (i - 1 + j - 1) + 1.0; end
        i = j + 1
    end
    return r
end

"erfc via the Numerical Recipes rational approximation (no SpecialFunctions dep)."
function _erfc(x::Float64)
    z = abs(x); t = 1.0 / (1.0 + 0.5z)
    ans = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
        t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
        t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    return x >= 0 ? ans : 2.0 - ans
end
_normal_sf(x) = 0.5 * _erfc(x / sqrt(2))

"""Fast DeLong (Sun & Xu 2014) for two correlated AUCs on the same labels.
Returns (aucA, aucB, diff=A-B, se, z, p, ci_lo, ci_hi). AUCs here are the
Mann-Whitney form used by the variance estimator."""
function delong(sA::Vector{Float64}, sB::Vector{Float64}, labels::Vector{Int})
    pos = labels .== 1; neg = labels .== 0
    m = count(pos); n = count(neg)
    preds = [sA, sB]; aucs = zeros(2); v01 = zeros(2, m); v10 = zeros(2, n)
    for r in 1:2
        X = preds[r][pos]; Y = preds[r][neg]
        tx = midrank(X); ty = midrank(Y); tz = midrank(vcat(X, Y))
        aucs[r] = (sum(tz[1:m]) - m * (m + 1) / 2) / (m * n)
        v01[r, :] = (tz[1:m] .- tx) ./ n
        v10[r, :] = 1.0 .- (tz[m+1:end] .- ty) ./ m
    end
    S = cov(v01') ./ m .+ cov(v10') ./ n
    diff = aucs[1] - aucs[2]
    var = S[1, 1] - S[1, 2] - S[2, 1] + S[2, 2]
    se = sqrt(var); z = diff / se; p = 2 * _normal_sf(abs(z))
    return aucs[1], aucs[2], diff, se, z, p, diff - 1.96se, diff + 1.96se
end

# --------------------------------------------------------------------- refit
"Re-solve the joint and leave-Kelliher-out couplings from the fetched raw data."
function refit_from_raw(orf2std, std2orf)
    regs_of = load_teufel_regs(NETWORK)
    mk(name, parser, path) = begin
        times, zser = parser(path, orf2std)
        Dataset(name, times, zser, key_factory(zser, std2orf, orf2std))
    end
    teufel   = mk("teufel",   parse_teufel,         WT)
    spellman = mk("spellman", parse_spellman_alpha, joinpath(EXT, "spellman_combined.txt"))
    pramila  = mk("pramila",  parse_pramila,        joinpath(EXT, "pramila.pcl"))
    orlando  = mk("orlando",  parse_orlando,        joinpath(EXT, "orlando_setA.pcl"))
    kelliher = mk("kelliher", parse_kelliher,       joinpath(EXT, "GSE80474_Scerevisiae_normalized.txt"))
    @info "re-solving joint LP from raw datasets (this takes a few seconds)"
    joint5 = joint_fit(Dataset[teufel, spellman, pramila, orlando, kelliher], regs_of; tau = TAU, w2 = W2)
    jointK = joint_fit(Dataset[teufel, spellman, pramila, orlando],            regs_of; tau = TAU, w2 = W2)
    jpath = "/tmp/joint_multidataset_alpha.csv"; lpath = "/tmp/loo_kelliher_alpha.csv"
    write_alpha_csv(jpath, joint5); write_alpha_csv(lpath, jointK)
    return jpath, lpath
end

# --------------------------------------------------------------------- main
function main(; refit::Bool = false)
    for f in (WT, NETWORK)
        isfile(f) || error("missing required file: $f")
    end
    for f in (SGD, GOLD)
        isfile(f) || error("missing fetched file: $f -- run transcription_multiplier/fetch_datasets.jl first")
    end

    orf2std, std2orf = build_name_maps(WT, SGD)
    gold, gtfs, gtgts = load_gold(GOLD, orf2std)

    joint_path = JOINT_CSV; loo_path = LOO_CSV
    if refit
        joint_path, loo_path = refit_from_raw(orf2std, std2orf)
    else
        for f in (TEUFEL_CSV, JOINT_CSV, LOO_CSV)
            isfile(f) || error("missing shipped coupling CSV: $f")
        end
    end

    teufel = load_edges(TEUFEL_CSV)
    joint  = load_edges(joint_path)
    loo    = load_edges(loo_path)

    # common rectangle: TFs and targets present in gold AND all three networks
    tfc = sort(collect(intersect(Set(t for (_, t, _) in teufel), Set(t for (_, t, _) in joint),
                                 Set(t for (_, t, _) in loo), gtfs)))
    tgc = sort(collect(intersect(Set(s for (s, _, _) in teufel), Set(s for (s, _, _) in joint),
                                 Set(s for (s, _, _) in loo), gtgts)))
    tfs = Set(tfc); tgs = Set(tgc)
    gpos = Set((tf, sub) for ((tf, sub), _) in gold if (tf in tfs) && (sub in tgs))
    cells = [(tf, sub) for tf in tfc for sub in tgc]
    labels = [c in gpos ? 1 : 0 for c in cells]

    amt = Dict((tf, sub) => a for (sub, tf, a) in teufel)
    amj = Dict((tf, sub) => a for (sub, tf, a) in joint)
    sT = [abs(get(amt, c, 0.0)) for c in cells]
    sJ = [abs(get(amj, c, 0.0)) for c in cells]

    aucT = auroc_trap(sT, labels)
    aucJ = auroc_trap(sJ, labels)
    _, _, diff, se, z, p, lo, hi = delong(sJ, sT, labels)

    @printf("\ncommon rectangle: %d TF x %d target = %d cells, %d gold positives\n",
            length(tfc), length(tgc), length(cells), sum(labels))
    @printf("\n  AUROC (trapezoid, common rectangle)\n")
    @printf("    Teufel (single-dataset) : %.3f   (headline 0.804)\n", aucT)
    @printf("    joint  (multi-dataset)  : %.3f   (headline 0.899)\n", aucJ)
    @printf("\n  paired DeLong test, joint vs Teufel\n")
    @printf("    diff = +%.3f   95%% CI [%.3f, %.3f]   z = %.2f   p = %.2e\n", diff, lo, hi, z, p)
    @printf("    (headline +0.071 [0.052, 0.091] p = 1.4e-12)\n")

    # assertions (loose tolerance so a future re-fit drift still reproduces)
    ok = true
    ok &= isapprox(aucT, ANCHOR_TEUFEL; atol = 0.01)
    ok &= isapprox(aucJ, ANCHOR_JOINT;  atol = 0.01)
    ok &= isapprox(diff, ANCHOR_DIFF;   atol = 0.01)
    ok &= (lo > 0) && isapprox(lo, ANCHOR_CI[1]; atol = 0.01) && isapprox(hi, ANCHOR_CI[2]; atol = 0.01)
    println("\n", ok ? "HEADLINE REPRODUCED (within tolerance)." :
                       "HEADLINE OUT OF TOLERANCE -- see numbers above.")
    return aucT, aucJ, diff, lo, hi, ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    refit = "--refit" in ARGS
    _, _, _, _, _, ok = main(; refit = refit)
    exit(ok ? 0 : 1)
end
