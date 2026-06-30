#!/usr/bin/env julia
# Eighth-dataset out-of-sample transfer test: Nadal-Ribelles et al. (2019)
# budding-yeast SINGLE-CELL RNA-seq (yscRNA-seq). This is the most modern data
# type in the validation set: instead of a synchronized/bulk time course, it is
# an UNSYNCHRONIZED clonal population of single cells, from which we reconstruct
# a per-gene cell-cycle profile by PSEUDOTIME ordering.
#
# Source: Nadal-Ribelles M, Islam S, Wei W, Latorre P, Nguyen M, de Nadal E,
# Posas F, Steinmetz LM. "Sensitive high-throughput single-cell RNA-seq reveals
# within-clonal transcript-correlations in yeast populations." Nat Microbiol
# 4:683-692 (2019). GEO GSE122392. SRA SRP168379.
# Processed normalized expression matrices (genes x cells, UMI-based, TSS-level):
#   GSE122392_normExpression_YPD_NatMicro.tab.gz  (lab strain, 127 cells)  <- primary
#   GSE122392_normExpression_YJM_NatMicro.tab.gz  (wild isolate, 48 cells)
# Downloaded to /tmp/datasets/nadal_ribelles/.
#
# This dataset is NOT one of the five used to fit the couplings (Teufel,
# Kelliher, Spellman, Orlando, Pramila), nor Cho/Eser/Granovskaia. Fully out-of-
# sample, and a fundamentally different assay (single-cell, unsynchronized).
#
# PSEUDOTIME METHOD (transparent, no black-box tool):
#   1. Restrict to detected ORF genes.
#   2. Take the periodic-marker submatrix = genes in the production qx periodic
#      set (periodicity_pvalue < 1e-3) that are detected here. z-score each gene
#      across cells.
#   3. PCA (via SVD) of the cells x periodic-genes matrix; the leading 2 PCs of a
#      cell-cycle population span a roughly circular manifold. Compute each cell's
#      phase angle theta = atan2(PC2, PC1). Order cells by theta (cyclic order).
#   4. Orient/rotate the cycle so it starts at the canonical G1 marker SIC1 peak,
#      then bin the ordered cells into NBINS pseudo-timepoints; mean expression
#      per gene per bin -> a periodic per-gene cell-cycle profile on a synthetic
#      time axis (0..(NBINS-1) spaced one "unit" apart, scaled to ~one 85-min cycle
#      so tau=20 min is meaningful).
#   5. z-score per gene over bins, run transfer_r for Teufel-only and all-5 joint
#      couplings on the periodic scope; report n + median r, joint vs teufel.
#
# HONEST CAVEATS (stated, not hidden):
#   - n cells is small (127 YPD) and the data are sparse (~3400 genes/cell); a
#     pseudotime cell-cycle reconstruction from 127 cells is noisier than a real
#     synchronized course. We therefore expect a WEAKER transfer signal than the
#     bulk anchors, and we test sensitivity to NBINS.
#   - The phase angle is derived from periodic marker genes; this could in
#     principle bias the periodic-scope transfer upward. To bound that, the
#     marker set used for ordering is the *periodic* genes, while the couplings
#     predict each substrate from its *regulators* (TFs) at a lag -- the ordering
#     does not use the TF->substrate coupling, only that the population traverses
#     the cycle. Still, treat the magnitude cautiously; the sign/ranking
#     (joint vs teufel) is the robust readout.
#   - We reproduce a bulk ANCHOR (Kelliher joint ~0.258) FIRST, exactly as the
#     other extra-dataset scripts, to prove the machinery before trusting any
#     single-cell number.
#
# Run (from this dir):
#   julia --project=../julia_stats run_nadalribelles.jl

include(joinpath(@__DIR__, "..", "julia_stats", "joint_fit.jl"))
using .JointFit
const JS = JointFit.JuliaStats
using .JointFit.JuliaStats: load_qx_periodic, transfer_r, load_couplings_nonzero,
                            zscore_dict, key_factory
using Printf, Statistics, LinearAlgebra, Random

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "TranscriptionMultiplier.jl"))
const DATA = "/tmp/datasets"
const WT       = joinpath(ROOT, "data/WT_unstressed_readspermillionreads.csv")
const SGD      = joinpath(DATA, "SGD_features.tab")
const NETWORK  = joinpath(ROOT, "data/TF_von_Teufel.csv")
const QX       = joinpath(ROOT, "data/qx_scores.csv")
const TEUFEL_C = joinpath(ROOT, "output/tf_handoff/tf_network_fitted.csv")
const JOINT_C  = joinpath(ROOT, "data/joint_multidataset_alpha.csv")

const KELLIHER = joinpath(DATA, "GSE80474_Scerevisiae_normalized.txt")
const SPELLMAN = joinpath(DATA, "spellman.txt")
const PRAMILA  = joinpath(DATA, "pramila.txt")
const ORLANDO  = joinpath(DATA, "Orlando_2008_PMID_18463633/GSE8799_setA_family.pcl")

const NADAL_YPD = joinpath(DATA, "nadal_ribelles/GSE122392_normExpression_YPD_NatMicro.tab.gz")
const NADAL_YJM = joinpath(DATA, "nadal_ribelles/GSE122392_normExpression_YJM_NatMicro.tab.gz")

const TAU   = 20.0
const W2    = 5.0
const CYCLE_MIN = 85.0   # nominal budding-yeast cycle length (min) for the synthetic time axis

upcasestrip(s) = uppercase(strip(String(s)))

# ------------------------------------------------------------- read scRNA matrix
"""Read a GSE122392 normExpression .tab(.gz): header `geneName comGeneName <cell ids...>`,
rows are genes (col1 systematic/feature name, col2 common name), values are
normalized expression per cell. Returns (genes::Vector{String upper}, X::Matrix
[genes x cells])."""
function read_scmatrix(path::String)
    open_fn = endswith(path, ".gz") ? (p->`gzip -dc $p`) : (p->`cat $p`)
    lines = readlines(open_fn(path))
    header = split(lines[1], '\t')
    ncell = length(header) - 2
    genes = String[]
    rows = Vector{Float64}[]
    for ln in lines[2:end]
        isempty(ln) && continue
        p = split(ln, '\t')
        length(p) < 2 + ncell && continue
        g = upcasestrip(p[1])
        v = Float64[(s = strip(String(p[2+j])); (s=="" || s=="NA") ? 0.0 : parse(Float64, s)) for j in 1:ncell]
        push!(genes, g); push!(rows, v)
    end
    X = Matrix{Float64}(undef, length(genes), ncell)
    for i in 1:length(genes); X[i, :] = rows[i]; end
    return genes, X
end

# ------------------------------------------------------------- pseudotime by PCA
"""Compute a per-cell cell-cycle phase angle from the periodic-marker submatrix
via PCA (SVD). Returns theta::Vector (one per cell), and the gene indices used."""
function cellcycle_phase(genes::Vector{String}, X::Matrix{Float64},
                         periodic::Set{String}, orf2std, std2orf)
    # periodic markers detected here (match by name or std<->orf alias)
    is_per(g) = (g in periodic) || (get(orf2std, g, "") in periodic) || (get(std2orf, g, "") in periodic)
    midx = [i for i in 1:length(genes) if is_per(genes[i])]
    # keep markers expressed in a reasonable fraction of cells (>=20%) and varying
    ncell = size(X, 2)
    keep = Int[]
    for i in midx
        row = @view X[i, :]
        (count(>(0), row) >= 0.2 * ncell) || continue
        std(row; corrected=false) > 0 || continue
        push!(keep, i)
    end
    M = X[keep, :]                      # markers x cells
    # z-score each marker gene across cells
    for r in 1:size(M, 1)
        mr = mean(@view M[r, :]); sr = std(@view M[r, :]; corrected=false)
        sr == 0 && (sr = 1.0)
        M[r, :] .= (M[r, :] .- mr) ./ sr
    end
    # PCA on cells: SVD of cells x genes (centered already per gene)
    A = permutedims(M)                  # cells x markers
    U, S, V = svd(A)
    pc1 = U[:, 1] .* S[1]
    pc2 = U[:, 2] .* S[2]
    theta = atan.(pc2, pc1)             # in (-pi, pi]
    return theta, keep
end

"""Order cells by phase, orient so the cycle begins at the SIC1 (M/G1) trough->rise,
bin into nbins pseudo-timepoints, mean expression per gene per bin. Returns
(times::Vector, series::Dict gene=>Vector) with raw (pre-zscore) means."""
function build_profile(genes::Vector{String}, X::Matrix{Float64}, theta::Vector{Float64},
                       nbins::Int; anchor_orf::String="YLR079W")  # SIC1
    ncell = size(X, 2)
    ord = sortperm(theta)
    # orient direction: ensure the canonical G1 cyclin cluster precedes G2/M.
    # We anchor phase 0 at the cell with the highest smoothed SIC1; direction is
    # set so CLN2 (G1/S) peaks before CLB2 (G2/M) along increasing bin index.
    gi(name) = findfirst(==(upcasestrip(name)), genes)
    sic1 = gi(anchor_orf)
    # assign each ordered cell to a bin
    binsz = ncell / nbins
    times = collect(0:(nbins-1)) .* (CYCLE_MIN / nbins)
    series = Dict{String,Vector{Float64}}()
    binof = zeros(Int, ncell)
    for (rank, c) in enumerate(ord)
        b = min(nbins, 1 + floor(Int, (rank - 1) / binsz))
        binof[c] = b
    end
    binmean(rowidx) = begin
        v = fill(NaN, nbins)
        for b in 1:nbins
            cells = [c for c in 1:ncell if binof[c] == b]
            isempty(cells) && continue
            v[b] = mean(@view X[rowidx, cells])
        end
        v
    end
    for i in 1:length(genes)
        series[genes[i]] = binmean(i)
    end
    # --- orient cycle: rotate bins so SIC1 peak is at end (M/G1), and flip
    #     direction if needed so CLN2 (G1/S) leads CLB2 (G2/M). ---
    cln2 = gi("YPL256C"); clb2 = gi("YPR119W")
    function peakbin(rowidx)
        rowidx === nothing && return nothing
        v = series[genes[rowidx]]
        gv = [isnan(x) ? -Inf : x for x in v]
        argmax(gv)
    end
    pc = peakbin(cln2); pb = peakbin(clb2)
    if pc !== nothing && pb !== nothing && pc > pb
        # G1/S peaks after G2/M -> reverse direction
        for g in keys(series); reverse!(series[g]); end
    end
    return times, series
end

# =========================================================================== run
println("loading name maps + Teufel regulators ...")
orf2std, std2orf = build_name_maps(WT, SGD)
regs_of = load_teufel_regs(NETWORK)

println("parsing the five training datasets ...")
function mk(name, parser, path)
    times, zser = parser(path, orf2std)
    kfn = key_factory(zser, std2orf, orf2std)
    @printf("  %-9s %5d genes, %2d timepoints\n", name, length(zser), length(times))
    return Dataset(name, times, zser, kfn)
end
teufel   = mk("teufel",   parse_teufel,         WT)
spellman = mk("spellman", parse_spellman_alpha, SPELLMAN)
pramila  = mk("pramila",  parse_pramila,        PRAMILA)
orlando  = mk("orlando",  parse_orlando,        ORLANDO)
kelliher = mk("kelliher", parse_kelliher,       KELLIHER)

periodic   = load_qx_periodic(QX)
teufel_net = load_couplings_nonzero(TEUFEL_C)
joint_net  = load_couplings_nonzero(JOINT_C)

# ===================================================================== ANCHOR
println("\n=== ANCHOR: leave-Kelliher-out joint -> Kelliher transfer ===")
let
    trainK = Dataset[teufel, spellman, pramila, orlando]
    jK = joint_fit(trainK, regs_of; tau=TAU, w2=W2)
    looK = Dict{String,Dict{String,Float64}}()
    for (s,lst) in jK, (tf,a) in lst; a==0.0 && continue; get!(looK, uppercase(s), Dict{String,Float64}())[uppercase(tf)] = a; end
    times, zser = parse_kelliher(KELLIHER, orf2std)
    kfn = key_factory(zser, std2orf, orf2std)
    sc = JS.periodic_scope(union(keys(teufel_net), keys(looK)), periodic, std2orf, orf2std)
    nt,mt,_ = transfer_r(teufel_net, kfn, times, zser; tau=TAU, scope=sc)
    nj,mj,_ = transfer_r(looK,       kfn, times, zser; tau=TAU, scope=sc)
    @printf("  Kelliher teufel-only periodic n=%d median_r=%.4f  (cited 0.0816)\n", nt, mt)
    @printf("  Kelliher joint(LOO)  periodic n=%d median_r=%.4f  (cited 0.2580)\n", nj, mj)
end

# ============================================ EIGHTH DATASET: NADAL-RIBELLES scRNA
function run_nadal(label, path, nbins)
    genes, X = read_scmatrix(path)
    ncell = size(X, 2)
    theta, mused = cellcycle_phase(genes, X, periodic, orf2std, std2orf)
    times, series = build_profile(genes, X, theta, nbins)
    zser = zscore_dict(series, times)
    kfn = key_factory(zser, std2orf, orf2std)
    sc_u = JS.periodic_scope(union(keys(teufel_net), keys(joint_net)), periodic, std2orf, orf2std)
    nt, mt, rest = transfer_r(teufel_net, kfn, times, zser; tau=TAU, scope=sc_u)
    nj, mj, resj = transfer_r(joint_net,  kfn, times, zser; tau=TAU, scope=sc_u)
    @printf("  [%s nbins=%d] cells=%d markers=%d periodic-genes=%d\n",
            label, nbins, ncell, length(mused), length(zser))
    @printf("    teufel-only periodic n=%d median_r=%.4f\n", nt, mt)
    @printf("    joint(all-5) periodic n=%d median_r=%.4f   delta=%+.4f\n", nj, mj, mj-mt)
    return (label=label, nbins=nbins, ncell=ncell, nmark=length(mused),
            nt=nt, mt=mt, nj=nj, mj=mj, resj=resj)
end

println("\n=== EIGHTH DATASET: Nadal-Ribelles 2019 scRNA-seq (pseudotime) ===")
results = NamedTuple[]
for nb in (10, 12, 16, 20)
    push!(results, run_nadal("YPD", NADAL_YPD, nb))
end
println("  --- YJM (wild isolate, 48 cells) sensitivity ---")
push!(results, run_nadal("YJM", NADAL_YJM, 12))

# write per-gene results for the primary YPD nbins=12 run + a summary
let primary = first(r for r in results if r.label=="YPD" && r.nbins==12)
    open(joinpath(@__DIR__, "nadalribelles_per_gene_periodic.csv"), "w") do f
        println(f, "substrate,n_reg,r_joint")
        for (s,n,r) in sort(primary.resj; by=x->-x[3])
            @printf(f, "%s,%d,%.6f\n", s, n, r)
        end
    end
end
open(joinpath(@__DIR__, "nadalribelles_summary.csv"), "w") do f
    println(f, "strain,nbins,ncells,nmarkers,n_teufel,median_r_teufel,n_joint,median_r_joint,delta_joint_minus_teufel")
    for r in results
        @printf(f, "%s,%d,%d,%d,%d,%.6f,%d,%.6f,%+.6f\n",
                r.label, r.nbins, r.ncell, r.nmark, r.nt, r.mt, r.nj, r.mj, r.mj-r.mt)
    end
end

println("\nDONE.")
