#!/usr/bin/env julia
# Tests for the cell-cycle transcription multiplier and the fit reproduction.
#   julia --project=. transcription_multiplier/runtests.jl
# Run from the repo root (uses ../data and ../output relative to this file).

using Test, CSV, DataFrames, Statistics
include(joinpath(@__DIR__, "refit.jl"))

const REPO = normpath(joinpath(@__DIR__, ".."))
const EXPR = joinpath(REPO, "data", "WT_unstressed_readspermillionreads.csv")
const NET  = joinpath(REPO, "data", "TF_von_Teufel.csv")
const COMMITTED = joinpath(REPO, "output", "tf_handoff", "alpha_edges.csv")
const FITTED = joinpath(@__DIR__, "data", "tf_network_fitted.csv")

@testset "tf-coupling-fit multiplier" begin

    # load once
    time_axis, genes, expr = load_rna_seq(EXPR)
    gidx = Dict(g => i for (i, g) in enumerate(genes))
    reg  = build_regulators(load_tf_network(NET), Set(genes))

    @testset "fit reproduces committed alpha exactly (JuMP/HiGHS)" begin
        res = fit_all(time_axis, gidx, expr, reg)
        mine = Dict{Tuple{String,String},Float64}()
        for (sub, lst) in res, (tf, a) in lst
            mine[(sub, tf)] = a
        end
        ref = Dict{Tuple{String,String},Float64}()
        for row in CSV.read(COMMITTED, DataFrame) |> eachrow
            ref[(String(row.substrate), String(row.tf))] = Float64(row.alpha)
        end
        @test Set(keys(mine)) == Set(keys(ref))
        maxdiff = maximum(abs(mine[k] - ref[k]) for k in keys(ref))
        @test maxdiff < 1e-6
        @info "max |Δalpha| vs committed handoff" maxdiff
    end

    # multiplier property tests use the fitted edges + ratios from the expression
    fitted = Dict{String,Vector{Tuple{String,Float64}}}()
    for row in CSV.read(FITTED, DataFrame) |> eachrow
        push!(get!(fitted, String(row.substrate), Tuple{String,Float64}[]),
              (String(row.tf), Float64(row.alpha)))
    end
    tfmean = Dict(g => mean(expr[gidx[g], :]) for g in genes)
    T = size(expr, 2)
    ratios(sub, k) = Float64[expr[gidx[tf], k] / tfmean[tf] for (tf, _) in fitted[sub]]
    alphas(sub) = Float64[a for (_, a) in fitted[sub]]
    Mseries(sub; q=1.0) = [multiplier(alphas(sub), ratios(sub, k); q=q) for k in 1:T]

    @testset "multiplier is 1 at the mean" begin
        for sub in ("CLN2", "VHR1", "RTS3", "HIS1")
            ones_ratio = fill(1.0, length(alphas(sub)))
            @test isapprox(multiplier(alphas(sub), ones_ratio), 1.0; atol=1e-12)
        end
    end

    @testset "mean-preserving over the cell cycle" begin
        for sub in ("CLN2", "VHR1", "RTS3", "YGP1")
            @test isapprox(mean(Mseries(sub)), 1.0; atol=0.06)
        end
    end

    @testset "committed TF means make <M> = 1 to machine precision" begin
        # Regression guard: tf_means.csv must hold the discrete mean of the SAME
        # NaN-interpolated 22-point trajectory the multiplier rides on, so the
        # deviation form is mean-preserving (<M_x>_t = 1) for EVERY gene, not just
        # the ones with complete TF data. Averaging TF means over only the non-NaN
        # timepoints (the prior bug) left <M>-1 up to 0.112 for ~394 genes.
        committed_means = Dict{String,Float64}()
        for row in CSV.read(joinpath(@__DIR__, "data", "tf_means.csv"), DataFrame) |> eachrow
            ismissing(row.tf_mean_rpm) || (committed_means[String(row.tf)] = Float64(row.tf_mean_rpm))
        end
        # ratios use the file's means (in RPM) against the interpolated trajectory.
        cratios(sub, k) = Float64[(expr[gidx[tf], k] / RPM_TO_MOLECULES) / committed_means[tf]
                                  for (tf, _) in fitted[sub]]
        cMseries(sub) = [multiplier(alphas(sub), cratios(sub, k)) for k in 1:T]
        maxdev = 0.0
        for sub in keys(fitted)
            all(haskey(committed_means, tf) for (tf, _) in fitted[sub]) || continue
            maxdev = max(maxdev, abs(mean(cMseries(sub)) - 1.0))
        end
        @info "max |<M>-1| over all genes (committed tf_means.csv)" maxdev
        @test maxdev < 1e-9
    end

    @testset "deviation form stays non-negative" begin
        for sub in ("VHR1", "RTS3", "YGP1", "LAP3", "CAR1")
            @test minimum(Mseries(sub)) >= -1e-9
        end
    end

    @testset "q scales amplitude, preserves mean" begin
        sub = "RTS3"
        @test all(isapprox.(Mseries(sub; q=0.0), 1.0; atol=1e-12))
        @test isapprox(mean(Mseries(sub; q=0.5)), 1.0; atol=0.06)
    end

    @testset "helper unit + edge cases" begin
        # linear_interp: clamp at edges, linear in the middle
        xp = [0.0, 10.0, 20.0]; fp = [1.0, 3.0, 2.0]
        @test linear_interp(-5, xp, fp) == 1.0          # below range -> first
        @test linear_interp(25, xp, fp) == 2.0          # above range -> last
        @test isapprox(linear_interp(5, xp, fp), 2.0)   # halfway 1->3
        # build_regulators: measured-measured only, dedup + sort
        e = [("A","X"), ("B","X"), ("A","X"), ("C","Y")]
        reg2 = build_regulators(e, Set(["A","B","X"]))
        @test reg2["X"] == ["A","B"]                    # deduped + sorted
        @test !haskey(reg2, "Y")                        # Y not measured -> dropped
        # multiplier: at the mean -> 1; no regulators -> 1; signed deviation
        @test multiplier([1.0, 1.0], [1.0, 1.0]) == 1.0
        @test multiplier(Float64[], Float64[]) == 1.0
        @test isapprox(multiplier([2.0, -1.0], [1.5, 1.0]), 1 + 1.0/3.0)
        @test isapprox(multiplier([2.0, -1.0], [1.5, 1.0]; q=0.0), 1.0)
        # loaded expression is finite and rescaled (NaN-interpolated, RPM->molecules)
        @test all(isfinite, expr)
    end

    # ----------------------------------------------------------------------
    # Joint multi-dataset headline anchor. Skips gracefully when the fetched
    # gold standard / name map (data/external/) are absent, so CI stays green
    # without third-party data. Run `julia transcription_multiplier/fetch_datasets.jl
    # --gold-only` to enable it. Scores the SHIPPED joint coupling CSVs (no LP
    # re-solve), so this is fast.
    @testset "joint fit reproduces headline AUROC + DeLong (needs data/external/)" begin
        ext  = joinpath(@__DIR__, "data", "external")
        need = [joinpath(ext, "SGD_features.tab"),
                joinpath(ext, "inferelator_gold_standard.tsv"),
                joinpath(@__DIR__, "data", "tf_network_fitted.csv"),
                joinpath(@__DIR__, "data", "joint_multidataset_alpha.csv"),
                joinpath(@__DIR__, "data", "loo_kelliher_alpha.csv")]
        if all(isfile, need)
            include(joinpath(@__DIR__, "joint_score.jl"))
            aucT, aucJ, diff, lo, hi, ok = Main.main(; refit = false)
            @test isapprox(aucT, 0.804; atol = 0.01)   # single-dataset (Teufel)
            @test isapprox(aucJ, 0.899; atol = 0.01)   # joint multi-dataset
            @test isapprox(diff, 0.071; atol = 0.01)   # paired DeLong delta
            @test lo > 0                               # CI excludes zero
            @test isapprox(lo, 0.052; atol = 0.01)
            @test isapprox(hi, 0.091; atol = 0.01)
            @test ok
        else
            absent = join([basename(f) for f in need if !isfile(f)], ", ")
            @info "skipping joint headline test; missing $absent. " *
                  "Run transcription_multiplier/fetch_datasets.jl --gold-only to enable."
            @test_skip true
        end
    end
end
