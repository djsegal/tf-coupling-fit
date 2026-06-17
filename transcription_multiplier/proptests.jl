#!/usr/bin/env julia
# Property / invariant tests for the cell-cycle transcription multiplier.
#
#   julia --project=. transcription_multiplier/proptests.jl        # from repo root
#   julia --project=. proptests.jl                                 # from package dir
#
# ADDITIVE: this file is independent of runtests.jl and does not modify it.
# It exercises the multiplier's *stated* invariants on randomized inputs:
#   (a) <M>_t == 1 to machine precision, for any alpha signs and any q_x
#   (b) M >= 0 (non-negativity) over a grid (clamped Dict method)
#   (c) monotonicity in q_x: larger q_x => larger deviation amplitude |M-1|
#   (d) numerical edge cases: all-zero alpha, single edge, repressor-only
#
# Prints a per-property PASS/FAIL summary and an "N/M properties passed" line.

using Test, Statistics, Random

# Prefer the module if it is loadable; otherwise include the canonical script.
# Either way we use the same `multiplier` methods (vector + Dict).
const _mult = let
    try
        @eval using TranscriptionMultiplier
        TranscriptionMultiplier.multiplier
    catch
        include(joinpath(@__DIR__, "multiplier.jl"))   # Dict method + load_handoff
        include(joinpath(@__DIR__, "refit.jl"))        # vector method (multiplier(alphas, ratios))
        multiplier
    end
end

Random.seed!(0xC0FFEE)

# --- helpers ---------------------------------------------------------------

# Build a random regulator set: k edges, alpha drawn from a sign-mixed dist.
function random_regs(k; allow_zero_sum = false)
    while true
        a = randn(k) .* (1.0 .+ 2 .* rand(k))      # mixed magnitudes & signs
        if !allow_zero_sum && sum(abs, a) > 1e-8
            return a
        elseif allow_zero_sum
            return a
        end
    end
end

# A random "trajectory" of TF levels for k TFs over T timepoints, plus the
# per-TF means computed FROM that trajectory (so ratios average to 1 per TF).
function random_trajectory(k, T)
    levels = abs.(randn(k, T)) .+ 0.1            # strictly positive TF amounts
    means  = vec(mean(levels; dims = 2))         # discrete mean per TF
    return levels, means
end

# Low-level vector multiplier value at a single timepoint (no clamp), via the
# vector method: multiplier(alphas::Vector, ratios::Vector; q).
Mvec(alphas, ratios; q = 1.0) = _mult(alphas, ratios; q = q)

results = Dict{String,Bool}()

@testset "multiplier invariants (property-based)" begin

    # ---------------------------------------------------------------- (a)
    # <M>_t == 1 to machine precision for ANY alpha signs and ANY q.
    # The deviation form is linear in ratios, so if each TF's ratio averages to
    # 1 across timepoints (means computed from the trajectory), <M>_t == 1.
    @testset "(a) mean-preserving <M>_t == 1" begin
        worst = 0.0
        for trial in 1:2000
            k = rand(1:8); T = rand(2:30)
            alphas = random_regs(k)
            levels, means = random_trajectory(k, T)
            q = rand() < 0.1 ? rand((-3.0, 0.0, 1e6)) : (2 * rand() - 1) * 5  # incl. weird q
            Ms = [Mvec(alphas, levels[:, t] ./ means; q = q) for t in 1:T]
            worst = max(worst, abs(mean(Ms) - 1.0))
        end
        @info "(a) worst |<M> - 1|" worst
        ok = worst < 1e-9
        results["(a) <M>_t == 1"] = ok
        @test ok
    end

    # ---------------------------------------------------------------- (b)
    # M >= 0 over the grid. The Dict-based public method clamps to max(0, M),
    # so non-negativity must hold for arbitrary (even pathological) inputs.
    @testset "(b) non-negativity (clamped Dict method)" begin
        worstmin = Inf
        for trial in 1:2000
            k = rand(1:8); T = rand(2:30)
            alphas = random_regs(k)
            levels, means = random_trajectory(k, T)
            tfs = ["TF$j" for j in 1:k]
            edges = Dict("S" => [(tfs[j], alphas[j]) for j in 1:k])
            meansd = Dict(tfs[j] => means[j] for j in 1:k)
            q = rand() * 10            # large q stresses the floor
            for t in 1:T
                lv = Dict(tfs[j] => levels[j, t] for j in 1:k)
                M = _mult("S", lv, edges, meansd; q = q, clamp = true)
                worstmin = min(worstmin, M)
            end
        end
        @info "(b) min clamped M over grid" worstmin
        ok = worstmin >= 0.0
        results["(b) M >= 0"] = ok
        @test ok
    end

    # ---------------------------------------------------------------- (c)
    # Monotonicity in q: larger q => larger deviation amplitude |M - 1|,
    # at every timepoint (using the unclamped vector form, which is exactly
    # linear in q: M - 1 = q * dev/sabs, so |M-1| is non-decreasing in q).
    @testset "(c) monotone amplitude in q_x" begin
        violations = 0
        for trial in 1:2000
            k = rand(1:8); T = rand(2:20)
            alphas = random_regs(k)
            levels, means = random_trajectory(k, T)
            qa, qb = sort(rand(2) .* 5)        # 0 <= qa <= qb
            for t in 1:T
                r = levels[:, t] ./ means
                amp_a = abs(Mvec(alphas, r; q = qa) - 1.0)
                amp_b = abs(Mvec(alphas, r; q = qb) - 1.0)
                amp_b + 1e-12 < amp_a && (violations += 1)
            end
        end
        @info "(c) monotonicity violations" violations
        ok = violations == 0
        results["(c) monotone in q_x"] = ok
        @test ok
    end

    # ---------------------------------------------------------------- (d)
    # Numerical edge cases.
    @testset "(d) edge cases" begin
        edgeok = true

        # all-zero alpha (vector form) -> M == 1 exactly (no regulation)
        v1 = Mvec([0.0, 0.0, 0.0], [5.0, 0.1, 2.0]; q = 3.0)
        @test v1 == 1.0
        edgeok &= (v1 == 1.0)

        # all-zero alpha (Dict form) -> M == 1 (sabs == 0 short-circuit)
        edges0 = Dict("S" => [("A", 0.0), ("B", 0.0)])
        means0 = Dict("A" => 2.0, "B" => 3.0)
        lv0 = Dict("A" => 10.0, "B" => 0.5)
        v2 = _mult("S", lv0, edges0, means0; q = 1.0)
        @test v2 == 1.0
        edgeok &= (v2 == 1.0)

        # empty edge list (no regulators) -> 1.0
        v3 = Mvec(Float64[], Float64[]; q = 1.0)
        @test v3 == 1.0
        edgeok &= (v3 == 1.0)

        # single activator edge: M = 1 + q*(r-1) since sabs == |alpha|
        for _ in 1:50
            a = abs(randn()) + 0.5; r = abs(randn()) + 0.1; q = rand() * 3
            expected = 1.0 + q * (r - 1.0)
            got = Mvec([a], [r]; q = q)
            ok = isapprox(got, expected; atol = 1e-12)
            @test ok
            edgeok &= ok
        end

        # repressor-only (all alpha < 0): higher TF (r>1) must DECREASE M.
        for _ in 1:50
            k = rand(1:5)
            a = -(abs.(randn(k)) .+ 0.5)        # strictly negative
            r = (abs.(randn(k)) .+ 1.05)        # all ratios > 1 (TF up)
            M = Mvec(a, r; q = 1.0)
            ok = M < 1.0
            @test ok
            edgeok &= ok
        end

        # single repressor at the mean (r==1) -> M == 1 regardless of q
        for _ in 1:50
            a = -(abs(randn()) + 0.5); q = rand() * 5
            got = Mvec([a], [1.0]; q = q)
            ok = isapprox(got, 1.0; atol = 1e-12)
            @test ok
            edgeok &= ok
        end

        results["(d) edge cases"] = edgeok
    end
end

# --- summary ---------------------------------------------------------------
println("\n--- multiplier invariant summary ---")
let npass = 0
    for k in sort(collect(keys(results)))
        status = results[k] ? "PASS" : "FAIL"
        println("  $status  $k")
        results[k] && (npass += 1)
    end
    println("$(npass)/$(length(results)) properties passed")
end
