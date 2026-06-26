#!/usr/bin/env julia
# EE6: type-stability / allocation spot check on the multiplier hot path.
#
#   julia --project=. perf_check.jl    # from repo root
#   julia --project=. perf_check.jl                             # from package dir
#
# ADDITIVE: read-only diagnostics; touches nothing else. Captures whether the
# core `multiplier` methods are type-stable (@code_warntype) and how much a
# representative call allocates (@allocated / @time).

using InteractiveUtils

const _loaded_module = try
    @eval using TranscriptionMultiplier
    true
catch
    include(joinpath(@__DIR__, "multiplier.jl"))
    include(joinpath(@__DIR__, "refit.jl"))
    false
end
const M = _loaded_module ? TranscriptionMultiplier.multiplier : multiplier

# ---------------------------------------------------------------- inputs ----
# Representative vector-form call (the inner loop in joint_score / Mseries):
alphas  = [2.0, -1.0, 0.5, 1.25]
ratios  = [1.1, 0.9, 1.3, 0.8]

# Representative Dict-form public call:
edges = Dict("CLN2" => [("SWI4", 2.0), ("MBP1", -1.0), ("SWI6", 0.5)])
means = Dict("SWI4" => 100.0, "MBP1" => 50.0, "SWI6" => 75.0)
levels = Dict("SWI4" => 120.0, "MBP1" => 40.0, "SWI6" => 80.0)

# warm up (compile)
M(alphas, ratios; q = 1.0)
M("CLN2", levels, edges, means; q = 1.0)

# ---------------------------------------------------- type stability --------
# Capture @code_warntype output as text and scan for instability markers
# (`::Any`, `::Union{...}`, red `Body` type). For these closed-form methods we
# expect a concrete `Float64` return and no boxed temporaries.
function warntype_text(f, argtypes)
    io = IOBuffer()
    code_warntype(io, f, argtypes)
    return String(take!(io))
end

function classify(label, txt)
    # A type-stable body has a concrete inferred return (Body::Float64) and no
    # `::Any` / abstract `Union` temporaries flagged by warntype.
    body_line = ""
    for ln in split(txt, '\n')
        if occursin("Body::", ln)
            body_line = strip(ln); break
        end
    end
    has_any   = occursin("::Any", txt)
    has_union = occursin(r"::Union\{", txt)
    stable    = occursin("Body::Float64", txt) && !has_any && !has_union
    println("  [$label] Body: ", isempty(body_line) ? "(not found)" : body_line)
    println("  [$label] contains ::Any   : ", has_any)
    println("  [$label] contains ::Union : ", has_union)
    println("  [$label] => ", stable ? "TYPE-STABLE" : "NOT type-stable (inspect below)")
    return stable
end

println("=== @code_warntype: vector method multiplier(::Vector, ::Vector; q) ===")
tv = warntype_text(M, (Vector{Float64}, Vector{Float64}))
stable_vec = classify("vector", tv)

println("\n=== @code_warntype: Dict method multiplier(::String, ::Dict, ::Dict, ::Dict; q) ===")
td = warntype_text(M, (String, typeof(levels), typeof(edges), typeof(means)))
stable_dict = classify("dict", td)

# ---------------------------------------------------- allocations -----------
# Steady-state allocations (exclude first-call compilation, already warmed).
println("\n=== allocations (@allocated, post-warmup) ===")
a_vec  = @allocated M(alphas, ratios; q = 1.0)
a_dict = @allocated M("CLN2", levels, edges, means; q = 1.0)
println("  vector method : $a_vec bytes / call")
println("  dict method   : $a_dict bytes / call")

println("\n=== @time of a representative call (vector method) ===")
@time M(alphas, ratios; q = 1.0)
println("=== @time of a representative call (dict method) ===")
@time M("CLN2", levels, edges, means; q = 1.0)

# Tight loop timing for a realistic batch (T=22 timepoints x many genes).
function loop_vec(alphas, ratios, n)
    s = 0.0
    for _ in 1:n
        s += M(alphas, ratios; q = 1.0)
    end
    return s
end
loop_vec(alphas, ratios, 1)   # warm
println("\n=== batch loop: 10_000 vector-method calls ===")
b = @allocated loop_vec(alphas, ratios, 10_000)
print("  @allocated over 10k calls: $b bytes; "); @time loop_vec(alphas, ratios, 10_000)

println("\n--- EE6 summary ---")
println("  vector method type-stable : ", stable_vec)
println("  dict method  type-stable  : ", stable_dict)
println("  vector method bytes/call  : ", a_vec)
println("  dict method  bytes/call   : ", a_dict)

# If anything is unstable, dump the full warntype text for inspection.
if !stable_vec
    println("\n[full @code_warntype, vector method]\n", tv)
end
if !stable_dict
    println("\n[full @code_warntype, dict method]\n", td)
end
