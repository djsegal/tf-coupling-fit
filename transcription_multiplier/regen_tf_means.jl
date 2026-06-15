#!/usr/bin/env julia
# Regenerate transcription_multiplier/data/tf_means.csv so each TF's mean is the
# discrete mean of the SAME NaN-interpolated 22-point trajectory the multiplier and
# refit evaluate on (refit.jl::load_rna_seq). This makes the deviation-form
# multiplier mean-preserving (<M_x>_t = 1) to machine precision.
#   julia --project=. transcription_multiplier/regen_tf_means.jl
using CSV, DataFrames
include(joinpath(@__DIR__, "refit.jl"))

const REPO = normpath(joinpath(@__DIR__, ".."))
const EXPR = joinpath(REPO, "data", "WT_unstressed_readspermillionreads.csv")
const MEANS = joinpath(@__DIR__, "data", "tf_means.csv")

# Interpolated, RPM->molecules expression on the 22-point grid (exact fit pipeline).
time_axis, genes, expr = load_rna_seq(EXPR)
gidx = Dict(g => i for (i, g) in enumerate(genes))

old = CSV.read(MEANS, DataFrame)
new = copy(old)
for (i, tf) in enumerate(old.tf)
    haskey(gidx, tf) || error("TF $tf not found in loaded expression")
    mean_molecules = sum(view(expr, gidx[tf], :)) / size(expr, 2)  # discrete mean (molecules)
    mean_rpm = mean_molecules / RPM_TO_MOLECULES
    # Full Float64 precision (not rounded to 6 digits) so the multiplier is
    # mean-preserving to machine precision when the CSV is read back.
    new.tf_mean_rpm[i] = mean_rpm
    new.tf_mean_molecules_per_cell[i] = mean_molecules
end
# n_timepoints_available / n_timepoints_total are raw-data metadata: leave untouched.
CSV.write(MEANS, new)
println("Rewrote $MEANS ($(nrow(new)) TFs) on the interpolated 22-point grid.")
