#!/usr/bin/env julia
# EE4: minimal Documenter.jl scaffold for TranscriptionMultiplier.
#
# Build (standalone) from the package directory:
#
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'   # once
#   julia --project=docs docs/make.jl
#
# Output goes to docs/build/. doctests in the docstrings are run as part of the
# build (doctest = true). This is a valid scaffold; you do not need to build the
# HTML to verify the docstrings exist.
#
# ADDITIVE: lives in its own environment (docs/Project.toml) and loads the
# package from the sibling project via LOAD_PATH; touches no existing files.

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using Documenter
using TranscriptionMultiplier

DocMeta.setdocmeta!(TranscriptionMultiplier, :DocTestSetup,
                    :(using TranscriptionMultiplier); recursive = true)

makedocs(
    sitename = "TranscriptionMultiplier.jl",
    modules  = [TranscriptionMultiplier],
    authors  = "Daniel J. Segal",
    pages    = ["Home" => "index.md"],
    # Run doctests during the build; fail on broken/missing cross-references.
    doctest  = true,
    checkdocs = :exports,
    format   = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
)
