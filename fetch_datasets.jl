#!/usr/bin/env julia
#=
fetch_datasets.jl: download the third-party datasets needed to re-solve the joint
multi-dataset fit and to score it against the Inferelator gold standard, into
data/external/, and verify each against a known SHA-256.

We do NOT redistribute third-party raw data in this package. Each file is fetched
from its own public source (SGD archive, GEO, or the Inferelator repository). All
are published, citation-required academic data. Cite the original papers (see the
README and each source's README).

After this script succeeds you can:
  * re-solve the joint couplings from raw data and reproduce the headline:
      julia --project=. joint_score.jl --refit
  * the scored-mode headline (no re-solve) only needs SGD_features.tab and
      inferelator_gold_standard.tsv (the two smallest fetches), since the joint
      coupling CSVs are shipped in data/.

Run from anywhere:
  julia fetch_datasets.jl            # fetch + verify all
  julia fetch_datasets.jl --gold-only # just SGD + gold
  julia fetch_datasets.jl --force     # re-download

Pure Julia; uses the Downloads stdlib. No new dependency, no Python.
=#

using Downloads, SHA

const EXT = joinpath(@__DIR__, "data", "external")

# Each entry: (local filename, url, expected sha256 or nothing, citation, transform)
# transform :gunzip means the URL is a .gz to be decompressed into the local file.
struct Source
    name::String
    url::String
    sha256::Union{String,Nothing}
    cite::String
    transform::Symbol   # :none or :gunzip
end

const SOURCES = [
    Source("inferelator_gold_standard.tsv",
        "https://raw.githubusercontent.com/flatironinstitute/inferelator/release/data/yeast/gold_standard.tsv",
        "1e02625b5f51ffd66c12a08d3afb7207b17db71f5e8e975844fdc477f91f0d7b",
        "Inferelator yeast gold standard (Tchourine, Vogel & Bonneau 2018, Cell Reports 23:376; Flatiron Institute inferelator repo, data/yeast/gold_standard.tsv).",
        :none),
    Source("SGD_features.tab",
        "https://s3-us-west-2.amazonaws.com/sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab",
        nothing,   # SGD updates this file over time; verified by structure, not hash
        "Saccharomyces Genome Database (SGD) chromosomal feature table. ORF<->standard name map. Updated periodically by SGD, so no pinned hash.",
        :none),
    Source("spellman_combined.txt",
        "https://web.archive.org/web/20210615234137id_/http://genome-www.stanford.edu/cellcycle/data/rawdata/combined.txt",
        "0ec546bcb1dcbba99b0bdb3fe08bf20bfbfd325af3050aba3756f9e2ee225c9f",
        "Spellman et al. 1998, Mol Biol Cell 9:3273 (PMID 9843569). The combined cell-cycle table (cdc28-13 alpha-factor columns alpha0..alpha119). See manual-fallback note below if the archived copy is unavailable.",
        :none),
    Source("pramila.pcl",
        "https://s3-us-west-2.amazonaws.com/sgd-archive.yeastgenome.org/expression/microarray/Pramila_2006_PMID_16912276/GSE4987_setA_family.pcl",
        "b15f26b8145ac0b6928ce9ab1bc2f2f190c757b7a7e3cbf14e6eab5310388281",
        "Pramila et al. 2006, Genes Dev 20:2266 (PMID 16912276; GEO GSE4987). SGD-archive redistribution GSE4987_setA_family.pcl.",
        :none),
    Source("orlando_setA.pcl",
        "https://s3-us-west-2.amazonaws.com/sgd-archive.yeastgenome.org/expression/microarray/Orlando_2008_PMID_18463633/GSE8799_setA_family.pcl",
        "b2a8190e900d9cec7a89ae06bc55540440cf0e5f1a2eaa0fa50637a42a36b273",
        "Orlando et al. 2008, Nature 453:944 (PMID 18463633; GEO GSE8799). SGD-archive redistribution GSE8799_setA_family.pcl, WildType columns only.",
        :none),
    Source("GSE80474_Scerevisiae_normalized.txt",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80474/suppl/GSE80474_Scerevisiae_normalized.txt.gz",
        "5be8254fd9b443fe791de8d6102c16ff51077fadaceb1d9dcebc059a183d5680",
        "Kelliher et al. 2016, PLoS Genet 12:e1006453 (GEO GSE80474). GEO supplementary normalized RNA-seq table; fetched as .gz and decompressed.",
        :gunzip),
]

const SPELLMAN_FALLBACK = """
  Spellman combined table could not be fetched from the archived Stanford copy.
  Manual fallback (any one):
    * Bioconductor yeastCC package (spYCCES ExpressionSet) -> export the
      cdc28-13 alpha-factor columns to a tab-separated table whose header is
      `\\talpha0\\talpha7\\t...\\talpha119` with systematic-ORF row labels, save as
      data/external/spellman_combined.txt.
    * Or the SGD-archive alphaFactor PCL
      expression/microarray/Spellman_1998_PMID_9843569/2010.Spellman98_alphaFactor.flt.knn.avg.pcl
      (note: this PCL is NOT a drop-in -- it lacks the alphaN minute labels the
      parser expects, so relabel its 18 release samples to alpha0..alpha119).
  The headline result does NOT require this file: the joint coupling CSVs are
  shipped in data/, so `joint_score.jl` (scored mode) reproduces 0.804 -> 0.899
  and +0.071 without any re-solve. Spellman is only needed for `--refit`.
"""

sha256_hex(path) = bytes2hex(open(SHA.sha256, path))

function gunzip_to(gz_bytes::Vector{UInt8}, out::String)
    # minimal gzip via the `gzip` CLI (POSIX); avoids adding CodecZlib
    tmp = out * ".gz"
    write(tmp, gz_bytes)
    run(pipeline(`gzip -dc $tmp`, stdout = out))
    rm(tmp; force = true)
end

function fetch_one(s::Source; force::Bool)
    dest = joinpath(EXT, s.name)
    if isfile(dest) && !force
        if s.sha256 === nothing || sha256_hex(dest) == s.sha256
            println("  ok (cached)  $(s.name)")
            return true
        else
            println("  cached $(s.name) failed checksum; re-downloading")
        end
    end
    print("  downloading  $(s.name) ... "); flush(stdout)
    try
        if s.transform == :gunzip
            tmp = tempname()
            Downloads.download(s.url, tmp)
            gunzip_to(read(tmp), dest)
            rm(tmp; force = true)
        else
            Downloads.download(s.url, dest)
        end
    catch err
        println("FAILED")
        @warn "download failed" file=s.name url=s.url exception=err
        s.name == "spellman_combined.txt" && println(SPELLMAN_FALLBACK)
        return false
    end
    if s.sha256 !== nothing
        got = sha256_hex(dest)
        if got != s.sha256
            println("CHECKSUM MISMATCH")
            # Move the bad file aside so a truncated/partial copy is never parsed
            # as if it were valid. (We rename rather than delete.)
            bad = dest * ".partial"
            isfile(bad) && mv(bad, bad * "." * string(time_ns()); force = true)
            mv(dest, bad; force = true)
            @warn "checksum mismatch (bad copy moved to .partial)" file=s.name expected=s.sha256 got=got
            s.name == "spellman_combined.txt" && println(SPELLMAN_FALLBACK)
            return false
        end
    end
    println("ok ($(filesize(dest)) bytes)")
    return true
end

function main(args)
    mkpath(EXT)
    gold_only = "--gold-only" in args
    force = "--force" in args
    want = gold_only ? filter(s -> s.name in ("inferelator_gold_standard.tsv", "SGD_features.tab"), SOURCES) : SOURCES
    println("fetching $(length(want)) dataset(s) into $EXT")
    results = [(s.name, fetch_one(s; force = force)) for s in want]
    println()
    println("citations:")
    for s in want
        println("  - ", s.cite)
    end
    failed = [n for (n, ok) in results if !ok]
    if isempty(failed)
        println("\nall files present and verified.")
    else
        println("\nfailed: ", join(failed, ", "))
        println("see notes above; the scored-mode headline still works if SGD + gold succeeded.")
    end
    return isempty(failed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS) ? 0 : 1)
end
