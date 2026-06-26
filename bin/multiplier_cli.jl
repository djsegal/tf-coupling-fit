#!/usr/bin/env julia
# EE5: tiny CLI for the cell-cycle transcription-rate multiplier.
#
# Given a substrate gene and the current TF levels, prints the multiplier Mₓ and
# the gated rate  k_x(t) = k_base * Mₓ. Reads the fitted handoff (edges + TF
# means) and the per-gene q_x from data/, with informative input validation.
#
# Usage (from the package directory):
#
#   julia --project=. bin/multiplier_cli.jl --gene CLN2 --tf SWI4=120,MBP1=40,SWI6=80
#   julia --project=. bin/multiplier_cli.jl --gene CLN2 --tf 120,40,80           # positional, regulator order
#   julia --project=. bin/multiplier_cli.jl --gene CLN2 --tf SWI4=120 --kbase 0.5 --q 1.0
#   julia --project=. bin/multiplier_cli.jl --gene CLN2 --list-regulators
#   julia --project=. bin/multiplier_cli.jl --help
#
# ADDITIVE: depends only on the package's public API (load_handoff, multiplier)
# plus CSV/DataFrames for the q_x table; touches no existing files.

const PKG_DIR = normpath(joinpath(@__DIR__, ".."))

# Load the package (preferred) or fall back to including the scripts directly.
try
    @eval using TranscriptionMultiplier
catch
    include(joinpath(PKG_DIR, "multiplier.jl"))
    include(joinpath(PKG_DIR, "refit.jl"))
end
using CSV, DataFrames

const PROG = "multiplier_cli.jl"

usage() = """
$PROG - cell-cycle transcription-rate multiplier CLI

Required:
  --gene NAME             substrate gene (e.g. CLN2)
  --tf  SPEC              TF levels, either
                            named:      SWI4=120,MBP1=40,SWI6=80
                            positional: 120,40,80   (in the substrate's regulator order)

Optional:
  --q VALUE               override q_x (cycling amplitude in [0,1]); default: from data/qx_scores.csv
  --kbase VALUE           base rate k_base; default 1.0 (so gated rate == M when omitted)
  --no-clamp              do not apply the max(0, M) floor
  --datadir DIR           data directory (default: the package data/)
  --list-regulators       print the substrate's (tf, alpha) regulators and exit
  --help, -h              show this message

Prints Mₓ and the gated rate k_x = k_base * Mₓ.
"""

die(msg) = (println(stderr, "error: ", msg); println(stderr, "\nrun `$PROG --help` for usage."); exit(2))

# ---- argument parsing ------------------------------------------------------
function parse_args(argv)
    opts = Dict{String,String}()
    flags = Set{String}()
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a in ("--help", "-h")
            push!(flags, "help"); i += 1
        elseif a == "--no-clamp"
            push!(flags, "no-clamp"); i += 1
        elseif a == "--list-regulators"
            push!(flags, "list-regulators"); i += 1
        elseif startswith(a, "--")
            key = a[3:end]
            i + 1 <= length(argv) || die("option --$key requires a value")
            opts[key] = argv[i+1]; i += 2
        else
            die("unexpected argument: $a")
        end
    end
    return opts, flags
end

# Parse the --tf spec into either named pairs or a positional vector.
function parse_tf(spec::AbstractString)
    spec = strip(spec)
    isempty(spec) && die("--tf was empty; give e.g. SWI4=120,MBP1=40 or 120,40")
    parts = split(spec, ',')
    named = Dict{String,Float64}()
    positional = Float64[]
    saw_named = false; saw_pos = false
    for p in parts
        p = strip(p)
        isempty(p) && continue
        if occursin('=', p)
            saw_named = true
            k, v = split(p, '='; limit = 2)
            k = strip(k); v = strip(v)
            isempty(k) && die("malformed --tf entry '$p' (empty TF name)")
            val = tryparse(Float64, v)
            val === nothing && die("TF level for '$k' is not a number: '$v'")
            named[String(k)] = val
        else
            saw_pos = true
            val = tryparse(Float64, p)
            val === nothing && die("TF level '$p' is not a number")
            push!(positional, val)
        end
    end
    saw_named && saw_pos &&
        die("--tf mixes named (TF=value) and positional entries; use one form")
    return saw_named ? (:named, named) : (:positional, positional)
end

# ---- main ------------------------------------------------------------------
function main(argv)
    opts, flags = parse_args(argv)
    if "help" in flags || isempty(argv)
        print(usage()); return 0
    end

    haskey(opts, "gene") || die("missing required --gene")
    gene = strip(opts["gene"])
    isempty(gene) && die("--gene was empty")

    datadir = get(opts, "datadir", (@isdefined(DATA_DIR) ? DATA_DIR : joinpath(PKG_DIR, "data")))

    # load_handoff has its own informative guards for missing dir / files.
    edges, means = load_handoff(datadir)

    haskey(edges, gene) ||
        die("unknown gene '$gene': not a substrate in the fitted network " *
            "($(joinpath(datadir, "tf_network_fitted.csv"))). " *
            "Genes with no regulators have multiplier == 1 by definition.")
    regs = edges[gene]   # Vector{Tuple{String,Float64}}

    if "list-regulators" in flags
        println("regulators of $gene (tf, alpha):")
        for (tf, a) in regs
            haskey(means, tf) || (println("  $tf\t$a\t(WARNING: no mean in tf_means.csv)"); continue)
            println("  $tf\t$a")
        end
        return 0
    end

    haskey(opts, "tf") || die("missing required --tf (e.g. --tf SWI4=120,MBP1=40)")
    kind, tfin = parse_tf(opts["tf"])

    # Build the tf_levels Dict the public multiplier expects, validating shape.
    reg_tfs = [tf for (tf, _) in regs]
    tf_levels = Dict{String,Float64}()
    if kind == :named
        # Warn about names that aren't regulators of this gene (likely typos).
        extra = setdiff(keys(tfin), Set(reg_tfs))
        isempty(extra) ||
            println(stderr, "warning: --tf names not among $gene's regulators (ignored): " *
                            join(sort(collect(extra)), ", "))
        for (tf, v) in tfin
            tf in reg_tfs && (tf_levels[tf] = v)
        end
        missing_tfs = [tf for tf in reg_tfs if !haskey(tf_levels, tf)]
        isempty(missing_tfs) ||
            die("no level given for regulator(s) of $gene: " * join(missing_tfs, ", ") *
                ". Provide them (e.g. $(missing_tfs[1])=<value>) or use --list-regulators.")
    else
        length(tfin) == length(reg_tfs) ||
            die("--tf has $(length(tfin)) value(s) but $gene has $(length(reg_tfs)) regulator(s) " *
                "($(join(reg_tfs, ", "))). Give one value per regulator in that order, " *
                "or use named form (TF=value).")
        for (tf, v) in zip(reg_tfs, tfin)
            tf_levels[tf] = v
        end
    end

    # q_x: from --q override, else look up the per-gene value in qx_scores.csv.
    local q::Float64
    if haskey(opts, "q")
        qv = tryparse(Float64, opts["q"])
        qv === nothing && die("--q is not a number: '$(opts["q"])'")
        q = qv
    else
        qx_path = joinpath(datadir, "qx_scores.csv")
        isfile(qx_path) ||
            die("missing q_x: $qx_path not found and no --q given. " *
                "Pass --q <value in [0,1]> explicitly.")
        qdf = CSV.read(qx_path, DataFrame)
        ("substrate" in names(qdf) && "q_x" in names(qdf)) ||
            die("$qx_path lacks 'substrate'/'q_x' columns; pass --q explicitly.")
        row = findfirst(==(gene), qdf.substrate)
        row === nothing &&
            die("missing q_x for gene '$gene' in $qx_path. Pass --q <value> explicitly.")
        qval = qdf.q_x[row]
        (ismissing(qval) || (qval isa Number && isnan(qval))) &&
            die("q_x for '$gene' is missing/NaN in $qx_path. Pass --q <value> explicitly.")
        q = Float64(qval)
    end

    do_clamp = !("no-clamp" in flags)
    kbase = 1.0
    if haskey(opts, "kbase")
        kv = tryparse(Float64, opts["kbase"])
        kv === nothing && die("--kbase is not a number: '$(opts["kbase"])'")
        kbase = kv
    end

    M = multiplier(gene, tf_levels, edges, means; q = q, clamp = do_clamp)
    rate = kbase * M

    println("gene        : $gene")
    println("regulators  : ", join([string(tf, "(", a, ")") for (tf, a) in regs], ", "))
    println("tf_levels   : ", join([string(tf, "=", tf_levels[tf]) for (tf, _) in regs], ", "))
    println("q_x         : $q", haskey(opts, "q") ? "  (from --q)" : "  (from qx_scores.csv)")
    println("clamp       : $do_clamp")
    println("k_base      : $kbase")
    println("------------------------------------------")
    println("M_x         : $M")
    println("gated rate  : $rate   (= k_base * M_x)")
    return 0
end

exit(main(ARGS))
