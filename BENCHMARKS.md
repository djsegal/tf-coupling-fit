# Benchmarks

Baseline median timings for the transcription-multiplier hot paths, captured
with [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl). These
are the numbers a future change is compared against, so regressions show up.

Reproduce (standalone, from the package directory):

```
julia --project=bench bench/benchmarks.jl --write
```

The suite times three hot paths:

1. the vector-form `multiplier` call (the inner loop of the M-series / scoring),
2. the Dict-form `multiplier` public call, and
3. a single-substrate L1 fit via `fit_substrate` (JuMP/HiGHS).

## Baseline (median over the BenchmarkTools sample)

Captured on Julia 1.10.5.

| hot path | median | min | memory | allocations |
|---|---|---|---|---|
| vector multiplier (`multiplier(::Vector, ::Vector; q)`) | 6.4 ns | 5.5 ns | 0 B | 0 |
| Dict multiplier (`multiplier(substrate, levels, edges, means; q)`) | 174.7 ns | 151.7 ns | 0 B | 0 |
| single-substrate L1 fit (`fit_substrate`) | 721.333 us | 687.709 us | 269.84 KiB | 4751 |

Notes:

- `multiplier_vector` / `multiplier_dict` are closed-form and allocation-light;
  they are the calls executed millions of times when scoring trajectories.
- `fit_substrate` is dominated by JuMP model construction + the HiGHS solve, so
  its time/allocations are orders of magnitude larger and naturally noisier.
- Times depend on the machine; treat the committed numbers as a *relative*
  regression guard rather than an absolute spec.
