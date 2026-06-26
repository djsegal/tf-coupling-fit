# Python adapter (downstream convenience)

The canonical implementation of this project is **Julia** (see `../multiplier.jl`,
`../refit.jl`, and the test suite `../runtests.jl`). All fitting, validation, and the
released data are produced and checked in Julia.

This directory is a **thin Python adapter** provided only so a group whose pipeline is
in Python (for example, a downstream whole-cell-model effort) can apply the released
couplings without a Julia dependency. It is not the implementation and carries no
analysis logic beyond evaluating the published formula on the released CSVs.

- `multiplier.py`: stdlib-only (`csv`); `load_handoff()` + `multiplier(substrate, tf_levels, edges, means; q, clamp)`.

It computes the identical deviation-form multiplier as the Julia code:

```
M_x(t) = 1 + q_x * ( sum_i alpha_i * (TF_i(t)/TF_i_mean - 1) ) / sum_i |alpha_i|
```

Use the released per-gene `q_x` (from `../data/qx_scores.csv`), not `q = 1`. If you can
run Julia, prefer the canonical `multiplier.jl`.
