"""
Python adapter for the cell-cycle transcription multiplier.

This is a DOWNSTREAM-CONVENIENCE ADAPTER, not the project's implementation. The
canonical code is Julia (../multiplier.jl, ../refit.jl, with the test suite in
../runtests.jl); this file exists only so a group whose pipeline is in Python can
apply the released couplings without a Julia dependency. It reads the same released
CSVs and computes the identical deviation-form multiplier:

    M_x(t) = 1 + q_x * ( sum_i alpha_i * (TF_i(t)/TF_i_mean - 1) ) / sum_i |alpha_i|

Stdlib only (csv), no third-party dependencies. Units cancel in TF_i(t)/TF_i_mean,
so feed TF levels in whatever unit your model stores.
"""
import csv
import os
from collections import defaultdict

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def load_handoff(data_dir=_DATA):
    """Return (edges, means): edges[substrate] = list of (tf, alpha); means[tf] = float."""
    edges = defaultdict(list)
    with open(os.path.join(data_dir, "tf_network_fitted.csv"), newline="") as fh:
        for row in csv.DictReader(fh):
            edges[row["substrate"]].append((row["tf"], float(row["alpha"])))
    means = {}
    with open(os.path.join(data_dir, "tf_means.csv"), newline="") as fh:
        r = csv.DictReader(fh)
        mean_col = "tf_mean_rpm" if "tf_mean_rpm" in r.fieldnames else r.fieldnames[1]
        tf_col = r.fieldnames[0]
        for row in r:
            means[row[tf_col]] = float(row[mean_col])
    return dict(edges), means


def multiplier(substrate, tf_levels, edges, means, q=1.0, clamp=True):
    """Deviation-form rate multiplier M_x for `substrate` given current `tf_levels`.

    tf_levels: dict {tf_name: current_level}. Regulators absent from tf_levels or
    `means` are skipped. Returns 1.0 if the substrate has no usable regulator.
    `q` is the per-gene cell-cycle weight in [0,1] (use the released q_x, not 1).
    `clamp=True` floors the multiplier at 0 (a rate cannot be negative).
    """
    num = 0.0
    denom = 0.0
    for tf, alpha in edges.get(substrate, []):
        m = means.get(tf)
        lvl = tf_levels.get(tf)
        if m is None or lvl is None or m == 0.0:
            continue
        num += alpha * (lvl / m - 1.0)
        denom += abs(alpha)
    if denom == 0.0:
        return 1.0
    val = 1.0 + q * num / denom
    return max(0.0, val) if clamp else val


if __name__ == "__main__":
    edges, means = load_handoff()
    demo = {tf: means[tf] for tf in list(means)[:5]}  # all regulators at their mean
    sub = next(iter(edges))
    print(f"M({sub}) at all-TFs-at-mean = {multiplier(sub, demo, edges, means):.6f} (expect ~1.0)")
