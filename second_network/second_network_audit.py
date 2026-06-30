#!/usr/bin/env python3
"""
SECOND-NETWORK demonstration of the prior-signs-only control vs the
expression-fitted coupling network (referee #1: "shown on only ONE network").

Network 1 (in the paper): yeast cell-cycle / Teufel, scored vs the Inferelator gold.
Network 2 (here):         DREAM5 Network 3 (E. coli), scored vs the DREAM5 gold.

Same audit logic as TranscriptionMultiplier.jl / reply/workshop/julia_stats
(score_gold + ranking_metrics): a TF x target rectangle is scored as a ranked
classifier against the gold standard.  Each cell carries the network's edge
magnitude -- |coef| for the fit, unit magnitude for the prior-signs control --
and lasso-zeroed candidates fall into the tied-zero block exactly as in the paper.
AUROC/AUPRC are the trapezoid forms of JuliaStats.ranking_metrics, cross-checked
against the tie-aware (Mann-Whitney) sklearn values (they agree within ~0.02,
as the paper documents).

DATA (all public; URLs below). Re-fetch with  python3 second_network_audit.py --fetch
  candidate prior (signed):  RegulonDB network_tf_gene  (curated TF->gene, +/- effect)
      https://raw.githubusercontent.com/Natpod/Synthetic-Biology-Network-Analysis-on-Regulon-DB-E-coli-Data/master/network_tf_gene.txt
      (a verbatim RegulonDB release file; ISO-8859-1; cols: TFid,TFname,geneid,genename,effect,evidence,confidence)
  expression to fit:         DREAM5 net3 compendium (805 microarray chips x 4511 genes)
  TF list / gene-id map:     net3_transcription_factors.tsv / net3_gene_ids.tsv
  independent gold standard: DREAM5 net3 gold (RegulonDB-derived; 2066 positive edges)
      DREAM5 files mirrored from Synapse syn2787209 at:
      https://github.com/DuaaAlawad/AGRN/tree/master/Dataset/DREAM5
      (training data/Network 3 - E. coli/*  and  test data/DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv)

The candidate prior (current RegulonDB) and the gold (DREAM5's RegulonDB-derived
benchmark) share curation lineage -- the exact scenario the yeast paper audits.

Set DATA_DIR to wherever the files live (default: alongside this script under data/).
"""
import os, sys, subprocess, urllib.request, urllib.parse
import numpy as np, pandas as pd
from sklearn.linear_model import LassoLarsIC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score

DATA_DIR = os.environ.get("DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))

AGRN = "https://raw.githubusercontent.com/DuaaAlawad/AGRN/master/Dataset/DREAM5"
FILES = {
    "net3_expression_data.tsv":        f"{AGRN}/training data/Network 3 - E. coli/net3_expression_data.tsv",
    "net3_transcription_factors.tsv":  f"{AGRN}/training data/Network 3 - E. coli/net3_transcription_factors.tsv",
    "net3_gene_ids.tsv":               f"{AGRN}/training data/Network 3 - E. coli/net3_gene_ids.tsv",
    "net3_gold.tsv":                   f"{AGRN}/test data/DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv",
    "regulondb_tf_gene.txt":           "https://raw.githubusercontent.com/Natpod/Synthetic-Biology-Network-Analysis-on-Regulon-DB-E-coli-Data/master/network_tf_gene.txt",
}

def fetch():
    os.makedirs(DATA_DIR, exist_ok=True)
    for fn, url in FILES.items():
        dest = os.path.join(DATA_DIR, fn)
        if os.path.exists(dest):
            continue
        print(f"fetching {fn} ...", file=sys.stderr)
        # percent-encode the path (the GitHub raw paths contain spaces) but keep
        # the URL delimiters; do NOT leave the space in `safe`, or urllib rejects it.
        req = urllib.request.Request(urllib.parse.quote(url, safe=":/?=&-."),
                                     headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as f:
            f.write(r.read())

# ---------------------------------------------------------------- ranking metrics
def trapz(y, x):
    return sum((x[i]-x[i-1])*(y[i]+y[i-1])/2.0 for i in range(1, len(x)))

def ranking_metrics(scores, labels):
    """Trapezoid AUROC + AUPRC, identical to JuliaStats.ranking_metrics
    (stable descending sort by (-score, index); prepend origin / prec[0])."""
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    n = len(scores)
    order = sorted(range(n), key=lambda i: (-scores[i], i))
    y = labels[order]; P = y.sum(); N = n - P
    if P == 0 or N == 0:
        return float('nan'), float('nan')
    tp = np.cumsum(y); fp = np.cumsum(1-y)
    tpr = tp/P; fpr = fp/N; prec = tp/(tp+fp); rec = tpr
    roc = trapz(np.r_[0.0, tpr], np.r_[0.0, fpr])
    pr  = trapz(np.r_[prec[0], prec], np.r_[0.0, rec])
    return roc, pr

# ---------------------------------------------------------------- load
def main():
    gid = pd.read_csv(f"{DATA_DIR}/net3_gene_ids.tsv", sep="\t"); gid.columns = ["gid", "name"]
    name2gid = {}
    for _, r in gid.iterrows():
        name2gid.setdefault(str(r["name"]).lower(), r["gid"])     # case-insensitive

    gold = pd.read_csv(f"{DATA_DIR}/net3_gold.tsv", sep="\t", header=None, names=["tf", "tgt", "lab"])
    gold_pos = set(zip(gold.loc[gold.lab == 1, "tf"], gold.loc[gold.lab == 1, "tgt"]))
    gold_tfs = set(gold["tf"]); gold_tgts = set(gold["tgt"])

    expr = pd.read_csv(f"{DATA_DIR}/net3_expression_data.tsv", sep="\t")
    X = {g: expr[g].to_numpy(float) for g in expr.columns}

    reg = pd.read_csv(f"{DATA_DIR}/regulondb_tf_gene.txt", sep="\t", comment="#",
                      header=None, engine="python", encoding="latin-1").dropna(subset=[1, 3])
    cand = {}
    for _, r in reg.iterrows():
        t = name2gid.get(str(r[1]).lower()); g = name2gid.get(str(r[3]).lower())
        if t and g and t != g:
            cand.setdefault((t, g), set()).add(str(r[4]).strip())
    sign_of = lambda s: (+1 if s == {"+"} else -1 if s == {"-"} else 0)
    cand_sign = {e: sign_of(s) for e, s in cand.items()}
    cand_edges = set(cand); cand_tfs = {t for t, _ in cand_edges}; cand_tgts = {g for _, g in cand_edges}
    print(f"RegulonDB candidate prior: {len(cand_edges)} edges, {len(cand_tfs)} TFs, {len(cand_tgts)} targets")
    print(f"DREAM5 gold: {len(gold_pos)} positives; {len(gold_pos & cand_edges)} are candidates "
          f"({100*len(gold_pos & cand_edges)/len(gold_pos):.0f}% of gold)")

    # ---- per-target L1 fit (LassoLarsIC/BIC; single-regulator -> std OLS slope)
    regs_of = {}
    for (t, g) in cand_edges:
        regs_of.setdefault(g, []).append(t)
    coef = {}
    for g, tfs in regs_of.items():
        if g not in X: continue
        tfs_m = [t for t in tfs if t in X]
        if not tfs_m: continue
        Xs = StandardScaler().fit_transform(np.column_stack([X[t] for t in tfs_m]))
        ys = X[g] - X[g].mean()
        if Xs.shape[1] == 1:
            betas = np.array([np.polyfit(Xs[:, 0], ys, 1)[0]])
        else:
            betas = LassoLarsIC(criterion="bic").fit(Xs, ys).coef_
        for t, b in zip(tfs_m, betas):
            coef[(t, g)] = float(b)
    corr = {}
    for (t, g) in coef:
        a, b = X[t]-X[t].mean(), X[g]-X[g].mean()
        d = np.sqrt((a*a).sum()*(b*b).sum()); corr[(t, g)] = float((a*b).sum()/d) if d > 0 else 0.0

    fit_mag = {e: abs(v) for e, v in coef.items()}
    fit_nz  = {e: v for e, v in fit_mag.items() if v > 0}
    corr_mag = {e: abs(v) for e, v in corr.items()}
    ctrl_mag = {e: 1.0 for e in cand_edges}

    def score(es, net_tfs, net_tgts, cells=None):
        if cells is None:
            tfb = sorted(net_tfs & gold_tfs); tgb = sorted(net_tgts & gold_tgts)
            cells = [(t, g) for t in tfb for g in tgb]
        sc = [es.get(c, 0.0) for c in cells]; lb = [1 if c in gold_pos else 0 for c in cells]
        roc, pr = ranking_metrics(sc, lb)
        mw = roc_auc_score(lb, sc) if 0 < sum(lb) < len(lb) else float('nan')
        ap = average_precision_score(lb, sc) if sum(lb) > 0 else float('nan')
        return dict(ncells=len(cells), npos=sum(lb), auroc=roc, auroc_mw=mw, auprc=pr, ap=ap), cells

    rows = []
    print("\n== NATIVE rectangles ==")
    for nm, es, tf, tg in [("fit_lasso", fit_mag, set(t for t,_ in fit_nz), set(g for _,g in fit_nz)),
                           ("fit_corr",  corr_mag, set(t for t,_ in corr_mag), set(g for _,g in corr_mag)),
                           ("prior_signs", ctrl_mag, cand_tfs, cand_tgts)]:
        m, _ = score(es, tf, tg)
        rows.append(dict(scope="native", network=nm, **m))
        print(f"  {nm:12s} AUROC(trapz)={m['auroc']:.3f} AUROC(MW)={m['auroc_mw']:.3f} "
              f"AUPRC={m['auprc']:.3f} [{m['ncells']} cells, {m['npos']} pos]")

    print("\n== COMMON rectangle (paired) ==")
    ctf = (set(t for t,_ in fit_nz) & set(t for t,_ in corr_mag) & cand_tfs) & gold_tfs
    ctg = (set(g for _,g in fit_nz) & set(g for _,g in corr_mag) & cand_tgts) & gold_tgts
    cells = [(t, g) for t in sorted(ctf) for g in sorted(ctg)]
    for nm, es in [("fit_lasso", fit_nz), ("fit_corr", corr_mag), ("prior_signs", ctrl_mag)]:
        m, _ = score(es, None, None, cells)
        rows.append(dict(scope="common", network=nm, **m))
        print(f"  {nm:12s} AUROC(trapz)={m['auroc']:.3f} AUROC(MW)={m['auroc_mw']:.3f} "
              f"AUPRC={m['auprc']:.3f} [{m['ncells']} cells, {m['npos']} pos]")

    print("\n== CANDIDATE-ONLY (isolates the estimator) ==")
    cc = sorted(cand_edges)
    for nm, es in [("fit_lasso", fit_mag), ("fit_corr", corr_mag), ("prior_signs", ctrl_mag)]:
        m, _ = score(es, None, None, cc)
        rows.append(dict(scope="candidate_only", network=nm, **m))
        print(f"  {nm:12s} AUROC(MW)={m['auroc_mw']:.3f} AUPRC={m['auprc']:.3f} "
              f"[{m['ncells']} candidate cells, {m['npos']} gold-pos]   (control MW=0.5 by construction)")

    rec = [(coef[e], cand_sign[e]) for e in coef if e in gold_pos and cand_sign[e] != 0]
    conc = np.mean([np.sign(c) == s for c, s in rec])
    act = [(c, s) for c, s in rec if s == +1]; rep = [(c, s) for c, s in rec if s == -1]
    print(f"\n== SIGN concordance (fit coef vs RegulonDB sign, {len(rec)} recovered signed edges) ==")
    print(f"  overall {conc:.3f} | activators {np.mean([np.sign(c)==s for c,s in act]):.3f} (n={len(act)})"
          f" | repressors {np.mean([np.sign(c)==s for c,s in rep]):.3f} (n={len(rep)})")

    pd.DataFrame(rows).to_csv(f"{os.path.dirname(os.path.abspath(__file__))}/second_network_results.csv", index=False)

if __name__ == "__main__":
    if "--fetch" in sys.argv or not os.path.exists(f"{DATA_DIR}/net3_gold.tsv"):
        fetch()
    main()
