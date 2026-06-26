# Known limitations

What the deliverable is *not*, stated plainly so it is used within its evidence.

- **Not a dynamic predictor.** Out of sample the q_x-weighted multiplier does not
  beat a flat rate on any gene set. Use it as a mean-preserving, bounded
  steady-state descriptor and safety constraint, not as a forecast of mRNA.

- **TF mRNA is a proxy for TF activity.** Many cell-cycle masters (SBF/MBF,
  forkheads, the SWI5/ACE2 pair) are timed by localization or phosphorylation;
  genome-wide only ~44% of periodic proteins have periodic mRNA. Where protein and
  mRNA diverge, the couplings track the mRNA program, not the active TF.

- **Repressor couplings give a sign, not a strength.** Repressor edges are poorly
  identified (only ~2% have bootstrap intervals excluding zero, versus ~55% of
  activators). Treat a repressor edge as a categorical sign correction, not a
  magnitude.

- **Structural corroboration is partly curation-aligned.** The fitted network's
  gold-standard AUROC (0.815) is measured on a candidate set that shares lineage
  with the gold standard; against independent binding references it falls to
  0.53–0.66. This is corroboration within a curated candidate set, not de-novo
  discovery.

- **One time course, ~18 effective timepoints.** The default couplings are fit to a
  single, noisy, partially missing course, which caps the recoverable network and
  leaves many multi-regulator genes over-parameterized.

- **Joint-fit transfer is course-dependent.** The cross-dataset transfer gain is
  real and generalizes in direction, but its magnitude varies by course and it is a
  loss on the Orlando elutriation course; it is specific to the alpha-factor
  RNA-seq / microarray axis rather than universal.

- **Not a complete TF list.** Only TFs in the Teufel 2019 supplementary network are
  considered; known regulators absent from it (e.g. MBP1 -> CLN2) are not captured.

- **Bulk, synchronized courses only.** The transfer signal does not extend to a
  single-cell pseudotime reconstruction, where population-averaged ordering scrambles
  the timed regulator–target lag the couplings encode.
