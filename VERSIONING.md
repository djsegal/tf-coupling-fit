# Data release versioning

The handoff (the fitted couplings + means + q_x + trust score) is a versioned data
product. This note states how releases are versioned and how columns may change.

## Versioning scheme
- The data product carries a semantic version `MAJOR.MINOR.PATCH`, recorded in
  `output/tf_handoff/metadata.json` (`data_version`) and in `CHANGELOG.md`.
  - PATCH: re-run with no value changes (provenance/regeneration only).
  - MINOR: added columns/files, or value changes within the documented method that do
    not change a column's meaning (e.g. the mean-preservation fix).
  - MAJOR: a column is renamed/removed or its meaning changes, or the fit method changes
    in a way that alters released alpha values.
- Every release appends a `CHANGELOG.md` entry (date, version, what changed) and
  refreshes `CHECKSUMS.sha256`.

## Canonical format and (optional) serialization
- The canonical, universal form is CSV in `data/` and `output/tf_handoff/`. CSV is the
  source of truth; nothing else is required to use the product.
- A single typed/binary artifact (Arrow) or a one-file JSON bundle is deliberately NOT a
  package dependency: it would add a heavy dep for marginal convenience over CSV. A
  downstream consumer who wants one can trivially serialize the three CSVs with their own
  Arrow/JSON tooling; a minimal example lives in the workshop, not in the shipped package.

## Deprecation policy (open items from the data-dictionary audit)
The next MINOR release should resolve, with a deprecation note, the naming carried for
backward compatibility:
- `q_x` and `q_rank` are identical (confirmed for all rows); keep one, document the other
  as a deprecated alias for one release, then drop it.
- `target` vs `substrate` and `source` vs `sources` are used inconsistently across CSVs;
  standardize on `substrate` and `sources`, aliasing the old names for one release.
- Complete the column documentation for every shipped CSV in the package README data
  dictionary (several files currently document the file but not all columns).
