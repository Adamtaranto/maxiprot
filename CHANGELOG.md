# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Changed

- **Breaking:** `maxiprot filter` takes the input GFF3 as a positional argument
  (was `--gff`), matching `maxiprot extract`.
- **Breaking:** `--max-gap` renamed to `--locus-pad` (alias `--locus-gap`).
- **Breaking:** `--score-mode` options replaced with the documented scoring model:
  `ms_cov_pos` (default), `AS`, `ms`, `pid_cov`, `pid_cov_len`, `length`,
  `linear`, `geom`, with weight flags `--w-pid --w-cov --w-len --w-pos --w-ms
  --w-AS` and `--length-metric {frac,aa}` (removed: `as_ms`, `ms_only`,
  `weighted`, `length_biased`).
- **Breaking:** the TSV summary is no longer written to stderr by default; it is
  only written when `--out-tsv` is given (`-` = stdout).
- **Breaking:** gate thresholds (`--min-cov`, `--min-pid`) now always exclude
  failing candidates from selection; without `--strict`, a locus with no
  passing candidate emits its best-scoring candidate flagged `passes=False`.
- Minimum supported Python is now 3.10.
- `maxiprot extract` uses Biopython for translation and reverse complement:
  all NCBI translation tables are supported and IUPAC ambiguity codes are
  preserved.

### Fixed

- `pid_aa` is now computed from the `cs:Z` alignment string (fallback
  `nmatch/alen`) instead of `ms/alen`, which used a DP score as a match count
  and could exceed 1.0.
- `cov_aa` is now `(qend - qstart) / qlen` instead of the gap-inclusive
  `alen / qlen`.
- Emitted CDS features synthesize a single shared `ID` per transcript instead
  of reading a (usually absent) `ID` from miniprot CDS lines, which produced
  `ID=nan` on real miniprot input.
- Synthesized features are converted from 0-based PAF to 1-based GFF3
  coordinates (previously off by one).
- Gene/mRNA ID collisions under `--emit-mode all_passing` are deduplicated.
- Inputs lacking `ms:i`/`AS:i`/`fs:i`/`st:i` tags no longer crash.
- Output GFF3 is coordinate-sorted and includes `##sequence-region` directives.
- `maxiprot extract` no longer drops CDS parts when CDS lines precede their
  mRNA line, and gene-to-mRNA links no longer depend on feature order.
- `rich` and `numpy` are declared as runtime dependencies (a clean install of
  `maxiprot filter` previously failed to import).

## [0.6.0] - 2025-09-01

### Added

- Set up pre-commit hooks

### Changed

- Restructure as package with single `maxiprot` entrypoint.
- Sub-modules `maxiprot filter` and `maxiprot extract`
- env yaml now pip installs maxiprot from github repo

[Unreleased]: https://github.com/Adamtaranto/maxiprot/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/Adamtaranto/maxiprot/releases/tag/v0.6.0
