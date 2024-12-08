# Changelog
## [0.2.4] - 2024-12-08
### Added
- Add `-g` to control the minimum number of unique marker genes (default: 1) required for a species to report its genome copies. Increase `-g` (1 -> 2) lowers recall (detection limit: 0.125 -> 0.25) but improves precision.


## [0.2.3] - 2024-12-01
### Fixed
- Fix a bug introduced in v0.2.2 causing ties not resolved properly. Results should be identical to v0.2.1.


## [0.2.2] - 2024-11-30
### Changed
- Avoid direct float comparison ([315b795](https://github.com/xinehc/melon/commit/315b795ebd7f5654822e0e9a12855c375cc11774)).


## [0.2.1] - 2024-11-01
### Changed
- Revert filtering criteria back to v0.1.6. The old criteria turn out to be helpful in rescuing certain chimeric reads in very rare settings.
### Fixed
- Reduce peak memory usage.


## [0.2.0] - 2024-06-28
### Added
- [***breaking***] Add decoy protein sequences (RefSeq fungi, protozoa, viral, plant, and human GRCh38/hg38) which effectively trap non-prokaryotic reads and prevent them from inflating total prokaryotic genome copy estimates if the pre-filtering module (default with `Kraken2`) is not enabled. Pre-filtering is no longer necessary even if samples are contaminated with human DNA or other common eukaryotes/viruses, unless the mean genome size of prokaryotes needs to be estimated. See [8918168](https://github.com/xinehc/melon-supplementary/commit/891816897bb3c82dcfff7ff44b45907593ba0eac) for more details. This function requires a database released on or after 2024-06-28.
### Changed
- Simplify filtering criteria for alignments.


## [0.1.6] - 2024-05-30
### Changed
- Prevent `extract_sequence` from loading all marker-containing reads into memory.
- Change `-F` to `--frameshift` and `max_iteration` to `max_iterations` for consistency.
- Switch from figshare to zenodo for better database versioning.


## [0.1.5] - 2024-04-26
### Fixed
- Fix a bug causing `tqdm` being disabled ([3bbd087](https://github.com/xinehc/melon/commit/3bbd087b8867e3167973a746af14f1fd797f9746)).


## [0.1.4] - 2024-04-26
### Changed
- Use `tqdm` for logging.
- Reduce peak memory usage by parsing PAF files on the fly.


## [0.1.3] - 2024-03-29
### Added
- Output both gap-compressed and gap-uncompressed (BLAST-like) identity.
- Refine output format.
### Changed
- Change alignment filtering criteria: make `AS` cutoff more stringent, drop `MS`. See [7cc6dbd](https://github.com/xinehc/melon/commit/7cc6dbd866027cf5c1adaa5c69ed7919d8630607) for details.
### Fixed
- Fix a bug causing total genome copies not being properly calculated with diamond>=2.1.9.


## [0.1.2] - 2023-12-20
### Added
- Add gap-compressed ANI to output.


## [0.1.1] - 2023-11-29
### Added
- Add options to control EM early stop.
### Changed
- Use `scipy.sparse` to reduce peak memory usage and computational time.
- Change default terminal condition of EM (`max_iteration`: 100 -> 1000; `epsilon`: 1e-5 -> 1e-10).
### Fixed
- Fix a bug causing chimeric reads not being aggregated.


## [0.1.0] - 2023-10-08
### Added
- Output a json file to indicate the lineage of processed reads.
### Changed
- Make databases indexed by default.
- Use only kraken's prediction for removal of non-prokaryotic reads.
- Change `--db_kraken` to `--db-kraken` for consistency.
- Change `sp[0-9]+` to ` sp[0-9]+` for consistency.
- Change `copies` to `copy` in output files for consistency.
### Fixed
- Prevent numpy from using all logical cores.


## [0.0.1] - 2023-09-19
### Added
- First release.
