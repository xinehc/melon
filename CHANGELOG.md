# Changelog
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
### Changed
- Change alignment filtering criteria: make `AS` cutoff more stringent, drop `MS`. See [7cc6dbd](https://github.com/xinehc/melon/commit/7cc6dbd866027cf5c1adaa5c69ed7919d8630607) for details.

### Fixed
- Fix a bug causing total genome copies not being properly calculated with diamond>=2.1.9.

### Added
- Output both gap-compressed and gap-uncompressed (BLAST-like) identity.
- Refine output format.

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
