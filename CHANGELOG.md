# Changelog

## [0.1.0] - 2023-10-08
### Added
- Output a json file to indicate the lineage of processed reads.

### Changed
- Make databases indexed by default.
- Use only kraken's prediction for removal of non-prokaryotic reads.
- Change `--db_kraken` to `--db-kraken` for consistency.
- Change `sp[0-9]+` to ` sp[0-9]+` for consistenty.
- Change `copies` to `copy` in output files for consistency.

### Fixed
- Prevent numpy from using all logical cores.

## [0.0.1] - 2023-09-19
### Added
- First release.
