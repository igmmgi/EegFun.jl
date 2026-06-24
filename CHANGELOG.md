# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - Unreleased

### Added
- **Automatic Dataset Download**: Added `download_eegfun_datasets()` leveraging `DataDeps.jl` to streamline caching and accessing tutorial datasets.
- **Data Loaders**: Added comprehensive test coverage for multiple importers including EEGLAB (`.set`), BrainVision (`.vhdr`), BioSemi (`.bdf`), and FieldTrip (`.mat`).
- **Data Loaders**: Added `bids_import.jl` and `xdf.jl` for extended native dataset support.
- **Representational Similarity Analysis (RSA)**: Core methods and cross-validation workflows integrated under `src/analysis/rsa/` with public API declarations (`rsa`, `compare_models`, `add_noise_ceiling!`, etc.).
- **Plotting & GUI**: Added new GUI interactivity helpers, updated databrowser tools, and improved spatial layout visualizations (`plot_layout_2d`, `plot_topography`).
- **Documentation**: Substantial DIATAXIS restructure. Refactored `demos/` directory into robust auto-generated `tutorials/`. Reintroduced all essential Walkthroughs (N170, Posner), pedagogical Signal Processing examples, and ICA teaching demos.

### Changed
- **API Improvements**: Exposed `create_eegfun_data` as the unified constructor, fully replacing legacy dataframe conversions. 
- **API Visibility**: Declared multiple internal core types (`ContinuousData`, `EpochData`, `Layout`, `StatsResult`) and pipeline utilities as `public` to prevent namespace clashes while making them accessible.
- **Quality Assurance**: Integrated `Aqua.jl` into the GitHub Actions CI pipeline for robust code quality testing.
- **Configuration**: Overhauled the internal `config.jl` structure and template generation to simplify analysis settings.

### Deprecated
- `create_eeg_dataframe` has been removed and fully replaced by `create_eegfun_data`.
- Pipeline v2 internal workflows were consolidated back into a unified `pipeline.jl`.

### Fixed
- Fixed documentation generation bugs that resulted in empty stubs for `demos/` paths after branch reorganization.
- Resolved ICA layout coordinate crash failures through the `has_valid_coordinates()` Spatial Validation Pattern.

---

*(Historical changes prior to 0.4.0 are tracked in Git history)*
