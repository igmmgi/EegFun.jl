# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to (well, attempts to :-)) [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2026-07-16

### Added

- zero-allocation virtual padding in time-frequency pipelines (`tf_morlet`, `tf_stft`, `tf_multitaper`) to reduce memory usage and improve speed.
- Added interactive "Topoplot 3D" plotting option to the Data Browser's selected region popup menu.
- Added keyboard interactivity (Up/Down arrow keys) to dynamically scale colormap limits in 3D topographic plots (`LScene`/`Mesh` support).

### Changed

- `resample` API has transitioned to using an explicit `target_rate` argument (e.g., target sample rate in Hz) instead of a downsampling `factor`.

### Deprecated

### Removed

- Removed old `factor`-based resampling logic and corresponding test validations.

### Fixed

- Fixed metadata preservation (`:time` and `condition` columns) edge cases during Time-Frequency data extraction.
- Fixed `plot_topography` contour scaling issues that previously resulted in white rendering artifacts for peak values outside of standard ranges.
