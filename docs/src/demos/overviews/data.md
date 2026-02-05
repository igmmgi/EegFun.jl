## Overview

This demo demonstrates loading EEG data and understanding the basic data structures.

### Supported File Formats

- **BioSemi (.bdf)**: High-resolution EEG recordings
- **BrainVision (.vhdr)**: BrainProducts format
- **EEGLAB (.set)**: Matlab-based EEG format
- **FieldTrip (.mat)**: Matlab-based format

### Core Data Structures

**ContinuousData:**
- Raw EEG time series
- Contains electrode data, triggers, sampling rate, metadata

**EpochData:**
- Segmented trials around events
- Organized by experimental conditions

**ErpData:**
- Averaged event-related potentials
- One waveform per condition
