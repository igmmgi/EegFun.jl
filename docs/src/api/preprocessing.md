# Preprocessing

This section documents functions for EEG data preprocessing including filtering, referencing, and artifact removal.

## Filtering

### `plot_filter_response(filter, sample_rate, filter_freq, transition_band)`
Plot the frequency response of a digital filter with ideal response overlay.

**Arguments:**
- `filter`: A digital filter object (FIR coefficients or DSP.jl filter)
- `sample_rate`: Sampling rate in Hz
- `filter_freq`: Cutoff frequency in Hz  
- `transition_band`: Width of transition band in Hz

**Returns:** Figure and Axis objects

## Rereferencing

### `rereference!(dat::DataFrame, channel_labels, reference_channels)`
Apply rereferencing to specified channels in a DataFrame.

**Arguments:**
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{Symbol}`: Names of channels to rereference
- `reference_channels`: Channels to use as reference

**Effects:** Modifies input data in-place by subtracting reference signal

### `calculate_reference(dat::DataFrame, reference_channels)`
Calculate reference signal from specified channels.

**Returns:** Vector containing the average of specified reference channels

## Automated Processing

### `preprocess_eeg_data(config::String)`
Preprocess EEG data according to the specified configuration file.

**Arguments:**
- `config::String`: Path to the configuration file in TOML format

### `make_output_filename(output_dir, input_file, suffix)`
Create an output filename from input file path with given suffix.

**Returns:** String - Full output filename path

## File Utilities

### `check_files_exist(conditions, filetype)`
Check if files exist for all given conditions with specified filetype.

**Returns:** Bool - true if all files exist, false otherwise

### `get_files(directory, files)`
Get file paths from directory matching specified patterns.

**Returns:** Vector of file paths 