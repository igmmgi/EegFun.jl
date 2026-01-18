"""
    Preprocessing Configuration Types

This module contains all configuration types specific to EEG data preprocessing,
including filtering, referencing, artifact detection, and ICA settings.
"""

"""
    FilterSection

Configuration for a single filter (highpass, lowpass, etc.).

# Fields
- `apply::Bool`: Whether to apply this filter
- `type::String`: Filter type ("hp"=highpass, "lp"=lowpass)
- `freq::Float64`: Cutoff frequency in Hz
- `func::String`: Filter function ("filt" or "filtfilt")
- `method::String`: Filter method ("iir" or "fir")
- `order::Int`: Filter order
"""
@kwdef struct FilterSection
    apply::Bool
    type::String
    freq::Float64
    func::String
    method::String
    order::Int
end

"""
    FilterConfig

Configuration for all filters used in preprocessing.

# Fields
- `highpass::FilterSection`: Highpass filter settings
- `lowpass::FilterSection`: Lowpass filter settings
- `ica_highpass::FilterSection`: Highpass filter for ICA data
- `ica_lowpass::FilterSection`: Lowpass filter for ICA data
"""
@kwdef struct FilterConfig
    highpass::FilterSection
    lowpass::FilterSection
    ica_highpass::FilterSection
    ica_lowpass::FilterSection
end

"""
    EogConfig

Configuration for EOG (Electrooculogram) channel calculation and detection.

This type contains all the parameters needed to configure EOG channel calculation
and artifact detection, including channel selections and detection criteria.

# Fields
- `vEOG_criterion::Float64`: Detection threshold for vertical EOG artifacts (in μV)
- `hEOG_criterion::Float64`: Detection threshold for horizontal EOG artifacts (in μV)
- `vEOG_channels::Vector{Vector{String}}`: Channel configuration for vertical EOG [channels1, channels2, output_channel]
- `hEOG_channels::Vector{Vector{String}}`: Channel configuration for horizontal EOG [channels1, channels2, output_channel]
"""
@kwdef struct EogConfig
    vEOG_criterion::Float64
    hEOG_criterion::Float64
    vEOG_channels::Vector{Vector{String}}
    hEOG_channels::Vector{Vector{String}}
end

"""
    EegConfig

Configuration for EEG-specific preprocessing settings.

# Fields
- `artifact_value_criterion::Int`: Threshold for artifact detection (μV)
- `extreme_value_criterion::Int`: Threshold for extreme value detection (μV)
"""
@kwdef struct EegConfig
    artifact_value_abs_criterion::Int
    extreme_value_abs_criterion::Int
end

"""
    IcaConfig

Configuration for Independent Component Analysis.

# Fields
- `apply::Bool`: Whether to apply ICA
- `percentage_of_data::Float64`: Percentage of data to use for ICA (0-100)
"""
@kwdef struct IcaConfig
    apply::Bool
    percentage_of_data::Float64
end

"""
    PreprocessConfig

Comprehensive configuration for EEG data preprocessing.

This type contains all the parameters needed to configure the complete preprocessing
pipeline, including filtering, referencing, artifact detection, and ICA settings.

# Fields
- `reference_channel::Symbol`: Reference channel for rereferencing
- `epoch_start::Float64`: Start time for epoch extraction (seconds)
- `epoch_end::Float64`: End time for epoch extraction (seconds)
- `filter::FilterConfig`: Filter configuration
- `eog::EogConfig`: EOG channel calculation and detection settings
- `eeg::EegConfig`: EEG-specific preprocessing settings
- `ica::IcaConfig`: ICA configuration settings
- `neighbour_criterion::Float64`: Distance criterion (in mm) for channel neighbour definition
"""
@kwdef struct PreprocessConfig
    reference_channel::Symbol
    epoch_start::Float64
    epoch_end::Float64
    filter::FilterConfig
    eog::EogConfig
    eeg::EegConfig
    ica::IcaConfig
    neighbour_criterion::Float64
end

# === CONSTRUCTORS ===

# Constructor from dictionary
function EogConfig(cfg::Dict)
    return EogConfig(
        vEOG_criterion = cfg["vEOG_criterion"],
        hEOG_criterion = cfg["hEOG_criterion"],
        vEOG_channels = cfg["vEOG_channels"],
        hEOG_channels = cfg["hEOG_channels"],
    )
end

function FilterSection(cfg::Dict)
    return FilterSection(
        apply = cfg["apply"],
        type = cfg["type"],
        freq = cfg["freq"],
        func = cfg["func"],
        method = cfg["method"],
        order = cfg["order"],
    )
end

function FilterConfig(cfg::Dict)
    return FilterConfig(
        highpass = FilterSection(cfg["highpass"]),
        lowpass = FilterSection(cfg["lowpass"]),
        ica_highpass = FilterSection(cfg["ica_highpass"]),
        ica_lowpass = FilterSection(cfg["ica_lowpass"]),
    )
end

function EegConfig(cfg::Dict)
    return EegConfig(
        artifact_value_abs_criterion = Int(cfg["artifact_value_abs_criterion"]),
        extreme_value_abs_criterion = Int(cfg["extreme_value_abs_criterion"]),
    )
end

function IcaConfig(cfg::Dict)
    return IcaConfig(apply = cfg["apply"], percentage_of_data = cfg["percentage_of_data"])
end

function PreprocessConfig(cfg::Dict)
    return PreprocessConfig(
        reference_channel = Symbol(cfg["reference_channel"]),
        epoch_start = cfg["epoch_start"],
        epoch_end = cfg["epoch_end"],
        filter = FilterConfig(cfg["filter"]),
        eog = EogConfig(cfg["eog"]),
        eeg = EegConfig(cfg["eeg"]),
        ica = IcaConfig(cfg["ica"]),
        neighbour_criterion = cfg["layout"]["neighbour_criterion"],
    )
end
