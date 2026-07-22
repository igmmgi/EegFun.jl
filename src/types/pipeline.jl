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
- `artifact_value_abs_criterion::Int`: Threshold for artifact detection (μV)
- `extreme_value_abs_criterion::Int`: Threshold for extreme value detection (μV)
- `artifact_interval_start::Union{Nothing,Float64}`: Start time for artifact rejection interval (optional)
- `artifact_interval_end::Union{Nothing,Float64}`: End time for artifact rejection interval (optional)
"""
@kwdef struct EegConfig
    artifact_value_abs_criterion::Int
    artifact_value_z_criterion::Float64
    extreme_value_abs_criterion::Int
    artifact_interval_start::Union{Nothing,Float64} = nothing
    artifact_interval_end::Union{Nothing,Float64} = nothing
end

"""
    IcaConfig

Configuration for Independent Component Analysis.

# Fields
- `apply::Bool`: Whether to apply ICA
- `percentage_of_data::Float64`: Percentage of data to use for ICA (0-100)
- `component_method::Symbol`: Method for identifying artifact components.
  - `:correlation` (default): Correlation-based detection (EOG, ECG, line noise, spatial kurtosis)


"""
@kwdef struct IcaConfig
    apply::Bool
    percentage_of_data::Float64
    component_method::Symbol = :correlation
end

"""
    CleanLineConfig

Configuration for CleanLine line-noise removal.

# Fields
- `apply::Bool`: Whether to apply CleanLine
- `line_frequencies::Vector{Float64}`: Target line noise frequencies (e.g., [50.0])
- `bandwidth::Float64`: Bandwidth for line frequency detection (Hz)
- `sliding_win_length::Float64`: Window length in seconds
- `sliding_win_step::Float64`: Window step in seconds
- `time_bandwidth::Float64`: Time-bandwidth product for tapers
- `k_tapers::Int`: Number of tapers
- `p_value::Float64`: Significance threshold for F-test
- `pad::Int`: FFT padding factor
"""
@kwdef struct CleanLineConfig
    apply::Bool
    line_frequencies::Vector{Float64}
    bandwidth::Float64
    sliding_win_length::Float64
    sliding_win_step::Float64
    time_bandwidth::Float64
    k_tapers::Int
    p_value::Float64
    pad::Int
end

"""
    ResampleConfig

Configuration for data resampling.

# Fields
- `apply::Bool`: Whether to apply resampling
- `target_rate::Int`: The target sampling rate in Hz
"""
@kwdef struct ResampleConfig
    apply::Bool
    target_rate::Int
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
- `cleanline::CleanLineConfig`: CleanLine line noise removal configuration
- `resample::ResampleConfig`: Resampling configuration
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
    cleanline::CleanLineConfig
    resample::ResampleConfig
    eog::EogConfig
    eeg::EegConfig
    ica::IcaConfig
    neighbour_criterion::Float64
    interactive_continuous::Bool
    interactive_ica::Bool
    interactive_epochs::Bool
end

# === CONSTRUCTORS ===

"""Construct `EogConfig` from a configuration dictionary."""
function EogConfig(cfg::Dict)
    return EogConfig(
        vEOG_criterion = cfg["vEOG_criterion"],
        hEOG_criterion = cfg["hEOG_criterion"],
        vEOG_channels = cfg["vEOG_channels"],
        hEOG_channels = cfg["hEOG_channels"],
    )
end

"""Construct `FilterSection` from a configuration dictionary."""
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

"""Construct `FilterConfig` from a configuration dictionary."""
function FilterConfig(cfg::Dict)
    return FilterConfig(
        highpass = FilterSection(cfg["highpass"]),
        lowpass = FilterSection(cfg["lowpass"]),
        ica_highpass = FilterSection(cfg["ica_highpass"]),
        ica_lowpass = FilterSection(cfg["ica_lowpass"]),
    )
end

"""Construct `EegConfig` from a configuration dictionary."""
function EegConfig(cfg::Dict)
    return EegConfig(
        artifact_value_abs_criterion = Int(cfg["artifact_value_abs_criterion"]),
        artifact_value_z_criterion = Float64(cfg["artifact_value_z_criterion"]),
        extreme_value_abs_criterion = Int(cfg["extreme_value_abs_criterion"]),
        artifact_interval_start = haskey(cfg, "artifact_interval_start") && !isnothing(cfg["artifact_interval_start"]) ? Float64(cfg["artifact_interval_start"]) : nothing,
        artifact_interval_end = haskey(cfg, "artifact_interval_end") && !isnothing(cfg["artifact_interval_end"]) ? Float64(cfg["artifact_interval_end"]) : nothing,
    )
end

"""Construct `IcaConfig` from a configuration dictionary."""
function IcaConfig(cfg::Dict)
    return IcaConfig(
        apply = cfg["apply"],
        percentage_of_data = cfg["percentage_of_data"],
        component_method = Symbol(get(cfg, "component_method", "correlation")),
    )
end

"""Construct `CleanLineConfig` from a configuration dictionary."""
function CleanLineConfig(cfg::Dict)
    return CleanLineConfig(
        apply = cfg["apply"],
        line_frequencies = Vector{Float64}(cfg["line_frequencies"]),
        bandwidth = Float64(cfg["bandwidth"]),
        sliding_win_length = Float64(cfg["sliding_win_length"]),
        sliding_win_step = Float64(cfg["sliding_win_step"]),
        time_bandwidth = Float64(cfg["time_bandwidth"]),
        k_tapers = Int(cfg["k_tapers"]),
        p_value = Float64(cfg["p_value"]),
        pad = Int(cfg["pad"]),
    )
end

"""Construct `ResampleConfig` from a configuration dictionary."""
function ResampleConfig(cfg::Dict)
    return ResampleConfig(
        apply = cfg["apply"],
        target_rate = Int(cfg["target_rate"]),
    )
end

"""Construct `PreprocessConfig` from a configuration dictionary."""
function PreprocessConfig(cfg::Dict)
    return PreprocessConfig(
        reference_channel = Symbol(cfg["reference_channel"]),
        epoch_start = cfg["epoch_start"],
        epoch_end = cfg["epoch_end"],
        filter = FilterConfig(cfg["filter"]),
        cleanline = CleanLineConfig(cfg["cleanline"]),
        resample = ResampleConfig(cfg["resample"]),
        eog = EogConfig(cfg["eog"]),
        eeg = EegConfig(cfg["eeg"]),
        ica = IcaConfig(cfg["ica"]),
        neighbour_criterion = cfg["layout"]["neighbour_criterion"],
        interactive_continuous = get(cfg, "interactive_continuous", false),
        interactive_ica = get(cfg, "interactive_ica", false),
        interactive_epochs = get(cfg, "interactive_epochs", false),
    )
end
