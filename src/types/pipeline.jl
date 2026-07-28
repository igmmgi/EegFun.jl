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
    apply::Bool = false
    type::String = "hp"
    freq::Float64 = 0.1
    func::String = "filtfilt"
    method::String = "iir"
    order::Int = 2
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
    highpass::FilterSection = FilterSection(apply = true, type = "hp", freq = 0.1)
    lowpass::FilterSection = FilterSection(apply = true, type = "lp", freq = 30.0)
    ica_highpass::FilterSection = FilterSection(apply = true, type = "hp", freq = 1.0)
    ica_lowpass::FilterSection = FilterSection(apply = true, type = "lp", freq = 30.0)
end

"""
    EogConfig

Configuration for EOG (Electrooculogram) channel calculation and detection.

# Fields
- `vEOG_criterion::Float64`: Detection threshold for vertical EOG artifacts (in μV)
- `hEOG_criterion::Float64`: Detection threshold for horizontal EOG artifacts (in μV)
- `vEOG_channels::Vector{Vector{String}}`: Channel configuration for vertical EOG [channels1, channels2, output_channel]
- `hEOG_channels::Vector{Vector{String}}`: Channel configuration for horizontal EOG [channels1, channels2, output_channel]
"""
@kwdef struct EogConfig
    vEOG_criterion::Float64 = 50.0
    hEOG_criterion::Float64 = 30.0
    vEOG_channels::Vector{Vector{String}} = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]]
    hEOG_channels::Vector{Vector{String}} = [["F9"], ["F10"], ["hEOG"]]
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
    artifact_value_abs_criterion::Int = 100
    artifact_value_z_criterion::Float64 = 0.0
    extreme_value_abs_criterion::Int = 500
    artifact_interval_start::Union{Nothing,Float64} = nothing
    artifact_interval_end::Union{Nothing,Float64} = nothing
end

"""
    IcaConfig

Configuration for Independent Component Analysis.

# Fields
- `apply::Bool`: Whether to apply ICA
- `percentage_of_data::Float64`: Percentage of data to use for ICA (0-100)
- `component_method::Symbol`: Method for identifying artifact components (:correlation)
"""
@kwdef struct IcaConfig
    apply::Bool = true
    percentage_of_data::Float64 = 100.0
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
    apply::Bool = false
    line_frequencies::Vector{Float64} = [50.0]
    bandwidth::Float64 = 2.0
    sliding_win_length::Float64 = 4.0
    sliding_win_step::Float64 = 2.0
    time_bandwidth::Float64 = 3.0
    k_tapers::Int = 5
    p_value::Float64 = 0.05
    pad::Int = 2
end

"""
    ResampleConfig

Configuration for data resampling.

# Fields
- `apply::Bool`: Whether to apply resampling
- `target_rate::Int`: The target sampling rate in Hz
"""
@kwdef struct ResampleConfig
    apply::Bool = false
    target_rate::Int = 512
end

"""
    PreprocessConfig

Comprehensive configuration for EEG data preprocessing.
"""
@kwdef struct PreprocessConfig
    reference_channel::Symbol = :avg
    epoch_start::Float64 = -1.0
    epoch_end::Float64 = 1.0
    filter::FilterConfig = FilterConfig()
    cleanline::CleanLineConfig = CleanLineConfig()
    resample::ResampleConfig = ResampleConfig()
    eog::EogConfig = EogConfig()
    eeg::EegConfig = EegConfig()
    ica::IcaConfig = IcaConfig()
    neighbour_criterion::Float64 = 0.25
    channel_repair_method::Symbol = :spherical_spline
    interactive_continuous::Bool = false
    interactive_ica::Bool = false
    interactive_epochs::Bool = false
end

"""
    InputFilesConfig

Configuration for pipeline input files and directories.
"""
@kwdef struct InputFilesConfig
    directory::String = "."
    raw_data_files::Union{Vector{String},String} = "\\.bdf"
    recursive::Bool = false
    layout_file::String = "biosemi72.csv"
    epoch_condition_file::String = ""
end

"""
    OutputFilesConfig

Configuration for pipeline output files and saving options.
"""
@kwdef struct OutputFilesConfig
    directory::String = "./preprocessed_files"
    save_continuous_data_raw::Bool = true
    save_continuous_data_corrected::Bool = true
    save_ica_data::Bool = true
    save_epoch_data_raw::Bool = true
    save_epoch_data_corrected::Bool = true
    save_epoch_data::Bool = true
    save_erp_data_raw::Bool = true
    save_erp_data_corrected::Bool = true
    save_erp_data::Bool = true
end

"""
    FilesConfig

Top-level configuration for pipeline files (input and output).
"""
@kwdef struct FilesConfig
    input::InputFilesConfig = InputFilesConfig()
    output::OutputFilesConfig = OutputFilesConfig()
end

# === GENERIC DICTIONARY CONVERTER ===

"""
    from_dict(::Type{T}, dict::Dict) where {T}

Generic helper to recursively instantiate `@kwdef` struct `T` from a string-keyed dictionary.
Automatically matches dictionary keys to struct field names, converts scalar types (e.g. `String` -> `Symbol`, `Int` -> `Float64`),
coerces vectors, and recursively constructs nested structs.
"""
function from_dict(::Type{T}, dict::Dict) where {T}
    kwargs = Dict{Symbol,Any}()
    for field in fieldnames(T)
        field_type = fieldtype(T, field)
        key_str = string(field)

        val = if haskey(dict, key_str)
            dict[key_str]
        elseif field == :neighbour_criterion && haskey(dict, "layout") && haskey(dict["layout"], "neighbour_criterion")
            dict["layout"]["neighbour_criterion"]
        elseif field == :channel_repair_method && haskey(dict, "channel_repair") && haskey(dict["channel_repair"], "method")
            dict["channel_repair"]["method"]
        elseif field == :channel_repair_method && haskey(dict, "repair_method")
            dict["repair_method"]
        else
            nothing
        end

        isnothing(val) && continue

        kwargs[field] = _coerce_value(field_type, val)
    end
    return T(; kwargs...)
end

function _coerce_value(::Type{T}, val) where {T}
    if val isa T
        return val
    elseif val isa Dict && !isabstracttype(T) && isstructtype(T)
        return from_dict(T, val)
    elseif T == Symbol && val isa String
        return Symbol(val)
    elseif T == Int && val isa Number
        return Int(val)
    elseif T == Float64 && val isa Number
        return Float64(val)
    elseif T == Union{Vector{String},String}
        if val isa String
            return val
        elseif val isa Vector
            return String.(val)
        else
            return string(val)
        end
    elseif T <: Vector && val isa Vector
        inner_type = T.parameters[1]
        if inner_type == Vector{String}
            return Vector{Vector{String}}([String.(v isa Vector ? v : [v]) for v in val])
        else
            return inner_type.(val)
        end
    else
        return val
    end
end

"""
    PipelineConfig

Unified top-level configuration for the preprocessing pipeline.
"""
@kwdef struct PipelineConfig
    files::FilesConfig = FilesConfig()
    preprocess::PreprocessConfig = PreprocessConfig()
end

# === DICTIONARY CONSTRUCTORS ===

FilterSection(cfg::Dict) = from_dict(FilterSection, cfg)
FilterConfig(cfg::Dict) = from_dict(FilterConfig, cfg)
EogConfig(cfg::Dict) = from_dict(EogConfig, cfg)
EegConfig(cfg::Dict) = from_dict(EegConfig, cfg)
IcaConfig(cfg::Dict) = from_dict(IcaConfig, cfg)
CleanLineConfig(cfg::Dict) = from_dict(CleanLineConfig, cfg)
ResampleConfig(cfg::Dict) = from_dict(ResampleConfig, cfg)
PreprocessConfig(cfg::Dict) = from_dict(PreprocessConfig, cfg)
InputFilesConfig(cfg::Dict) = from_dict(InputFilesConfig, cfg)
OutputFilesConfig(cfg::Dict) = from_dict(OutputFilesConfig, cfg)
FilesConfig(cfg::Dict) = from_dict(FilesConfig, cfg)
PipelineConfig(cfg::Dict) = from_dict(PipelineConfig, cfg)
