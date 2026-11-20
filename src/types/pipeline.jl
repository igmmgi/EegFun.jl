# Typed configuration structs and dictionary coercion for preprocessing pipelines.

######################################
######### Types & Constants ##########
######################################

"""
    FilterSection

Settings for a single high-pass or low-pass filter stage (`apply`, `type`, `freq`, `func`, `method`, `order`).

# Fields
- `apply::Bool = false`: Whether to apply this filter stage
- `type::Symbol = :hp`: Filter passband type (`:hp` or `:lp`)
- `freq::Float64 = 0.1`: Cutoff frequency in Hz
- `func::Symbol = :filtfilt`: Filtering function (`:filtfilt` zero-phase or `:filt` causal)
- `method::Symbol = :iir`: Filter design family (`:iir` or `:fir`)
- `order::Int = 2`: Filter order
"""
@kwdef struct FilterSection
    apply::Bool = false
    type::Symbol = :hp
    freq::Float64 = 0.1
    func::Symbol = :filtfilt
    method::Symbol = :iir
    order::Int = 2
end

"""
    FilterConfig

High-pass and low-pass filter settings for main preprocessing and ICA preparation.

# Fields
- `highpass::FilterSection`: Main preprocessing high-pass filter settings
- `lowpass::FilterSection`: Main preprocessing low-pass filter settings
- `ica_highpass::FilterSection`: High-pass filter settings applied before ICA training
- `ica_lowpass::FilterSection`: Low-pass filter settings applied before ICA training
"""
@kwdef struct FilterConfig
    highpass::FilterSection = FilterSection(apply = true, type = :hp, freq = 0.1, order = 1)
    lowpass::FilterSection = FilterSection(apply = false, type = :lp, freq = 30.0, order = 3)
    ica_highpass::FilterSection = FilterSection(apply = true, type = :hp, freq = 1.0, order = 1)
    ica_lowpass::FilterSection = FilterSection(apply = false, type = :lp, freq = 30.0, order = 3)
end

"""
    EogConfig

Vertical and horizontal EOG bipolar channel definitions and artifact thresholds (μV).

# Fields
- `vEOG_criterion::Float64 = 50.0`: Vertical EOG step/onset threshold in μV
- `hEOG_criterion::Float64 = 30.0`: Horizontal EOG step/onset threshold in μV
- `vEOG_channels::Vector{Vector{String}}`: `[superior_channels, inferior_channels, [output_label]]` specification for bipolar vEOG
- `hEOG_channels::Vector{Vector{String}}`: `[left_channels, right_channels, [output_label]]` specification for bipolar hEOG
"""
@kwdef struct EogConfig
    vEOG_criterion::Float64 = 50.0
    hEOG_criterion::Float64 = 30.0
    vEOG_channels::Vector{Vector{String}} = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]]
    hEOG_channels::Vector{Vector{String}} = [["F9"], ["F10"], ["hEOG"]]
end

"""
    EegConfig

EEG rereferencing and amplitude/z-score artifact rejection thresholds (μV).

# Fields
- `rereference::Union{Symbol,Vector{Symbol}} = :avg`: Reference channel(s) or `:avg`
- `artifact_value_abs_criterion::Float64 = 100.0`: Epoch absolute amplitude rejection threshold in μV
- `artifact_value_z_criterion::Float64 = 0.0`: Epoch z-score rejection threshold (`0.0` disables z-score rejection)
- `extreme_value_abs_criterion::Float64 = 500.0`: Continuous extreme-amplitude threshold in μV
- `artifact_interval_start::Union{Nothing,Float64} = nothing`: Optional start time (s) for epoch artifact detection
- `artifact_interval_end::Union{Nothing,Float64} = nothing`: Optional end time (s) for epoch artifact detection
"""
@kwdef struct EegConfig
    rereference::Union{Symbol,Vector{Symbol}} = :avg
    artifact_value_abs_criterion::Float64 = 100.0
    artifact_value_z_criterion::Float64 = 0.0
    extreme_value_abs_criterion::Float64 = 500.0
    artifact_interval_start::Union{Nothing,Float64} = nothing
    artifact_interval_end::Union{Nothing,Float64} = nothing
end

"""
    IcaConfig

ICA decomposition and artifact component identification settings for preprocessing pipelines.

# Fields
- `apply::Bool = true`: Whether to run ICA decomposition and component removal
- `percentage_of_data::Float64 = 100.0`: Percentage of clean samples subsampled for ICA training
- `component_method::Symbol = :correlation`: Automated artifact IC selection method
"""
@kwdef struct IcaConfig
    apply::Bool = true
    percentage_of_data::Float64 = 100.0
    component_method::Symbol = :correlation
end

"""
    CleanLineConfig

Multi-taper CleanLine sinusoidal line-noise removal settings.

# Fields
- `apply::Bool = false`: Whether to apply CleanLine sinusoidal noise removal
- `line_frequencies::Vector{Float64} = [50.0]`: Line-noise frequencies in Hz to test and remove
- `bandwidth::Float64 = 2.0`: Spectral search bandwidth around each line frequency in Hz
- `sliding_win_length::Float64 = 4.0`: Sliding window length in seconds
- `sliding_win_step::Float64 = 2.0`: Sliding window step size in seconds
- `time_bandwidth::Float64 = 3.0`: DPSS time-half-bandwidth product (`NW`)
- `k_tapers::Int = 5`: Number of Slepian tapers
- `p_value::Float64 = 0.05`: F-test significance threshold for sinusoidal component removal
- `pad::Int = 2`: Zero-padding factor exponent for FFT resolution
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

Resampling settings (`apply` and `target_rate` in Hz).

# Fields
- `apply::Bool = false`: Whether to resample continuous data
- `target_rate::Int = 512`: Target sampling frequency in Hz
"""
@kwdef struct ResampleConfig
    apply::Bool = false
    target_rate::Int = 512
end

"""
    AutoRejectConfig

Adaptive cross-validated Autoreject artifact rejection and channel repair settings.

# Fields
- `apply::Bool = false`: Whether to run Autoreject on extracted epochs
- `mode::Symbol = :local`: Autoreject mode (`:global` threshold or `:local` per-channel repair + rejection)
- `n_folds::Int = 10`: Number of cross-validation folds
- `n_candidates::Int = 30`: Number of candidate peak-to-peak thresholds evaluated
- `consensus_perc::Float64 = 0.2`: Consensus fraction limit (`κ`) for bad channels per epoch
- `max_repair_perc::Float64 = 0.1`: Maximum fraction of channels repaired in a single epoch
- `repair_method::Symbol = :spherical_spline`: Spatial interpolation method used during repair
- `seed::Union{Int,Nothing} = 42`: Random seed for fold partitioning (`nothing` for non-deterministic)
"""
@kwdef struct AutoRejectConfig
    apply::Bool = false
    mode::Symbol = :local
    n_folds::Int = 10
    n_candidates::Int = 30
    consensus_perc::Float64 = 0.2
    max_repair_perc::Float64 = 0.1
    repair_method::Symbol = :spherical_spline
    seed::Union{Int,Nothing} = 42
end

"""
    PreprocessConfig

Complete EEG preprocessing configuration combining filtering, CleanLine, resampling, EOG, ICA, Autoreject, and repair settings.

# Fields
- `reference_channel::Union{Symbol,Vector{Symbol}} = :avg`: Target reference channel(s) or `:avg`
- `epoch_start::Float64 = -1.0`: Epoch start time relative to trigger in seconds
- `epoch_end::Float64 = 1.0`: Epoch end time relative to trigger in seconds
- `filter::FilterConfig`: High-pass and low-pass filter settings
- `cleanline::CleanLineConfig`: CleanLine sinusoidal line-noise removal settings
- `resample::ResampleConfig`: Resampling settings
- `eog::EogConfig`: Bipolar EOG derivation and threshold settings
- `eeg::EegConfig`: EEG amplitude/z-score rejection settings
- `ica::IcaConfig`: ICA decomposition settings
- `autoreject::AutoRejectConfig`: Autoreject settings
- `neighbor_criterion::Float64 = 0.35`: Spatial distance criterion for channel neighbor graph construction
- `channel_repair_method::Symbol = :spherical_spline`: Channel interpolation method (`:spherical_spline` or `:neighbor_interpolation`)
- `interactive_continuous::Bool = false`: Open interactive databrowser during continuous preprocessing
- `interactive_ica::Bool = false`: Open interactive ICA component review GUI
- `interactive_epochs::Bool = false`: Open interactive epoch rejection GUI
"""
@kwdef struct PreprocessConfig
    reference_channel::Union{Symbol,Vector{Symbol}} = :avg
    epoch_start::Float64 = -1.0
    epoch_end::Float64 = 1.0
    filter::FilterConfig = FilterConfig()
    cleanline::CleanLineConfig = CleanLineConfig()
    resample::ResampleConfig = ResampleConfig()
    eog::EogConfig = EogConfig()
    eeg::EegConfig = EegConfig()
    ica::IcaConfig = IcaConfig()
    autoreject::AutoRejectConfig = AutoRejectConfig()
    neighbor_criterion::Float64 = 0.35
    channel_repair_method::Symbol = :spherical_spline
    interactive_continuous::Bool = false
    interactive_ica::Bool = false
    interactive_epochs::Bool = false
end

"""
    InputFilesConfig

Pipeline input directory, file pattern, layout file, and epoch condition file settings.

# Fields
- `directory::String = "."`: Directory containing raw EEG recordings
- `raw_data_files::Union{Vector{String},String} = "\\.bdf"`: Regex pattern or explicit list of raw filenames
- `recursive::Bool = false`: Whether to search subdirectories recursively
- `layout_file::String = "biosemi72.csv"`: Electrode layout filename or path
- `epoch_condition_file::String = ""`: Path to TOML/Julia epoch condition definition file
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

Pipeline output directory and intermediate/final stage saving flags.

# Fields
- `directory::String = "./preprocessed_files"`: Output directory for saved `.jld2` artifacts and logs
- `save_continuous_data_raw::Bool = false`: Save raw imported continuous data
- `save_continuous_data_corrected::Bool = false`: Save artifact-corrected continuous data
- `save_ica_data::Bool = true`: Save `InfoIca` decomposition results
- `save_epoch_data_raw::Bool = false`: Save pre-rejection extracted epochs
- `save_epoch_data_corrected::Bool = false`: Save ICA/channel-repaired epochs prior to final rejection
- `save_epoch_data::Bool = true`: Save final clean `EpochData`
- `save_erp_data_raw::Bool = false`: Save condition-averaged ERPs before artifact rejection
- `save_erp_data_corrected::Bool = false`: Save condition-averaged ERPs after repair
- `save_erp_data::Bool = true`: Save final clean `ErpData`
"""
@kwdef struct OutputFilesConfig
    directory::String = "./preprocessed_files"
    save_continuous_data_raw::Bool = false
    save_continuous_data_corrected::Bool = false
    save_ica_data::Bool = true
    save_epoch_data_raw::Bool = false
    save_epoch_data_corrected::Bool = false
    save_epoch_data::Bool = true
    save_erp_data_raw::Bool = false
    save_erp_data_corrected::Bool = false
    save_erp_data::Bool = true
end

"""
    FilesConfig

Input and output file configuration for preprocessing pipelines.

# Fields
- `input::InputFilesConfig`: Input file discovery and layout configuration
- `output::OutputFilesConfig`: Output directory and stage-saving flags
"""
@kwdef struct FilesConfig
    input::InputFilesConfig = InputFilesConfig()
    output::OutputFilesConfig = OutputFilesConfig()
end

"""
    PipelineConfig

Top-level pipeline configuration bundling `files::FilesConfig` and `preprocess::PreprocessConfig`.

# Fields
- `files::FilesConfig`: Input/output file settings
- `preprocess::PreprocessConfig`: Signal processing and artifact handling settings
"""
@kwdef struct PipelineConfig
    files::FilesConfig = FilesConfig()
    preprocess::PreprocessConfig = PreprocessConfig()
end

######################################
########### User Functions ###########
######################################

FilterSection(cfg::AbstractDict) = _from_dict(FilterSection, cfg)

FilterConfig(cfg::AbstractDict) = _from_dict(FilterConfig, cfg)

EogConfig(cfg::AbstractDict) = _from_dict(EogConfig, cfg)

EegConfig(cfg::AbstractDict) = _from_dict(EegConfig, cfg)

IcaConfig(cfg::AbstractDict) = _from_dict(IcaConfig, cfg)

CleanLineConfig(cfg::AbstractDict) = _from_dict(CleanLineConfig, cfg)

ResampleConfig(cfg::AbstractDict) = _from_dict(ResampleConfig, cfg)

AutoRejectConfig(cfg::AbstractDict) = _from_dict(AutoRejectConfig, cfg)

PreprocessConfig(cfg::AbstractDict) = _from_dict(PreprocessConfig, cfg)

InputFilesConfig(cfg::AbstractDict) = _from_dict(InputFilesConfig, cfg)

OutputFilesConfig(cfg::AbstractDict) = _from_dict(OutputFilesConfig, cfg)

FilesConfig(cfg::AbstractDict) = _from_dict(FilesConfig, cfg)

PipelineConfig(cfg::AbstractDict) = _from_dict(PipelineConfig, cfg)

######################################
######### Internal Functions #########
######################################

@generated function _from_dict(::Type{T}, dict::AbstractDict) where {T}
    fnames = fieldnames(T)
    exprs = Expr[]
    for f in fnames
        ft = fieldtype(T, f)
        key = string(f)
        def_expr = :(def.$f)
        if f === :neighbor_criterion && ft === Float64
            push!(exprs, :(_resolve_neighbor_criterion(dict, $def_expr)))
        elseif f === :channel_repair_method && ft === Symbol
            push!(exprs, :(_resolve_repair_method(dict, $def_expr)))
        elseif Nothing <: ft
            push!(exprs, :(haskey(dict, $key) ? _coerce_value($ft, dict[$key]) : $def_expr))
        else
            push!(exprs, quote
                if haskey(dict, $key)
                    v = dict[$key]
                    isnothing(v) ? $def_expr : _coerce_value($ft, v)
                else
                    $def_expr
                end
            end)
        end
    end
    return quote
        def = T()
        return T($(exprs...))
    end
end

_coerce_value(::Type{Symbol}, val) = val isa Symbol ? val : Symbol(val)

_coerce_value(::Type{Int}, val) = val isa Int ? val : Int(val)

_coerce_value(::Type{Float64}, val) = val isa Float64 ? val : Float64(val)

_coerce_value(::Type{Union{Nothing,Int}}, val) = isnothing(val) ? nothing : Int(val)

_coerce_value(::Type{Union{Nothing,Float64}}, val) = isnothing(val) ? nothing : Float64(val)

_coerce_value(::Type{Union{Symbol,Vector{Symbol}}}, val) = val isa AbstractVector ? Symbol[Symbol(x) for x in val] : Symbol(val)

_coerce_value(::Type{Union{Vector{String},String}}, val) = val isa AbstractVector ? String[String(x) for x in val] : string(val)

_coerce_value(::Type{Vector{Vector{String}}}, val) =
    val isa Vector{Vector{String}} ? val :
    Vector{String}[v isa AbstractVector ? String[String(x) for x in v] : String[String(v)] for v in val]

_coerce_value(::Type{Vector{E}}, val) where {E} = val isa Vector{E} ? val : E[x isa E ? x : _coerce_value(E, x) for x in val]

function _coerce_value(::Type{T}, val) where {T}
    val isa T && return val
    (val isa AbstractDict && !isabstracttype(T) && isstructtype(T)) && return _from_dict(T, val)
    return val
end

function _resolve_neighbor_criterion(dict::AbstractDict, default::Float64)::Float64
    haskey(dict, "neighbor_criterion") && !isnothing(dict["neighbor_criterion"]) && return Float64(dict["neighbor_criterion"])
    lay = get(dict, "layout", nothing)
    if lay isa AbstractDict
        haskey(lay, "neighbor_criterion") && !isnothing(lay["neighbor_criterion"]) && return Float64(lay["neighbor_criterion"])
    end
    return default
end

function _resolve_repair_method(dict::AbstractDict, default::Symbol)::Symbol
    if haskey(dict, "channel_repair_method") && !isnothing(dict["channel_repair_method"])
        return _coerce_value(Symbol, dict["channel_repair_method"])
    end
    cr = get(dict, "channel_repair", nothing)
    if cr isa AbstractDict && haskey(cr, "method") && !isnothing(cr["method"])
        return _coerce_value(Symbol, cr["method"])
    end
    if haskey(dict, "repair_method") && !isnothing(dict["repair_method"])
        return _coerce_value(Symbol, dict["repair_method"])
    end
    return default
end
