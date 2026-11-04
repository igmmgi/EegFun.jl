"""
    FilterSection

Settings for a single high-pass or low-pass filter stage (`apply`, `type`, `freq`, `func`, `method`, `order`).
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
"""
@kwdef struct IcaConfig
    apply::Bool = true
    percentage_of_data::Float64 = 100.0
    component_method::Symbol = :correlation
end

"""
    CleanLineConfig

Multi-taper CleanLine sinusoidal line-noise removal settings.
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
"""
@kwdef struct ResampleConfig
    apply::Bool = false
    target_rate::Int = 512
end

"""
    AutoRejectConfig

Adaptive cross-validated Autoreject artifact rejection and channel repair settings.
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
"""
@kwdef struct FilesConfig
    input::InputFilesConfig = InputFilesConfig()
    output::OutputFilesConfig = OutputFilesConfig()
end

"""
    PipelineConfig

Top-level pipeline configuration bundling `files::FilesConfig` and `preprocess::PreprocessConfig`.
"""
@kwdef struct PipelineConfig
    files::FilesConfig = FilesConfig()
    preprocess::PreprocessConfig = PreprocessConfig()
end

# === DICTIONARY DESERIALIZATION ===

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
    (val isa AbstractDict && !isabstracttype(T) && isstructtype(T)) && return from_dict(T, val)
    return val
end

function _resolve_neighbor_criterion(dict::AbstractDict, default::Float64)::Float64
    haskey(dict, "neighbor_criterion") && !isnothing(dict["neighbor_criterion"]) && return Float64(dict["neighbor_criterion"])
    if haskey(dict, "neighbour_criterion")
        @warn "'neighbour_criterion' is deprecated; use 'neighbor_criterion' instead." _module = nothing _file = nothing _line = nothing
        val = dict["neighbour_criterion"]
        return isnothing(val) ? default : Float64(val)
    end
    lay = get(dict, "layout", nothing)
    if lay isa AbstractDict
        haskey(lay, "neighbor_criterion") && !isnothing(lay["neighbor_criterion"]) && return Float64(lay["neighbor_criterion"])
        if haskey(lay, "neighbour_criterion")
            @warn "'neighbour_criterion' is deprecated; use 'neighbor_criterion' instead." _module = nothing _file = nothing _line = nothing
            val = lay["neighbour_criterion"]
            return isnothing(val) ? default : Float64(val)
        end
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

@generated function from_dict(::Type{T}, dict::AbstractDict) where {T}
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

FilterSection(cfg::AbstractDict) = from_dict(FilterSection, cfg)
FilterConfig(cfg::AbstractDict) = from_dict(FilterConfig, cfg)
EogConfig(cfg::AbstractDict) = from_dict(EogConfig, cfg)
EegConfig(cfg::AbstractDict) = from_dict(EegConfig, cfg)
IcaConfig(cfg::AbstractDict) = from_dict(IcaConfig, cfg)
CleanLineConfig(cfg::AbstractDict) = from_dict(CleanLineConfig, cfg)
ResampleConfig(cfg::AbstractDict) = from_dict(ResampleConfig, cfg)
AutoRejectConfig(cfg::AbstractDict) = from_dict(AutoRejectConfig, cfg)
PreprocessConfig(cfg::AbstractDict) = from_dict(PreprocessConfig, cfg)
InputFilesConfig(cfg::AbstractDict) = from_dict(InputFilesConfig, cfg)
OutputFilesConfig(cfg::AbstractDict) = from_dict(OutputFilesConfig, cfg)
FilesConfig(cfg::AbstractDict) = from_dict(FilesConfig, cfg)
PipelineConfig(cfg::AbstractDict) = from_dict(PipelineConfig, cfg)
