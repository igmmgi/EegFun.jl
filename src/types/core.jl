"""
    EegFunData

Abstract supertype for all primary EegFun data and result structures.
"""
abstract type EegFunData end

"""
    EegData <: EegFunData

Abstract supertype for EEG data containers (continuous, epoched, ERP, spectrum, and time-frequency).
"""
abstract type EegData <: EegFunData end

"""
    SingleDataFrameEeg <: EegData

Abstract supertype for EEG data stored in a single `DataFrame` (e.g., `ContinuousData`, `ErpData`, `SpectrumData`, `TimeFreqData`).
"""
abstract type SingleDataFrameEeg <: EegData end

"""
    MultiDataFrameEeg <: EegData

Abstract supertype for trial-level EEG data stored across a `Vector{DataFrame}` (e.g., `EpochData`, `TimeFreqEpochData`).
"""
abstract type MultiDataFrameEeg <: EegData end

"""
    AnalysisInfo

Preprocessing and analysis metadata attached to EEG data containers.

# Fields
- `reference::Symbol = :none`: Reference scheme (`:avg`, `:mastoid`, `:none`, or channel label)
- `hp_filter::Float64 = 0.0`: High-pass cutoff in Hz (`0.0` if unfiltered)
- `lp_filter::Float64 = 0.0`: Low-pass cutoff in Hz (`0.0` if unfiltered)
- `sample_rate::Float64 = 0.0`: Sample rate in Hz
- `n_ica_components_removed::Int = 0`: Number of ICA components removed
- `n_channels_repaired::Int = 0`: Number of channels repaired
- `repaired_channels::Vector{Symbol} = Symbol[]`: Labels of repaired channels
"""
@kwdef mutable struct AnalysisInfo
    reference::Symbol = :none
    hp_filter::Float64 = 0.0
    lp_filter::Float64 = 0.0
    sample_rate::Float64 = 0.0
    n_ica_components_removed::Int = 0
    n_channels_repaired::Int = 0
    repaired_channels::Vector{Symbol} = Symbol[]
end

"""
    AnalysisSettings

Interactive preprocessing settings recorded from the databrowser.

# Fields
- `hp_filter::Float64 = 0.0`: High-pass filter cutoff in Hz
- `lp_filter::Float64 = 0.0`: Low-pass filter cutoff in Hz
- `reference::Symbol = :none`: Reference scheme
- `repaired_channels::Vector{Symbol} = Symbol[]`: Repaired channel labels
- `repair_method::Symbol = :none`: Channel repair method
- `selected_regions::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[]`: Selected time intervals in seconds
- `removed_ica_components::Vector{Int} = Int[]`: Removed ICA component indices
"""
@kwdef struct AnalysisSettings
    hp_filter::Float64 = 0.0
    lp_filter::Float64 = 0.0
    reference::Symbol = :none
    repaired_channels::Vector{Symbol} = Symbol[]
    repair_method::Symbol = :none
    selected_regions::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[]
    removed_ica_components::Vector{Int} = Int[]
end

"""
    Neighbors

Spatial neighbor channels, distances, and interpolation weights for a single electrode.

# Fields
- `channels::Vector{Symbol}`: Neighboring channel labels
- `distances::Vector{Float64}`: Distances to each neighbor in normalized coordinate units
- `weights::Vector{Float64}`: Spatial interpolation weights
"""
struct Neighbors
    channels::Vector{Symbol}
    distances::Vector{Float64}
    weights::Vector{Float64}
end

"""
    Layout

Electrode positions (polar, 2D, and/or 3D Cartesian) and spatial neighbor graph.

# Fields
- `data::DataFrame`: Electrode coordinate table (`:label`, `:inc`, `:azi`, and optional Cartesian columns)
- `neighbors::Union{Nothing,Dict{Symbol,Neighbors}}`: Neighbor graph keyed by electrode label
- `criterion::Union{Nothing,Float64}`: Distance criterion used to compute neighbors
- `criterion_type::Union{Nothing,Symbol}`: Coordinate space used for neighbors (`:xy` or `:xyz`)
"""
mutable struct Layout
    data::DataFrame
    neighbors::Union{Nothing,Dict{Symbol,Neighbors}}
    criterion::Union{Nothing,Float64}
    criterion_type::Union{Nothing,Symbol}
end

Layout() = Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing, nothing)

function Base.copy(layout::Layout)::Layout
    neighbors_copy = if isnothing(layout.neighbors)
        nothing
    else
        Dict{Symbol,Neighbors}(k => Neighbors(copy(v.channels), copy(v.distances), copy(v.weights)) for (k, v) in layout.neighbors)
    end
    return Layout(copy(layout.data; copycols = true), neighbors_copy, layout.criterion, layout.criterion_type)
end

# Backward-compatible JLD2 deserialization for legacy files saved with UK spelling (:neighbours / Neighbours)
_convert_legacy_jld2_neighbors(::Nothing) = nothing

function _convert_legacy_jld2_neighbors(d::JLD2.SerializedDict)
    out = Dict{Symbol,Neighbors}()
    for p in d.kvvec
        out[Symbol(p.first)] = p.second isa Neighbors ? p.second : convert(Neighbors, p.second)
    end
    return out
end

function _convert_legacy_jld2_neighbors(d)
    out = Dict{Symbol,Neighbors}()
    for (k, v) in d
        out[Symbol(k)] = v isa Neighbors ? v : convert(Neighbors, v)
    end
    return out
end

Base.convert(::Type{Neighbors}, x::Union{JLD2.ReconstructedStatic,JLD2.ReconstructedMutable}) =
    Neighbors(Symbol.(x.channels), Float64.(x.distances), Float64.(x.weights))

function Base.convert(::Type{Layout}, x::JLD2.ReconstructedMutable)
    nb_raw = hasproperty(x, :neighbors) ? x.neighbors : (hasproperty(x, :neighbours) ? x.neighbours : nothing)
    return Layout(x.data, _convert_legacy_jld2_neighbors(nb_raw), x.criterion, x.criterion_type)
end

JLD2.rconvert(::Type{Neighbors}, nt::NamedTuple) = Neighbors(Symbol.(nt.channels), Float64.(nt.distances), Float64.(nt.weights))

function JLD2.rconvert(::Type{Layout}, nt::NamedTuple)
    nb_raw = haskey(nt, :neighbors) ? nt.neighbors : (haskey(nt, :neighbours) ? nt.neighbours : nothing)
    return Layout(nt.data, _convert_legacy_jld2_neighbors(nb_raw), nt.criterion, nt.criterion_type)
end

const _LEGACY_JLD2_TYPEMAP = Dict{String,Any}("EegFun.Neighbours" => Neighbors, "EegFun.Layout" => JLD2.Upgrade(Layout))

"""
    ContinuousData <: SingleDataFrameEeg

Continuous EEG recording with layout and preprocessing metadata.

# Fields
- `file::String`: Source filename
- `data::DataFrame`: Continuous samples (`:time`, metadata, channel columns, and optional extra columns)
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Sample rate in Hz
- `analysis_info::AnalysisInfo`: Preprocessing metadata
"""
mutable struct ContinuousData <: SingleDataFrameEeg
    file::String
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
end

"""
    ErpData <: SingleDataFrameEeg

Condition-averaged event-related potential (ERP) data.

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Condition label
- `data::DataFrame`: Averaged ERP samples
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Sample rate in Hz
- `analysis_info::AnalysisInfo`: Preprocessing metadata
- `n_epochs::Int64`: Number of epochs averaged
"""
mutable struct ErpData <: SingleDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
    n_epochs::Int64
end

"""
    EpochData <: MultiDataFrameEeg

Single-trial epoched EEG data for one experimental condition.

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Condition label
- `data::Vector{DataFrame}`: Per-epoch `DataFrame`s
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Sample rate in Hz
- `analysis_info::AnalysisInfo`: Preprocessing metadata
"""
mutable struct EpochData <: MultiDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data::Vector{DataFrame}
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
end

"""
    BaselineInfo

Baseline correction parameters applied to time-frequency data (`method` and `(start, stop)` `window` in seconds).
"""
@kwdef struct BaselineInfo
    method::Symbol
    window::Tuple{Float64,Float64}
end

"""
    TimeFreqData <: SingleDataFrameEeg

Condition-averaged time-frequency power and phase data. Accessing `.data` aliases `.data_power`.

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Condition label
- `data_power::DataFrame`: Power values (`:time`, `:freq`, and channel columns)
- `data_phase::DataFrame`: Phase values in radians (`:time`, `:freq`, and channel columns)
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Original sample rate in Hz
- `method::Symbol`: Spectral estimation method
- `baseline::Union{BaselineInfo,Nothing}`: Applied baseline correction, if any
- `analysis_info::AnalysisInfo`: Preprocessing metadata
"""
mutable struct TimeFreqData <: SingleDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data_power::DataFrame
    data_phase::DataFrame
    layout::Layout
    sample_rate::Int64
    method::Symbol
    baseline::Union{BaselineInfo,Nothing}
    analysis_info::AnalysisInfo
end

Base.getproperty(tf_data::TimeFreqData, name::Symbol) = name === :data ? getfield(tf_data, :data_power) : getfield(tf_data, name)
Base.propertynames(tf::TimeFreqData, private::Bool = false) = (fieldnames(typeof(tf))..., :data)

function Base.copy(tf_data::TimeFreqData)::TimeFreqData
    return TimeFreqData(
        tf_data.file,
        tf_data.condition,
        tf_data.condition_name,
        copy(tf_data.data_power; copycols = true),
        copy(tf_data.data_phase; copycols = true),
        copy(tf_data.layout),
        tf_data.sample_rate,
        tf_data.method,
        tf_data.baseline,
        copy(tf_data.analysis_info),
    )
end

"""
    TimeFreqEpochData <: MultiDataFrameEeg

Single-trial time-frequency power and phase data. Accessing `.data` aliases `.data_power`.

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Condition label
- `data_power::Vector{DataFrame}`: Per-trial power `DataFrame`s
- `data_phase::Vector{DataFrame}`: Per-trial phase `DataFrame`s (radians)
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Original sample rate in Hz
- `method::Symbol`: Spectral estimation method
- `baseline::Union{BaselineInfo,Nothing}`: Applied baseline correction, if any
- `analysis_info::AnalysisInfo`: Preprocessing metadata
"""
mutable struct TimeFreqEpochData <: MultiDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data_power::Vector{DataFrame}
    data_phase::Vector{DataFrame}
    layout::Layout
    sample_rate::Int64
    method::Symbol
    baseline::Union{BaselineInfo,Nothing}
    analysis_info::AnalysisInfo
end

Base.getproperty(tf_data::TimeFreqEpochData, name::Symbol) = name === :data ? getfield(tf_data, :data_power) : getfield(tf_data, name)
Base.propertynames(tf::TimeFreqEpochData, private::Bool = false) = (fieldnames(typeof(tf))..., :data)

function Base.copy(tf_data::TimeFreqEpochData)::TimeFreqEpochData
    return TimeFreqEpochData(
        tf_data.file,
        tf_data.condition,
        tf_data.condition_name,
        [copy(df; copycols = true) for df in tf_data.data_power],
        [copy(df; copycols = true) for df in tf_data.data_phase],
        copy(tf_data.layout),
        tf_data.sample_rate,
        tf_data.method,
        tf_data.baseline,
        copy(tf_data.analysis_info),
    )
end

"""
    SpectrumData <: SingleDataFrameEeg

Frequency-domain power spectral density data (`:freq` and electrode channel columns).

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Condition label
- `data::DataFrame`: Power spectrum values (μV²/Hz)
- `layout::Layout`: Electrode layout
- `sample_rate::Int64`: Original sample rate in Hz
- `method::Symbol`: Spectral estimation method (`:welch`, `:multitaper`, etc.)
- `analysis_info::AnalysisInfo`: Preprocessing metadata
"""
mutable struct SpectrumData <: SingleDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    method::Symbol
    analysis_info::AnalysisInfo
end

function Base.copy(spectrum_data::SpectrumData)::SpectrumData
    return SpectrumData(
        spectrum_data.file,
        spectrum_data.condition,
        spectrum_data.condition_name,
        copy(spectrum_data.data; copycols = true),
        copy(spectrum_data.layout),
        spectrum_data.sample_rate,
        spectrum_data.method,
        copy(spectrum_data.analysis_info),
    )
end

abstract type AbstractSelection end
struct TimeSelection <: AbstractSelection
    start::Float64
    stop::Float64
end
struct AllSelection <: AbstractSelection end

"""
    Interval

Type alias for time/sample interval specifications: `TimeSelection` (via `times(start, stop)`),
`AllSelection`, `(start, stop)` tuples, numeric vectors, sample predicate `Function`s, or `nothing`.
"""
const Interval = Union{AbstractSelection,Tuple{Real,Real},AbstractVector{<:Real},Function,Nothing}

"""
    EpochCondition

Trigger-matching and timing rules defining an experimental condition for epoch extraction.

# Fields
- `name::String`: Condition name
- `trigger_sequences::Vector{Vector{Any}}`: Trigger sequence(s) to match (integers, `:any`, ranges, or sets)
- `reference_index::Int = 1`: 1-based position within the sequence aligned to `t = 0`
- `timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing`: Sequence index pairs constrained by `min_interval`/`max_interval`
- `min_interval::Union{Nothing,Float64} = nothing`: Minimum inter-trigger interval in seconds
- `max_interval::Union{Nothing,Float64} = nothing`: Maximum inter-trigger interval in seconds
- `mask_before_trigger::Union{Nothing,Int} = nothing`: Only match sequences before this trigger code
- `mask_after_trigger::Union{Nothing,Int} = nothing`: Only match sequences after this trigger code
- `mask_triggers::Union{Nothing,Vector{Int}} = nothing`: Trigger codes defining search boundaries
- `mask_between_triggers::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing`: Trigger code pairs bounding valid search regions
"""
@kwdef struct EpochCondition
    name::String
    trigger_sequences::Vector{Vector{Any}} = Vector{Any}[]
    reference_index::Int = 1
    timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing
    min_interval::Union{Nothing,Float64} = nothing
    max_interval::Union{Nothing,Float64} = nothing
    mask_before_trigger::Union{Nothing,Int} = nothing
    mask_after_trigger::Union{Nothing,Int} = nothing
    mask_triggers::Union{Nothing,Vector{Int}} = nothing
    mask_between_triggers::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing

    function EpochCondition(
        name::String,
        trigger_sequences,
        reference_index::Int = 1,
        timing_pairs = nothing,
        min_interval = nothing,
        max_interval = nothing,
        mask_before_trigger = nothing,
        mask_after_trigger = nothing,
        mask_triggers = nothing,
        mask_between_triggers = nothing,
    )
        norm_seqs = if trigger_sequences isa Number || trigger_sequences isa Symbol || trigger_sequences isa AbstractRange
            Vector{Any}[Any[trigger_sequences]]
        elseif isempty(trigger_sequences)
            Vector{Any}[]
        elseif all(x -> x isa AbstractVector && !(x isa AbstractRange), trigger_sequences)
            Vector{Any}[collect(Any, s) for s in trigger_sequences]
        else
            Vector{Any}[collect(Any, trigger_sequences)]
        end
        min_iv = isnothing(min_interval) ? nothing : Float64(min_interval)
        max_iv = isnothing(max_interval) ? nothing : Float64(max_interval)
        return new(
            name,
            norm_seqs,
            reference_index,
            timing_pairs,
            min_iv,
            max_iv,
            mask_before_trigger,
            mask_after_trigger,
            mask_triggers,
            mask_between_triggers,
        )
    end
end

"""
    IcaPrms

Algorithm hyperparameters for Infomax, Picard, and AMICA ICA decomposition.
"""
@kwdef mutable struct IcaPrms
    l_rate::Float64 = 0.001
    max_iter::Int = 512
    w_change::Float64 = 1e-6
    anneal_deg::Float64 = 60.0
    anneal_step::Float64 = 0.9
    blowup::Float64 = 1e15
    blowup_fac::Float64 = 0.8
    max_weight::Float64 = 1e8
    restart_factor::Float64 = 0.9
    degconst::Float64 = 180.0 / π
    default_stop::Float64 = 1e-6
    use_gpu::Bool = false
    picard_m::Int = 7
    picard_lambda_min::Float64 = 0.01
    picard_ls_tries::Int = 10
    amica_lrate::Float64 = 0.1
    amica_max_iter::Int = 2000
    amica_m::Int = 3
end

"""
    InfoIca <: EegFunData

ICA decomposition result storing unmixing/mixing/sphering matrices, component variances, and layout.

# Fields
- `filename::String`: Source filename
- `unmixing::Matrix{Float64}`: Unmixing matrix (`sources = unmixing * data`)
- `mixing::Matrix{Float64}`: Mixing matrix (`data = mixing * sources`)
- `sphere::Matrix{Float64}`: Sphering matrix
- `variance::Vector{Float64}`: Proportion of variance explained per component
- `scale::Float64`: Data scaling factor
- `mean::Vector{Float64}`: Per-channel mean vector
- `ica_label::Vector{Symbol}`: Component labels (`:IC1`, `:IC2`, ...)
- `removed_activations::Dict{Int,Matrix{Float64}}`: Removed component activations keyed by epoch
- `layout::Layout`: Electrode layout
- `is_sub_gaussian::Vector{Bool}`: Sub-Gaussian component flags
"""
struct InfoIca <: EegFunData
    filename::String
    unmixing::Matrix{Float64}
    mixing::Matrix{Float64}
    sphere::Matrix{Float64}
    variance::Vector{Float64}
    scale::Float64
    mean::Vector{Float64}
    ica_label::Vector{Symbol}
    removed_activations::Dict{Int,Matrix{Float64}}
    layout::Layout
    is_sub_gaussian::Vector{Bool}
end

"""
    ErpMeasurementsResult <: EegFunData

Tabular ERP amplitude/latency measurements bundled with their `analysis_type` and time intervals.
"""
struct ErpMeasurementsResult <: EegFunData
    data::DataFrame
    analysis_type::String
    analysis_interval::Union{Interval,Nothing}
    baseline_interval::Union{Interval,Nothing}
end

"""
    TriggerInfo

Wrapper around a trigger-count `DataFrame` for formatted display.
"""
struct TriggerInfo
    data::DataFrame
end

# === DISPLAY & FILENAME HELPERS ===

function Base.show(io::IO, result::ErpMeasurementsResult)
    println(io, "ErpMeasurementsResult")
    println(io, "  Analysis type: $(getfield(result, :analysis_type))")
    println(io, "  Analysis interval: $(getfield(result, :analysis_interval))")
    println(io, "  Baseline interval: $(getfield(result, :baseline_interval))")
    println(io, "\nResults:")
    show(io, getfield(result, :data))
end

function Base.show(io::IO, cond::EpochCondition)
    println(io, "EpochCondition: \"$(cond.name)\"")
    n_seq = length(cond.trigger_sequences)
    println(io, "├─ Trigger sequences: ($n_seq sequence(s))")
    for (i, seq) in enumerate(cond.trigger_sequences)
        prefix = i == n_seq ? "└─" : "├─"
        println(io, "│  $prefix $(join(seq, " → "))")
    end

    max_len = isempty(cond.trigger_sequences) ? 0 : maximum(length, cond.trigger_sequences)
    max_len > 1 && println(io, "├─ Reference index: $(cond.reference_index) (t=0)")

    has_timing = !isnothing(cond.timing_pairs) || !isnothing(cond.min_interval) || !isnothing(cond.max_interval)
    if has_timing
        println(io, "├─ Timing constraints:")
        if !isnothing(cond.timing_pairs)
            pairs_str = join(("($p1 → $p2)" for (p1, p2) in cond.timing_pairs), ", ")
            println(io, "│  ├─ Timing pairs: $pairs_str")
        end
        !isnothing(cond.min_interval) && println(io, "│  ├─ Min interval: $(cond.min_interval) s")
        !isnothing(cond.max_interval) && println(io, "│  └─ Max interval: $(cond.max_interval) s")
    end

    if !isnothing(cond.mask_after_trigger) || !isnothing(cond.mask_before_trigger)
        println(io, "$(has_timing ? "└─" : "├─") Search boundaries:")
        if !isnothing(cond.mask_after_trigger)
            println(io, "   $(isnothing(cond.mask_before_trigger) ? "└─" : "├─") After trigger: $(cond.mask_after_trigger)")
        end
        !isnothing(cond.mask_before_trigger) && println(io, "   └─ Before trigger: $(cond.mask_before_trigger)")
    end
end

function _format_layout_rows(df::DataFrame)
    out = DataFrame()
    for col in names(df)
        vals = df[!, col]
        out[!, col] = eltype(vals) <: Number ? [@sprintf("%.2f", v) for v in vals] : vals
    end
    return out
end

function Base.show(io::IO, layout::Layout)
    n_electrodes = nrow(layout.data)
    has_2d = has_2d_coords(layout)
    has_3d = has_3d_coords(layout)
    has_neigh = has_neighbors(layout)

    println(io, "Layout ($n_electrodes channels)")
    println(io, "2D coords: $(has_2d ? "✓" : "✗"), 3D coords: $(has_3d ? "✓" : "✗"), Neighbors: $(has_neigh ? "✓" : "✗")")
    if has_neigh
        avg_nb = round(average_number_of_neighbors(layout.neighbors); digits = 1)
        println(io, "Criterion: $(layout.criterion) ($(layout.criterion_type)), Avg neighbors: $avg_nb")
    end
    println(io)

    if n_electrodes <= 10
        PrettyTables.pretty_table(io, _format_layout_rows(layout.data); alignment = :r)
    else
        first_rows = _format_layout_rows(layout.data[1:3, :])
        last_rows = _format_layout_rows(layout.data[(end-2):end, :])
        ellipsis_row = DataFrame([col => ["..."] for col in names(first_rows)])
        PrettyTables.pretty_table(io, vcat(first_rows, ellipsis_row, last_rows); alignment = :r)
        println(io, "\n[showing first 3 and last 3 of $n_electrodes electrodes]")
    end
end

Base.show(io::IO, ::MIME"text/plain", layout::Layout) = show(io, layout)

function Base.show(io::IO, neighbors_dict::Dict{Symbol,Neighbors})
    isempty(neighbors_dict) && return print(io, "Neighbors Dictionary (empty)")
    avg_nb = round(average_number_of_neighbors(neighbors_dict); digits = 1)
    print(io, "Neighbors Dictionary ($(length(neighbors_dict)) electrodes, avg neighbors: $avg_nb)")
end

function Base.show(io::IO, ::MIME"text/plain", neighbors_dict::Dict{Symbol,Neighbors})
    n_electrodes = length(neighbors_dict)
    avg_nb = round(average_number_of_neighbors(neighbors_dict); digits = 1)
    println(io, "Neighbors Dictionary ($n_electrodes electrodes): Average neighbors per electrode: $avg_nb\n")

    if n_electrodes <= 6
        for (electrode, nb) in neighbors_dict
            _format_electrode(io, electrode, nb)
        end
    else
        entries = collect(neighbors_dict)
        for i = 1:3
            _format_electrode(io, entries[i][1], entries[i][2])
        end
        println(io, "⋮")
        for i = (n_electrodes-2):n_electrodes
            _format_electrode(io, entries[i][1], entries[i][2])
        end
        println(io, "[showing first 3 and last 3 of $n_electrodes electrodes]")
    end
end

filename(
    dat::Union{
        BiosemiDataFormat.BiosemiData,
        EuropeanDataFormat.EdfData,
        BrainVisionDataFormat.BrainVisionData,
        FunctionalImageFormat.FifData,
        ExtensibleDataFormat.XdfData,
    },
)::String = basename_without_ext(dat.filename)
filename(dat::EegData)::String = dat.file

function Base.show(io::IO, dat::MultiDataFrameEeg)
    n_rows = isempty(dat.data) ? 0 : nrow(dat.data[1])
    println(io, "File: $(filename(dat))")
    println(io, "Type: $(typeof(dat))")
    println(io, "Condition $(condition_number(dat)): $(condition_name(dat))")
    println(io, "Size: $(n_epochs(dat)) (epoch) x $(n_rows) (rows) x $(n_channels(dat)) (columns)")
    println(io, "Labels: ", _print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::ContinuousData)
    println(io, "File: $(filename(dat))")
    println(io, "Type: $(typeof(dat))")
    println(io, "Labels: ", _print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::SingleDataFrameEeg)
    println(io, "File: $(filename(dat))")
    println(io, "Type: $(typeof(dat))")
    println(io, "Condition $(condition_number(dat)): $(condition_name(dat))")
    println(io, "Labels: ", _print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::AnalysisInfo)
    hp, lp = filter_info(dat)
    print(io, "Reference: ", reference(dat), ", HP Filter: $hp, LP Filter: $lp")
end

function Base.show(io::IO, ica::InfoIca)
    n_components = length(ica.ica_label)
    n_ch = length(ica.layout.data.label)

    println(io, "InfoIca Result")
    println(io, "├─ Filename: $(ica.filename)")
    println(io, "├─ Components: $n_components")
    println(io, "├─ Channels: $n_ch")
    println(io, "├─ Scale: $(round(ica.scale; digits = 3))")
    println(io, "├─ Top 5 variance explained:")
    for i = 1:min(5, n_components)
        println(io, "│  $(ica.ica_label[i]): $(round(ica.variance[i] * 100; digits = 2))%")
    end
    n_components > 5 && println(io, "│  ... and $(n_components - 5) more components")
    println(io, "├─ Matrix sizes:")
    println(io, "│  ├─ Unmixing: $(size(ica.unmixing))")
    println(io, "│  ├─ Mixing: $(size(ica.mixing))")
    println(io, "│  └─ Sphere: $(size(ica.sphere))")
    println(io, "└─ Channel labels: $(join(@view(ica.layout.data.label[1:min(5, n_ch)]), ", "))$(n_ch > 5 ? " ..." : "")")
end

Base.show(io::IO, ::MIME"text/plain", ica::InfoIca) = show(io, ica)

function Base.copy(ica::InfoIca)::InfoIca
    return InfoIca(
        ica.filename,
        copy(ica.unmixing),
        copy(ica.mixing),
        copy(ica.sphere),
        copy(ica.variance),
        ica.scale,
        copy(ica.mean),
        copy(ica.ica_label),
        Dict{Int,Matrix{Float64}}(k => copy(v) for (k, v) in ica.removed_activations),
        copy(ica.layout),
        copy(ica.is_sub_gaussian),
    )
end
