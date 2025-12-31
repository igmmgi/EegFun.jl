# === ABSTRACT TYPES ===
"""
    EegData

Abstract supertype for all EEG data structures.

This type represents the base interface for EEG data objects, providing a common
interface for different data formats including continuous data, event-related
potentials, and epoch-based data.
"""
abstract type EegData end

"""
    SingleDataFrameEeg

Abstract type for EEG data stored in a single DataFrame.

This type represents EEG data where all samples are stored in a single DataFrame,
typically used for continuous data or event-related potentials that have been
averaged into a single time series.
"""
abstract type SingleDataFrameEeg <: EegData end

"""
    MultiDataFrameEeg

Abstract type for EEG data stored across multiple DataFrames.

This type represents EEG data where samples are distributed across multiple
DataFrames, typically used for epoch-based data where each epoch is stored
separately.
"""
abstract type MultiDataFrameEeg <: EegData end

# === ANALYSIS TYPES ===
"""
    AnalysisInfo

Stores metadata about data preprocessing and analysis parameters.

This type contains information about how the EEG data has been processed,
including filtering parameters, reference information, and other analysis
settings that affect the interpretation of the data.

# Fields
- `reference::Symbol`: Reference type used (e.g., :avg, :mastoid, :none)
- `hp_filter::Float64`: High-pass filter cutoff in Hz (0.0 if none)
- `lp_filter::Float64`: Low-pass filter cutoff in Hz (0.0 if none)
"""
@kwdef mutable struct AnalysisInfo
    reference::Symbol = :none
    hp_filter::Float64 = 0.0
    lp_filter::Float64 = 0.0
end

"""
    AnalysisSettings

Stores the final analysis settings applied by the user in the databrowser.

# Fields
- `hp_filter::Float64`: High-pass filter frequency (0.0 if not applied)
- `lp_filter::Float64`: Low-pass filter frequency (0.0 if not applied)
- `reference::Symbol`: Reference type used (:avg, :mastoid, :none, or channel name)
- `repaired_channels::Vector{Symbol}`: List of channels that were repaired
- `repair_method::Symbol`: Repair method used for all repaired channels
- `selected_regions::Vector{Tuple{Float64,Float64}}`: Time regions selected by user
- `removed_ica_components::Vector{Int}`: ICA components that were removed
"""
struct AnalysisSettings
    hp_filter::Float64
    lp_filter::Float64
    reference::Symbol
    repaired_channels::Vector{Symbol}
    repair_method::Symbol  # Single repair method for all repaired channels
    selected_regions::Vector{Tuple{Float64,Float64}}
    removed_ica_components::Vector{Int}
end

AnalysisSettings() = AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])

# === LAYOUT TYPES ===
"""
    Neighbours

Stores spatial neighbor information for an electrode in layout-based operations.

This type contains information about neighboring electrodes for spatial
interpolation, artifact detection, and other layout-dependent operations.
The distances and weights are typically calculated based on electrode
positions in 2D or 3D space.

# Fields
- `channels::Vector{Symbol}`: List of neighboring channel labels
- `distances::Vector{Float64}`: Distances to each neighbor in mm
- `weights::Vector{Float64}`: Interpolation weights for each neighbor
"""
struct Neighbours
    channels::Vector{Symbol}
    distances::Vector{Float64}
    weights::Vector{Float64}
end

"""
    Layout

Stores electrode layout information and spatial relationships.

This type contains the complete electrode layout information including
electrode positions in various coordinate systems (polar, 2D Cartesian,
3D Cartesian) and neighbor relationships for spatial operations.

# Fields
- `data::DataFrame`: DataFrame containing layout information with metadata groups
- `neighbours::Union{Nothing, OrderedDict{Symbol, Neighbours}}`: Dictionary of neighbours for each electrode
- `criterion::Union{Nothing, Float64}`: Distance criterion for neighbour calculation in mm
"""
mutable struct Layout
    data::DataFrame
    neighbours::Union{Nothing,OrderedDict{Symbol,Neighbours}}
    criterion::Union{Nothing,Float64}
end

"""
    Base.copy(layout::Layout) -> Layout

Create a copy of Layout with copied DataFrame and shared immutable fields.
The data DataFrame is copied with `copycols=true` to ensure independence,
while neighbours and criterion are shared since they are immutable.
"""
function Base.copy(layout::Layout)::Layout
    return Layout(
        copy(layout.data, copycols = true),
        layout.neighbours,  # Share since OrderedDict and Neighbours are immutable
        layout.criterion,    # Share since Float64 is immutable
    )
end

# === EEG DATA TYPES ===
"""
    ContinuousData

Stores continuous EEG data with associated layout and analysis information.

This type represents continuous EEG recordings where all time points are
stored in a single DataFrame. The data typically includes time series
for each electrode channel along with metadata columns.

# Fields
- `file::String`: Source filename
- `data::DataFrame`: DataFrame containing continuous data (without file column)
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the data in Hz
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct ContinuousData <: SingleDataFrameEeg
    file::String
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
end


"""
    ErpData

Stores event-related potential data with associated layout and analysis information.

This type represents averaged event-related potentials where multiple epochs
have been averaged together into a single time series. The data includes
the averaged ERP waveform for each electrode channel.

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Name of the condition
- `data::DataFrame`: DataFrame containing averaged ERP data (without condition/condition_name/n_epochs columns)
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the data in Hz
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
- `n_epochs::Int64`: Number of epochs that were averaged together
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
    EpochData

Stores epoch-based EEG data with associated layout and analysis information.

This type represents EEG data organized into individual epochs, where each
epoch is stored as a separate DataFrame. This format is useful for
event-related potential analysis and other epoch-based processing.

# Fields
- `file::String`: Source filename (constant across all epochs)
- `condition::Int64`: Condition number (constant across all epochs)
- `condition_name::String`: Name of the condition (constant across all epochs)
- `data::Vector{DataFrame}`: Vector of DataFrames, one for each epoch (without condition/condition_name/file columns; epoch column remains for original numbering)
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the data in Hz
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
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

Stores information about baseline correction applied to time-frequency data.

# Fields
- `method::Symbol`: Baseline method (`:db`, `:percent`, `:relchange`)
- `window::Tuple{Float64,Float64}`: Baseline time window (start, stop) in seconds
"""
@kwdef struct BaselineInfo
    method::Symbol
    window::Tuple{Float64,Float64}
end

"""
    TimeFreqData

Stores time-frequency analysis results for a single condition/average.

This type represents time-frequency power data where all time-frequency points
are stored in a single DataFrame with columns for time, frequency, and each
electrode channel. Suitable for averaged TF representations.

# Fields
- `file::String`: Source filename
- `condition::Int64`: Condition number
- `condition_name::String`: Name of the condition
- `data::DataFrame`: DataFrame with columns: time, freq, [electrode channels...]
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the original data in Hz
- `method::Symbol`: Analysis method (`:wavelet`, `:superlet`, `:multitaper`, `:spectrum`)
- `baseline::Union{BaselineInfo,Nothing}`: Baseline correction information (if applied)
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct TimeFreqData <: SingleDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    method::Symbol
    baseline::Union{BaselineInfo,Nothing}
    analysis_info::AnalysisInfo
end

"""
    Base.copy(tf_data::TimeFreqData) -> TimeFreqData

Create a copy of TimeFreqData with copied DataFrame and layout.
The data DataFrame is copied with `copycols=true` to ensure independence.
"""
function Base.copy(tf_data::TimeFreqData)::TimeFreqData
    return TimeFreqData(
        tf_data.file,
        tf_data.condition,
        tf_data.condition_name,
        copy(tf_data.data, copycols=true),
        copy(tf_data.layout),
        tf_data.sample_rate,
        tf_data.method,
        tf_data.baseline,  # BaselineInfo is immutable, can share
        tf_data.analysis_info,  # AnalysisInfo is small, can share or copy if needed
    )
end

"""
    TimeFreqEpochData

Stores time-frequency analysis results with individual trials preserved.

This type represents time-frequency power data organized into individual trials,
where each trial is stored as a separate DataFrame. Each DataFrame contains
columns for time, frequency, and each electrode channel.

# Fields
- `file::String`: Source filename (constant across all trials)
- `condition::Int64`: Condition number (constant across all trials)
- `condition_name::String`: Name of the condition (constant across all trials)
- `data::Vector{DataFrame}`: Vector of DataFrames, one per trial (columns: time, freq, [electrodes...])
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the original data in Hz
- `method::Symbol`: Analysis method (`:wavelet`, `:superlet`, `:multitaper`, `:spectrum`)
- `baseline::Union{BaselineInfo,Nothing}`: Baseline correction information (if applied)
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct TimeFreqEpochData <: MultiDataFrameEeg
    file::String
    condition::Int64
    condition_name::String
    data::Vector{DataFrame}
    layout::Layout
    sample_rate::Int64
    method::Symbol
    baseline::Union{BaselineInfo,Nothing}
    analysis_info::AnalysisInfo
end

"""
    AbstractInterval

Abstract type for interval specifications used in baseline correction and other time-based operations.
"""
abstract type AbstractInterval end

"""
    IntervalTime

Defines a time interval with start and end times.

This type represents a time window or interval, typically used for
defining analysis windows, artifact rejection periods, or other
time-based operations.

# Fields
- `start::Float64`: Start time of the interval in seconds
- `stop::Float64`: End time of the interval in seconds
"""
@kwdef struct IntervalTime <: AbstractInterval
    start::Float64
    stop::Float64
end

# Convert tuples to IntervalTime for convenience
IntervalTime(t::Tuple{Real,Real}) = IntervalTime(start = Float64(t[1]), stop = Float64(t[2]))

"""
    IntervalIndex

Defines a time interval using sample indices.

This type represents a time window or interval using sample indices rather
than time values. This is useful for operations that work directly with
sample positions, such as data slicing, artifact detection, or epoch
extraction.

# Fields
- `start::Int`: Start sample index (1-based)
- `stop::Int`: End sample index (inclusive)
"""
@kwdef struct IntervalIndex <: AbstractInterval
    start::Int
    stop::Int
end

const BaselineInterval = Union{AbstractInterval,Nothing}

"""
    EpochCondition

Defines parameters for extracting epochs for a specific experimental condition.

This type specifies the criteria for identifying and extracting epochs
from continuous EEG data based on trigger sequences and timing constraints.
It supports complex trigger patterns and timing relationships between events.

# Fields
- `name::String`: Descriptive condition name for identification
- `trigger_sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}`: Trigger sequences to match (e.g., [[1, 2, 3]], [[1, :any, 3]], [[1:5], [10:15]])
- `reference_index::Int`: Which trigger position is t=0 (1-based, default: 1)
- `timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}}`: Which trigger pairs to apply min/max intervals to (optional, default: nothing)
- `min_interval::Union{Nothing,Float64}`: Minimum time between specified trigger pairs in seconds (optional, default: nothing)
- `max_interval::Union{Nothing,Float64}`: Maximum time between specified trigger pairs in seconds (optional, default: nothing)
- `after::Union{Nothing,Int}`: Only search for sequences after this trigger value (optional, default: nothing)
- `before::Union{Nothing,Int}`: Only search for sequences before this trigger value (optional, default: nothing)
"""
@kwdef struct EpochCondition
    name::String
    trigger_sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}
    reference_index::Int = 1
    timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing
    min_interval::Union{Nothing,Float64} = nothing
    max_interval::Union{Nothing,Float64} = nothing
    after::Union{Nothing,Int} = nothing
    before::Union{Nothing,Int} = nothing
end


# === ICA TYPES ===
"""
    IcaPrms

Parameters for Independent Component Analysis (ICA) decomposition.

This type contains all the parameters needed to configure ICA decomposition
of EEG data, including learning rates, convergence criteria, and algorithm
specific settings.

# Fields
- `l_rate::Float64`: Learning rate for the ICA algorithm
- `max_iter::Int`: Maximum number of iterations for convergence
- `w_change::Float64`: Weight change threshold for convergence
- `anneal_deg::Float64`: Annealing degree for temperature scheduling
- `anneal_step::Float64`: Annealing step size for temperature updates
- `blowup::Float64`: Blowup factor for numerical stability
- `blowup_fac::Float64`: Blowup factor for algorithm stability
- `max_weight::Float64`: Maximum allowed weight value
- `restart_factor::Float64`: Factor for restarting stuck algorithms
- `degconst::Float64`: Degree constant for spherical coordinates
- `default_stop::Float64`: Default stopping criterion threshold
"""
mutable struct IcaPrms
    l_rate::Float64
    max_iter::Int
    w_change::Float64
    anneal_deg::Float64
    anneal_step::Float64
    blowup::Float64
    blowup_fac::Float64
    max_weight::Float64
    restart_factor::Float64
    degconst::Float64
    default_stop::Float64
end

"""
    InfoIca

Stores ICA analysis results and decomposition information.

This type contains the complete results of ICA decomposition including
the unmixing and mixing matrices, component statistics, and metadata
about the decomposition process.

# Fields
- `unmixing::Matrix{Float64}`: Unmixing matrix (sources = unmixing × data)
- `mixing::Matrix{Float64}`: Mixing matrix (data = mixing × sources)
- `sphere::Matrix{Float64}`: Sphering matrix for data preprocessing
- `variance::Vector{Float64}`: Variance explained by each component
- `scale::Float64`: Scaling factor applied to the data
- `mean::Vector{Float64}`: Mean vector subtracted from the data
- `ica_label::Vector{Symbol}`: Component labels (e.g., [:IC1, :IC2, ...])
- `layout::Layout`: Layout information for the ICA components (contains channel labels)
- `removed_activations::OrderedDict{Int, Matrix{Float64}}`: Removed component activations by epoch
- `kurtosis_signs::Vector{Bool}`: Boolean vector indicating if each component is sub-Gaussian (true = sub-Gaussian, false = super-Gaussian). For regular Infomax, all components are super-Gaussian (all false).
"""
struct InfoIca
    unmixing::Matrix{Float64}
    mixing::Matrix{Float64}
    sphere::Matrix{Float64}
    variance::Vector{Float64}
    scale::Float64
    mean::Vector{Float64}
    ica_label::Vector{Symbol}
    removed_activations::OrderedDict{Int,Matrix{Float64}}
    layout::Layout  # Layout information for the ICA components
    is_sub_gaussian::Vector{Bool}  # true = sub-Gaussian, false = super-Gaussian
end

"""
    ErpMeasurementsResult

Stores ERP measurement results with metadata about the analysis.

This type contains the measurement results DataFrame along with metadata
about the analysis type and window used, allowing the plotting function to
automatically use the correct parameters.

# Fields
- `data::DataFrame`: DataFrame containing measurement results
- `analysis_type::String`: Type of measurement (e.g., "mean_amplitude", "max_peak_latency")
- `analysis_window::Function`: Analysis window predicate function used
- `baseline_window::Function`: Baseline window predicate function used (if any)
"""
struct ErpMeasurementsResult
    data::DataFrame
    analysis_type::String
    analysis_window::Function
    analysis_window_desc::String
    baseline_window::Union{Function,Nothing}
    baseline_window_desc::String
end

# Make it show as a DataFrame for convenience
function Base.show(io::IO, result::ErpMeasurementsResult)
    println(io, "ErpMeasurementsResult")
    println(io, "  Analysis type: $(getfield(result, :analysis_type))")
    println(io, "  Analysis window: $(getfield(result, :analysis_window_desc))")
    println(io, "  Baseline window: $(getfield(result, :baseline_window_desc))")
    println(io, "\nResults:")
    show(io, getfield(result, :data))
end

"""
    TriggerInfo

Stores trigger count information in a DataFrame format with pretty printing support.

This type wraps a DataFrame containing trigger counts and provides automatic
pretty table formatting when displayed.

# Fields
- `data::DataFrame`: The underlying DataFrame with trigger count information
"""
struct TriggerInfo
    data::DataFrame
end

# === CONFIGURATION TYPES ===
# Preprocessing configuration types are now in preprocess_config.jl

# === DISPLAY FUNCTIONS ===
function Base.show(io::IO, layout::Layout)
    n_electrodes = size(layout.data, 1)
    has_2d = has_2d_coords(layout)
    has_3d = has_3d_coords(layout)
    has_neigh = has_neighbours(layout)

    println(io, "Layout ($n_electrodes channels)")
    println(
        io,
        "2D coords: $(has_2d ? "✓" : "✗"), 3D coords: $(has_3d ? "✓" : "✗"), Neighbours: $(has_neigh ? "✓" : "✗")",
    )

    if has_neigh
        avg_neighbours = average_number_of_neighbours(layout.neighbours)
        println(io, "Criterion: $(layout.criterion), Avg neighbours: $(round(avg_neighbours, digits=1))")
    end
    println(io)

    # Format numeric columns for display
    display_data = copy(layout.data)
    for col in names(display_data)
        if eltype(display_data[!, col]) <: Number
            display_data[!, col] = [@sprintf("%.2f", val) for val in display_data[!, col]]
        end
    end

    # Don't print too much noise!
    if n_electrodes <= 10
        PrettyTables.pretty_table(io, display_data, alignment = :r)
    else
        # Show first 3, ellipsis, last 3
        first_rows = display_data[1:3, :]
        last_rows = display_data[(end-2):end, :]

        # Create ellipsis row more concisely
        ellipsis_row = DataFrame([col => ["..."] for col in names(display_data)])
        combined_data = vcat(first_rows, ellipsis_row, last_rows)

        PrettyTables.pretty_table(io, combined_data, alignment = :r)
        println(io, "\n[showing first 3 and last 3 of $n_electrodes electrodes]")
    end
end

# Compact display for arrays
function Base.show(io::IO, ::MIME"text/plain", layout::Layout)
    show(io, layout)
end

# Custom show method for neighbours OrderedDict
function Base.show(io::IO, neighbours_dict::OrderedDict{Symbol,Neighbours})
    show(io, MIME"text/plain"(), neighbours_dict)
end

filename(dat::BiosemiDataFormat.BiosemiData)::String = basename_without_ext(dat.filename)
filename(dat::ContinuousData)::String = dat.file
filename(dat::ErpData)::String = dat.file
filename(dat::EpochData)::String = dat.file
filename(dat::TimeFreqData)::String = dat.file
filename(dat::SingleDataFrameEeg)::String = dat.data.file[1]  # Fallback for other SingleDataFrameEeg types
filename(dat::MultiDataFrameEeg)::String = dat.data[1].file[1]  # Fallback for other MultiDataFrameEeg types



function Base.show(io::IO, ::MIME"text/plain", neighbours_dict::OrderedDict{Symbol,Neighbours})
    n_electrodes = length(neighbours_dict)
    avg_neighbours = average_number_of_neighbours(neighbours_dict)

    println(
        io,
        "Neighbours Dictionary ($n_electrodes electrodes): Average neighbours per electrode: $(round(avg_neighbours, digits=1))",
    )
    println(io)

    # Don't print too much noise!
    if n_electrodes <= 6
        for (electrode, neighbours) in neighbours_dict
            _format_electrode(io, electrode, neighbours)
        end
    else
        entries = collect(neighbours_dict)
        # Show first 3, ellipsis, last 3
        for entry in entries[1:3]
            _format_electrode(io, entry[1], entry[2])
        end
        println(io, "⋮")
        for entry in entries[(end-2):end]
            _format_electrode(io, entry[1], entry[2])
        end
        println(io, "[showing first 3 and last 3 of $n_electrodes electrodes]")
    end
end


function Base.show(io::IO, dat::MultiDataFrameEeg)
    println(io, "File: $(filename(dat))")
    println(io, "Type: $(typeof(dat))")
    println(io, "Condition $(condition_number(dat)): $(condition_name(dat))")
    println(
        io,
        "Size: $(n_epochs(dat)) (epoch) x $(nrow(meta_data(dat))) (rows) x $(length(channel_labels(dat))) (columns)",
    )
    println(io, "Labels: ", print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::ContinuousData)
    println(io, "File: $(filename(dat))")
    println(io, "Type: $(typeof(dat))")
    println(io, "Labels: ", print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end


function Base.show(io::IO, dat::SingleDataFrameEeg)
    println(io, "File: $(filename(dat))")
    println(io, "Type: $(typeof(dat))")
    println(io, "Condition $(condition_number(dat)): $(condition_name(dat))")
    println(io, "Labels: ", print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end


function Base.show(io::IO, dat::TimeFreqData)
    n_times = length(unique(dat.data.time))
    n_freqs = length(unique(dat.data.freq))
    freqs = unique(dat.data.freq)
    println(io, "File: $(filename(dat))")
    println(io, "Type: TimeFreqData ($(dat.method))")
    println(io, "Condition $(condition_number(dat)): $(condition_name(dat))")
    println(io, "Size: $(n_times) times × $(n_freqs) freqs × $(length(channel_labels(dat))) channels")
    println(io, "Freq range: $(minimum(freqs)) - $(maximum(freqs)) Hz")
    if dat.baseline !== nothing
        println(io, "Baseline: $(dat.baseline.method) ($(dat.baseline.window[1]) to $(dat.baseline.window[2]) s)")
    end
    println(io, "Labels: ", print_vector(channel_labels(dat)))
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::TimeFreqEpochData)
    n_trials = length(dat.data)
    n_times = n_trials > 0 ? length(unique(dat.data[1].time)) : 0
    n_freqs = n_trials > 0 ? length(unique(dat.data[1].freq)) : 0
    freqs = n_trials > 0 ? unique(dat.data[1].freq) : Float64[]
    println(io, "File: $(filename(dat))")
    println(io, "Type: TimeFreqEpochData ($(dat.method))")
    println(io, "Condition $(condition_number(dat)): $(condition_name(dat))")
    println(io, "Size: $(n_trials) trials × $(n_times) times × $(n_freqs) freqs × $(length(channel_labels(dat))) channels")
    !isempty(freqs) && println(io, "Freq range: $(minimum(freqs)) - $(maximum(freqs)) Hz")
    if dat.baseline !== nothing
        println(io, "Baseline: $(dat.baseline.method) ($(dat.baseline.window[1]) to $(dat.baseline.window[2]) s)")
    end
    println(io, "Labels: ", print_vector(channel_labels(dat)))
    println(io, "Sample Rate: ", sample_rate(dat))
end


function Base.show(io::IO, dat::AnalysisInfo)
    print(io, "Reference: ", reference(dat), ", ")
    print(io, "HP Filter: $(filter_info(dat)[1]), LP Filter: $(filter_info(dat)[2])")
end

########### ICA ##############
# Custom display method for InfoIca
function Base.show(io::IO, ica::InfoIca)
    n_components = length(ica.ica_label)
    n_channels = length(ica.layout.data.label)

    println(io, "InfoIca Result")
    println(io, "├─ Components: $n_components")
    println(io, "├─ Channels: $n_channels")
    println(io, "├─ Scale: $(round(ica.scale, digits=3))")
    println(io, "├─ Top 5 variance explained:")

    # Show top 5 components with their variance
    for i = 1:min(5, n_components)
        var_pct = round(ica.variance[i] * 100, digits = 2)
        println(io, "│  $(ica.ica_label[i]): $(var_pct)%")
    end

    if n_components > 5
        println(io, "│  ... and $(n_components - 5) more components")
    end

    println(io, "├─ Matrix sizes:")
    println(io, "│  ├─ Unmixing: $(size(ica.unmixing))")
    println(io, "│  ├─ Mixing: $(size(ica.mixing))")
    println(io, "│  └─ Sphere: $(size(ica.sphere))")

    println(
        io,
        "└─ Channel labels: $(join(ica.layout.data.label[1:min(5, n_channels)], ", "))$(n_channels > 5 ? " ..." : "")",
    )
end

# Compact display for arrays
function Base.show(io::IO, ::MIME"text/plain", ica::InfoIca)
    show(io, ica)
end

# Custom copy method for InfoIca
function Base.copy(ica::InfoIca)::InfoIca
    return InfoIca(
        copy(ica.unmixing),
        copy(ica.mixing),
        copy(ica.sphere),
        copy(ica.variance),
        ica.scale,
        copy(ica.mean),
        copy(ica.ica_label),
        copy(ica.removed_activations),
        ica.layout,  # Layout is shared, no need to copy
        copy(ica.is_sub_gaussian),
    )
end
