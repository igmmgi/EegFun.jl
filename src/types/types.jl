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

# === LAYOUT TYPES ===
"""
    Neighbours

Stores spatial neighbor information for an electrode in layout-based operations.

This type contains information about neighboring electrodes for spatial
interpolation, artifact detection, and other layout-dependent operations.
The distances and weights are typically calculated based on electrode
positions in 2D or 3D space.

# Fields
- `electrodes::Vector{Symbol}`: List of neighboring electrode labels
- `distances::Vector{Float64}`: Distances to each neighbor in mm
- `weights::Vector{Float64}`: Interpolation weights for each neighbor
"""
struct Neighbours
    electrodes::Vector{Symbol}
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
- `data::DataFrame`: DataFrame containing continuous data with metadata groups
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the data in Hz
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct ContinuousData <: SingleDataFrameEeg
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
- `data::DataFrame`: DataFrame containing averaged ERP data with metadata groups
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the data in Hz
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
- `n_epochs::Int64`: Number of epochs that were averaged together
"""
mutable struct ErpData <: SingleDataFrameEeg
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
- `data::Vector{DataFrame}`: Vector of DataFrames, one for each epoch
- `layout::Layout`: Layout object containing electrode positioning information
- `sample_rate::Int64`: Sample rate of the data in Hz
- `analysis_info::AnalysisInfo`: Analysis information and preprocessing metadata
"""
mutable struct EpochData <: MultiDataFrameEeg
    data::Vector{DataFrame}
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
end

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
@kwdef struct IntervalTime
    start::Float64
    stop::Float64
end


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
@kwdef struct IntervalIndex
    start::Int
    stop::Int
end


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
end

# === CONFIGURATION TYPES ===
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

# Constructor from dictionary
function EogConfig(cfg::Dict)
    return EogConfig(
        vEOG_criterion = cfg["vEOG_criterion"],
        hEOG_criterion = cfg["hEOG_criterion"],
        vEOG_channels = cfg["vEOG_channels"],
        hEOG_channels = cfg["hEOG_channels"]
    )
end

"""
    PreprocessConfig

Comprehensive configuration for EEG data preprocessing.

This type contains all the parameters needed to configure the complete preprocessing
pipeline, including filtering, referencing, artifact detection, and ICA settings.

# Fields
- `reference_channel::String`: Reference channel for rereferencing
- `filter::Dict`: Filter configuration dictionary
- `eog::EogConfig`: EOG channel calculation and detection settings
- `eeg::Dict`: EEG-specific preprocessing settings (artifact detection, etc.)
- `ica::Dict`: ICA configuration settings
"""
@kwdef struct PreprocessConfig
    reference_channel::String
    filter::Dict
    eog::EogConfig
    eeg::Dict
    ica::Dict
end

# Constructor from dictionary
function PreprocessConfig(cfg::Dict)
    return PreprocessConfig(
        reference_channel = cfg["reference_channel"],
        filter = cfg["filter"],
        eog = EogConfig(cfg["eog"]),
        eeg = cfg["eeg"],
        ica = cfg["ica"]
    )
end

# === DISPLAY FUNCTIONS ===
function Base.show(io::IO, layout::Layout)
    n_electrodes = size(layout.data, 1)
    has_2d = has_2d_coords(layout)
    has_3d = has_3d_coords(layout)
    has_neigh = has_neighbours(layout)

    println(io, "Layout ($n_electrodes channels)")
    println(io, "2D coords: $(has_2d ? "✓" : "✗"), 3D coords: $(has_3d ? "✓" : "✗"), Neighbours: $(has_neigh ? "✓" : "✗")")
    
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
filename(dat::SingleDataFrameEeg)::String = dat.data.file[1]
filename(dat::MultiDataFrameEeg)::String = dat.data[1].file[1]



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


function Base.show(io::IO, dat::EegData)
    println(io, "Type: $(typeof(dat))")
    println(
        io,
        "Size: $(n_epochs(dat)) (epoch) x $(nrow(meta_data(dat))) (rows) x $(length(channel_labels(dat))) (columns)",
    )
    println(io, "Labels: ", print_vector(channel_labels(dat)))
    println(io, "Duration: ", duration(dat), " S")
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
    )
end
