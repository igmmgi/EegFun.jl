# Abstract types
abstract type EegData end
abstract type SingleDataFrameEeg <: EegData end
abstract type MultiDataFrameEeg <: EegData end

"""
    AnalysisInfo

Basic information about data preprocessing.

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
    Neighbours

Stores neighbor information for an electrode in layout-based operations.

# Fields
- `electrodes::Vector{Symbol}`: List of neighboring electrode labels
- `distances::Vector{Float64}`: Distances to each neighbor
- `weights::Vector{Float64}`: Interpolation weights for each neighbor
"""
struct Neighbours
    electrodes::Vector{Symbol}
    distances::Vector{Float64}
    weights::Vector{Float64}
end

mutable struct Layout
    data::DataFrame
    neighbours::Union{Nothing, OrderedDict{Symbol, Neighbours}}
    criterion::Union{Nothing, Float64}
end


# Concrete types with basic fields
mutable struct ContinuousData <: SingleDataFrameEeg
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
end
 
mutable struct ErpData <: SingleDataFrameEeg
    data::DataFrame
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
    n_epochs::Int64
end
 
mutable struct EpochData <: MultiDataFrameEeg
    data::Vector{DataFrame}
    layout::Layout
    sample_rate::Int64
    analysis_info::AnalysisInfo
end

mutable struct CoordXY
    x::Any
    y::Any
end

mutable struct Coord
    coord::Array{CoordXY}
end


# Helper methods for neighbour management
has_neighbours(layout::Layout) = !isnothing(layout.neighbours)
criterion(layout::Layout) = layout.criterion

# Helper function for average neighbours calculation
function _average_number_of_neighbours(neighbours_dict::OrderedDict{Symbol, Neighbours})
    if isempty(neighbours_dict)
        return 0.0
    end
    
    total_neighbours = sum(length(neighbours.electrodes) for neighbours in values(neighbours_dict))
    return total_neighbours / length(neighbours_dict)
end

# Get neighbours with automatic computation if needed (mutating)
function get_neighbours_xy!(layout::Layout, distance_criterion::Real)
    if !has_neighbours(layout) || layout.criterion != distance_criterion
        layout.neighbours = get_electrode_neighbours_xy(layout, distance_criterion)
        layout.criterion = distance_criterion
    end
    return nothing
end

function get_neighbours_xyz!(layout::Layout, distance_criterion::Real)
    if !has_neighbours(layout) || layout.criterion != distance_criterion
        layout.neighbours = get_electrode_neighbours_xyz(layout, distance_criterion)
        layout.criterion = distance_criterion
    end
    return nothing
end

# Getter functions that return the neighbours
neighbours(layout::Layout) = layout.neighbours

# Clear neighbours (useful when layout data changes)
function clear_neighbours!(layout::Layout)
    layout.neighbours = nothing
    layout.criterion = nothing
end

# Display method for Layout
function Base.show(io::IO, layout::Layout)
    # Show metadata first
    n_electrodes = size(layout.data, 1)
    has_2d = has_2d_coords(layout)
    has_3d = has_3d_coords(layout)
    has_neigh = has_neighbours(layout)
    
    println(io, "Layout ($n_electrodes electrodes)")
    println(io, "2D coords: $(has_2d ? "✓" : "✗"), 3D coords: $(has_3d ? "✓" : "✗"), Neighbours: $(has_neigh ? "✓" : "✗")")
    if has_neigh
        avg_neighbours = _average_number_of_neighbours(layout.neighbours)
        println(io, "Criterion: $(layout.criterion), Avg neighbours: $(round(avg_neighbours, digits=1))")
    end
    println(io)
    
    # Format the data for display
    display_data = copy(layout.data)

    # Format numeric columns with different precision based on column type
    for col in names(display_data)
        if eltype(display_data[!, col]) <: Number
            # Coordinate columns - always show exactly 2 decimal places
            display_data[!, col] = [@sprintf("%.2f", val) for val in display_data[!, col]]
        end
    end
    
    # Create a combined view with first and last rows
    if n_electrodes <= 10
        # Show all rows if 10 or fewer
        PrettyTables.pretty_table(io, display_data, 
            header=names(display_data),
            alignment=:r,  # Right align all columns
            crop=:none
        )
    else
        # Show first 5 and last 5 rows with ellipsis
        first_rows = display_data[1:5, :]
        last_rows = display_data[end-4:end, :]
        
        # Create ellipsis row
        ellipsis_row = DataFrame()
        for col in names(display_data)
            ellipsis_row[!, col] = ["..."]
        end
        
        # Combine the data
        combined_data = vcat(first_rows, ellipsis_row, last_rows)
        
        PrettyTables.pretty_table(io, combined_data, 
            header=names(display_data),
            alignment=:r,  # Right align all columns
            crop=:none
        )
        
        println(io, "\n[showing first 5 and last 5 of $n_electrodes electrodes]")
    end
end

# Compact display for arrays
function Base.show(io::IO, ::MIME"text/plain", layout::Layout)
    show(io, layout)
end

# Custom show method for neighbours OrderedDict
function Base.show(io::IO, neighbours_dict::OrderedDict{Symbol, Neighbours})
    # Use the text/plain MIME type to ensure our custom method is used
    show(io, MIME"text/plain"(), neighbours_dict)
end

# Helper function to format a single electrode entry
function _format_electrode(io, electrode, neighbours)
    n_neighbours = length(neighbours.electrodes)
    avg_distance = round(mean(neighbours.distances), digits=1)
    
    println(io, "$(rpad(string(electrode), 6)): $(n_neighbours) neighbours (avg dist: $(avg_distance)mm)")
    if n_neighbours > 0
        neighbour_details = []
        for j in 1:n_neighbours
            neighbour = string(neighbours.electrodes[j])
            distance = round(neighbours.distances[j], digits=1)
            weight = round(neighbours.weights[j], digits=3)
            push!(neighbour_details, "$(neighbour) ($(distance),$(weight))")
        end
        neighbour_list = join(neighbour_details, ", ")
        println(io, "    Neighbours: $neighbour_list")
    end
    println(io)
end

function Base.show(io::IO, ::MIME"text/plain", neighbours_dict::OrderedDict{Symbol, Neighbours})
    n_electrodes = length(neighbours_dict)
    avg_neighbours = average_number_of_neighbours(neighbours_dict)
    
    println(io, "Neighbours Dictionary ($n_electrodes electrodes): Average neighbours per electrode: $(round(avg_neighbours, digits=1))")
    println(io)
    
    # Show entries based on size
    if n_electrodes <= 6 
        # Show all entries
        for (electrode, neighbours) in neighbours_dict
            _format_electrode(io, electrode, neighbours)
        end
    else
        entries = collect(neighbours_dict)
        # First/last 3 entries
        for i in 1:3
            _format_electrode(io, entries[i][1], entries[i][2])
        end
        println(io, "⋮")
        for i in (n_electrodes-2):n_electrodes
            _format_electrode(io, entries[i][1], entries[i][2])
        end
        
        println(io, "[showing first 3 and last 3 of $n_electrodes electrodes]")
    end
end

struct IntervalIdx
    interval_start::Int
    interval_end::Int
end

struct IntervalTime
    interval_start::Float64
    interval_end::Float64
end

"""
    EpochCondition

Defines parameters for extracting epochs for a specific experimental condition.

# Fields
- `name::String`: Descriptive condition name
- `trigger_sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}`: Trigger sequences to match (e.g., [[1, 2, 3]], [[1, :any, 3]], [[1:5], [10:15]])
- `reference_index::Int`: Which trigger position is t=0 (1-based, default: 1)
- `timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}}`: Which trigger pairs to apply min/max intervals to (optional, default: nothing)
- `min_interval::Union{Nothing,Float64}`: Minimum time between specified trigger pairs (optional, default: nothing)
- `max_interval::Union{Nothing,Float64}`: Maximum time between specified trigger pairs (optional, default: nothing)
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




# Basic information functions right with the types
channels(dat::EegData) = dat.layout.data.label

# For SingleDataFrameEeg (ContinuousData, ErpData)
all_channels(dat::SingleDataFrameEeg) = propertynames(dat.data)
metadata_columns(dat::SingleDataFrameEeg) = all_channels(dat)[1:findfirst(col -> col in channels(dat), all_channels(dat)) - 1]
extra_channels(dat::SingleDataFrameEeg) = setdiff(propertynames(data(dat)), [metadata_columns(dat); channels(dat)])

# For MultiDataFrameEeg (EpochData)
all_channels(dat::MultiDataFrameEeg) = propertynames(dat.data[1])  # Use first epoch as reference
metadata_columns(dat::MultiDataFrameEeg) = all_channels(dat)[1:findfirst(col -> col in channels(dat), all_channels(dat)) - 1]
extra_channels(dat::MultiDataFrameEeg) = setdiff(propertynames(data(dat)), [metadata_columns(dat); channels(dat)])


times(dat::SingleDataFrameEeg) = dat.data.time
times(dat::MultiDataFrameEeg) = first(dat.data).time  # assume all epochs are the same
sample_rate(dat::ContinuousData) = dat.sample_rate
sample_rate(dat::EegData) = dat.sample_rate
sample_rate(dat::DataFrame) = Int(1 / mean(diff(dat.time)))
reference(dat::EegData) = dat.analysis_info.reference
reference(dat::AnalysisInfo) = dat.reference
filter_info(dat::AnalysisInfo) = [dat.hp_filter, dat.lp_filter]
data(dat::SingleDataFrameEeg) = dat.data # single data frame
data(dat::MultiDataFrameEeg) = to_data_frame(dat) # single data frame with all epochs

# Timepoints, channel, epoch, and duration 
samples(dat::SingleDataFrameEeg) = dat.data.sample
samples(dat::MultiDataFrameEeg, epoch::Int) = dat.data[epoch].sample
n_samples(dat::SingleDataFrameEeg) = nrow(dat.data)
n_samples(dat::MultiDataFrameEeg) = nrow(first(dat.data))
n_channels(dat::EegData) = length(channels(dat))
n_epochs(dat::SingleDataFrameEeg) = 1
n_epochs(dat::MultiDataFrameEeg) = length(dat.data)
duration(dat::EegData) = last(times(dat)) - first(times(dat))

n_average(dat::ErpData) = dat.n_epochs

# channel information
has_channels(dat::EegData, chans::Vector{Symbol}) = all(in(channels(dat)), chans)
common_channels(dat1::EegData, dat2::EegData) = intersect(channels(dat1), channels(dat2))


function Base.show(io::IO, dat::EegData)
    println(io, "Type: $(typeof(dat))")
    println(io, "Size: $(n_epochs(dat)) (epoch) x $(n_samples(dat)) (rows) x $(n_channels(dat)) (columns)")
    println(io, "Labels: ", print_vector(channels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::AnalysisInfo)
    print(io, "Reference: ", reference(dat), ", ")
    print(io, "HP Filter: $(filter_info(dat)[1]), LP Filter: $(filter_info(dat)[2])")
end

########### ICA ##############
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

struct InfoIca
    unmixing::Matrix{Float64}
    mixing::Matrix{Float64}
    sphere::Matrix{Float64}
    variance::Vector{Float64}
    scale::Float64
    mean::Vector{Float64}
    ica_label::Vector{Symbol}
    data_label::Vector{Symbol}
    removed_activations::OrderedDict{Int, Matrix{Float64}}
end

# Custom display method for InfoIca
function Base.show(io::IO, ica::InfoIca)
    n_components = length(ica.ica_label)
    n_channels = length(ica.data_label)
    
    println(io, "InfoIca Result")
    println(io, "├─ Components: $n_components")
    println(io, "├─ Channels: $n_channels")
    println(io, "├─ Scale: $(round(ica.scale, digits=3))")
    println(io, "├─ Top 5 variance explained:")
    
    # Show top 5 components with their variance
    for i in 1:min(5, n_components)
        var_pct = round(ica.variance[i] * 100, digits=2)
        println(io, "│  $(ica.ica_label[i]): $(var_pct)%")
    end
    
    if n_components > 5
        println(io, "│  ... and $(n_components - 5) more components")
    end
    
    println(io, "├─ Matrix sizes:")
    println(io, "│  ├─ Unmixing: $(size(ica.unmixing))")
    println(io, "│  ├─ Mixing: $(size(ica.mixing))")
    println(io, "│  └─ Sphere: $(size(ica.sphere))")
    
    println(io, "└─ Channel labels: $(join(ica.data_label[1:min(5, n_channels)], ", "))$(n_channels > 5 ? " ..." : "")")
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
        copy(ica.data_label),
        copy(ica.removed_activations)
    )
end
