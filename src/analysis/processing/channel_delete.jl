"""
    channel_delete!(dat::EegData, channels_to_delete::Union{Symbol, Vector{Symbol}})

Delete channels from EEG data and layout.
Modifies the data in place by removing both the data columns and the layout rows.

# Arguments
- `dat::EegData`: The EEG data object to modify
- `channels_to_delete::Union{Symbol, Vector{Symbol}}`: Channel(s) to delete

# Returns
- `nothing` (modifies the data in place)

# Examples
```julia
channel_delete!(dat, :Fp1)                 # Delete a single channel
channel_delete!(dat, [:Fp1, :Fp2, :Diode]) # Delete multiple channels
```

# Notes
- Only channels that exist in the data will be deleted
- Updates both the data columns and the layout rows
- Clears any cached neighbour information in the layout since channels have changed
"""
function channel_delete!(dat::EegData, channels_to_delete::Vector{Symbol})

    # Check which channels exist in the data columns (includes extra channels)
    existing_data_channels = Set(all_labels(dat))
    channels_in_data = intersect(existing_data_channels, Set(channels_to_delete))

    # Check which channels exist in the layout (only for layout removal)
    existing_layout_channels = Set(dat.layout.data.label)
    channels_in_layout = intersect(existing_layout_channels, Set(channels_to_delete))

    if isempty(channels_in_data)
        @info "channel_delete!: No channels found to delete in data"
        return nothing
    end

    # Remove channels from layout (only if they exist in layout)
    if !isempty(channels_in_layout)
        filter!(:label => ∉(channels_in_layout), dat.layout.data)

        # Clear neighbours if they exist (since channels have changed)
        if has_neighbours(dat.layout)
            @info "channel_delete!: Clearing neighbours since channels have changed"
            clear_neighbours!(dat.layout)
        end
    end

    # Remove channels from data columns using multiple dispatch
    _delete_data_columns!(dat, channels_in_data)

    @info "channel_delete!: Deleted $(length(channels_in_data)) channel(s): $(join(string.(channels_in_data), ", "))"
    return nothing
end
channel_delete!(dat::EegData, channel::Symbol) = channel_delete!(dat, [channel])

function channel_delete!(dat::Vector{<:EegData}, channels_to_delete::Union{Symbol,Vector{Symbol}})
    for d in dat
        channel_delete!(d, channels_to_delete)
    end
    return dat
end

# Multiple dispatch for different EEG data types
function _delete_data_columns!(df::DataFrame, channels_to_delete::Set{Symbol})
    cols_to_remove = [c for c in propertynames(df) if c in channels_to_delete]
    if !isempty(cols_to_remove)
        select!(df, Not(cols_to_remove))
    end
    return nothing
end

function _delete_data_columns!(dat::SingleDataFrameEeg, channels_to_delete::Set{Symbol})
    _delete_data_columns!(dat.data, channels_to_delete)
end

function _delete_data_columns!(dat::MultiDataFrameEeg, channels_to_delete::Set{Symbol})
    isempty(dat.data) && return nothing
    cols_to_remove = [c for c in propertynames(dat.data[1]) if c in channels_to_delete]
    if !isempty(cols_to_remove)
        for df in dat.data
            select!(df, Not(cols_to_remove))
        end
    end
    return nothing
end

"""
    channel_delete(dat::EegData, channels_to_delete::Union{Symbol, Vector{Symbol}})

Create a copy of EEG data with specified channels deleted.

# Arguments
- `dat::EegData`: The EEG data object to copy and modify
- `channels_to_delete::Union{Symbol, Vector{Symbol}}`: Channel(s) to delete

# Returns
- `EegData`: A new EEG data object with channels deleted

# Examples
```julia
new_dat = channel_delete(dat, :Fp1)                 # Delete a single channel
new_dat = channel_delete(dat, [:Fp1, :Fp2, :Diode]) # Delete multiple channels
```

# Notes
- Only channels that exist in the data will be deleted
- Updates both the data columns and the layout rows
- Clears any cached neighbour information in the layout since channels have changed
- The original data is not modified
"""
function channel_delete(dat::EegData, channels_to_delete::Vector{Symbol})
    new_dat = copy(dat)
    channel_delete!(new_dat, channels_to_delete)
    return new_dat
end
channel_delete(dat::EegData, channel::Symbol) = channel_delete(dat, [channel])

function channel_delete(dat::Vector{<:EegData}, channels_to_delete::Union{Symbol,Vector{Symbol}})
    return [channel_delete(d, channels_to_delete) for d in dat]
end
