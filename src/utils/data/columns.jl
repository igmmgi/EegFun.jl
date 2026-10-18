# === COLUMN IDENTIFICATION SYSTEM ===
"""
    _get_cols_by_group(dat::EegData, group::Symbol) -> Vector{Symbol}

Get columns by group type for EegData objects using layout-based identification.

# Arguments
- `dat::EegData`: The EEG data object
- `group::Symbol`: The group type (:metadata, :channels, :extra)

# Returns
- `Vector{Symbol}`: Column names of the specified group

# Group Types
- `:channels`: EEG channel columns (intersection of layout labels and DataFrame columns)
- `:metadata`: System columns (all columns before first layout label)
- `:extra`: Derived columns (all columns after last layout label)
"""
function _get_cols_by_group(dat::EegData, group::Symbol)

    if group ∉ (:channels, :metadata, :extra)
        @minimal_error "Unknown group type: $group"
    end

    labels = all_labels(dat)
    layout_channels = dat.layout.data.label

    if group == :channels
        return intersect(layout_channels, labels)
    elseif group == :metadata
        isempty(layout_channels) && return Symbol[]
        valid_channels = intersect(layout_channels, labels)
        isempty(valid_channels) && return Symbol[]
        first_channel_idx = findfirst(==(valid_channels[1]), labels)
        return labels[1:(first_channel_idx-1)]
    elseif group == :extra
        isempty(layout_channels) && return Symbol[]
        valid_channels = intersect(layout_channels, labels)
        isempty(valid_channels) && return Symbol[]
        last_channel_idx = findlast(==(valid_channels[end]), labels)
        return labels[(last_channel_idx+1):end]
    end
end


"""
    _get_cols_by_group(dat::Vector{<:EegData}, group::Symbol) -> Vector{Symbol}

Get columns by group type for a vector of EegData objects.
Delegates to the first element (assumes all have the same structure).

# Arguments
- `dat::Vector{<:EegData}`: Vector of EEG data objects
- `group::Symbol`: The group type (:metadata, :channels, :extra)

# Returns
- `Vector{Symbol}`: Column names of the specified group
"""
function _get_cols_by_group(dat::Vector{<:EegData}, group::Symbol)
    isempty(dat) && return Symbol[]
    return _get_cols_by_group(first(dat), group)
end
