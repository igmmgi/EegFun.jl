# === CHANNEL RENAMING UTILITIES ===
"""
    rename_channel!(dat::EegData, rename_dict::Dict{Symbol, Symbol})

Rename channels in EEG data using a dictionary mapping old names to new names.
Modifies the data in place by updating both the data columns and the layout.

# Arguments
- `dat::EegData`: The EEG data object to modify
- `rename_dict::Dict{Symbol, Symbol}`: Dictionary mapping old channel names to new names

# Returns
- `nothing` (modifies the data in place)

# Examples
```julia
# Rename Fp1 to Fpz and Fp2 to Fpz
rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
rename_channel!(dat, rename_dict)

# Rename a single channel
rename_channel!(dat, Dict(:Cz => :Cz_new))
```

# Notes
- Only channels that exist in the data will be renamed
- If multiple channels would be renamed to the same name, an error is thrown to prevent duplicates
- Updates both the data columns and the layout labels
- Clears any cached neighbour information in the layout since channel names have changed
- Properly handles swaps (e.g., Dict(:A => :B, :B => :A) correctly exchanges the channels)
"""
function rename_channel!(dat::EegData, rename_dict::Dict{Symbol,Symbol})
    # Capture original channel names before renaming layout
    original_channels = Set(dat.layout.data.label)

    # First rename channels in the layout
    rename_channel!(dat.layout, rename_dict)

    # Get the list of channels that were actually renamed in the layout
    layout_channels = dat.layout.data.label
    channels_to_rename = keys(rename_dict)
    channels_found = intersect(original_channels, channels_to_rename)

    if isempty(channels_found)
        @info "rename_channel!: No channels found to rename in data"
        return dat
    end

    # Now rename the corresponding data columns using multiple dispatch
    _rename_data_columns!(dat, rename_dict, original_channels)

    @info "rename_channel!: Renamed $(length(channels_found)) channels in data and layout"
    return dat
end

"""Rename DataFrame columns using a two-phase swap-safe approach."""
function _rename_data_columns!(df::DataFrame, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    # Check for potential duplicate names before applying any renames
    final_names = Symbol[]
    for (old_name, new_name) in rename_dict
        if old_name ∈ existing_channels && old_name ∈ propertynames(df)
            push!(final_names, new_name)
        end
    end

    # Check for duplicates in final names
    if length(final_names) != length(unique(final_names))
        duplicate_names = filter(x -> count(==(x), final_names) > 1, unique(final_names))
        @minimal_error "Cannot rename channels to duplicate names: $(join(duplicate_names, ", "))"
    end

    # Apply the renaming with proper swap handling
    # First, collect all the final rename mappings to avoid interference
    final_renames = Dict{Symbol,Symbol}()  # old_name => final_name

    for (old_name, new_name) in rename_dict
        if old_name ∈ existing_channels && old_name ∈ propertynames(df)
            final_renames[old_name] = new_name
        end
    end

    # Nothing to do
    isempty(final_renames) && return

    # Two-phase rename to avoid collisions (e.g., swaps)
    # Phase 1: rename old names to unique temporary names
    temp_renames = Pair{Symbol,Symbol}[]
    used_names = Set(propertynames(df))  # currently present in df
    union!(used_names, Set(values(final_renames)))  # also avoid targeting final names

    for (old_name, new_name) in final_renames
        # Propose a unique temporary name
        base_tmp = Symbol(string(new_name), "__tmp__")
        tmp = base_tmp
        counter = 1
        while tmp ∈ used_names
            tmp = Symbol(string(base_tmp), "_", counter)
            counter += 1
        end
        push!(temp_renames, old_name => tmp)
        push!(used_names, tmp)
    end

    # Execute phase 1 renames
    for p in temp_renames
        rename!(df, p)
    end

    # Phase 2: rename temporary names to final names
    for (old_name, new_name) in final_renames
        # Find the temporary assigned to this old_name
        tmp = only(last.(filter(p -> first(p) == old_name, temp_renames)))
        rename!(df, tmp => new_name)
    end
end

function _rename_data_columns!(dat::SingleDataFrameEeg, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    _rename_data_columns!(dat.data, rename_dict, existing_channels)
end

function _rename_data_columns!(dat::MultiDataFrameEeg, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    # Use broadcasting to apply the DataFrame method to all DataFrames
    _rename_data_columns!.(dat.data, Ref(rename_dict), Ref(existing_channels))
end




"""
    rename_channel(dat::EegData, rename_dict::Dict{Symbol, Symbol})

Create a renamed copy of EEG data using a dictionary mapping old names to new names.

# Arguments
- `dat::EegData`: The EEG data object to rename
- `rename_dict::Dict{Symbol, Symbol}`: Dictionary mapping old channel names to new names

# Returns
- `EegData`: A new EEG data object with renamed channels

# Examples
```julia
# Rename Fp1 to Fpz and Fp2 to Fpz
rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
new_dat = rename_channel(dat, rename_dict)
```

# Notes
- Only channels that exist in the data will be renamed
- If multiple channels would be renamed to the same name, an error is thrown to prevent duplicates
- Updates both the data columns and the layout labels
- The original data is not modified
- Properly handles swaps (e.g., Dict(:A => :B, :B => :A) correctly exchanges the channels)
"""
function rename_channel(dat::EegData, rename_dict::Dict{Symbol,Symbol})
    return rename_channel!(copy(dat), rename_dict)
end
