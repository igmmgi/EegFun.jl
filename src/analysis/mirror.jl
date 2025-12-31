"""
Data mirroring for reducing edge artifacts.

This module provides functions to mirror (reverse and append) EEG data before and/or
after epochs. This is commonly used before filtering to reduce edge artifacts by
creating smooth transitions at the epoch boundaries.
"""

#=============================================================================
    CORE MIRRORING FUNCTIONS
=============================================================================#

"""
    mirror!(dat::EpochData, side::Symbol = :both)::Nothing

Mirror epoched data in-place by appending reversed data sections.

Data mirroring adds reversed (time-flipped) copies of the data before and/or after
each epoch. This reduces edge artifacts when filtering by creating smooth continuity
at epoch boundaries.

# Arguments
- `dat::EpochData`: Epoched data to mirror
- `side::Symbol`: Which side(s) to mirror - `:pre`, `:post`, or `:both` (default: `:both`)

# Effects
- Modifies epochs in-place, extending their length
- Updates time vectors to reflect the extended data
- Original data is preserved in the middle section

# Examples
```julia
using eegfun, JLD2

# Load epochs
epochs = load("participant_1_epochs.jld2", "epochs")

# Mirror on both sides (recommended for filtering)
mirror!(epochs, :both)

# Now filter the data (edges are protected)
filter!(epochs, 1, filter_type = "hp")

# Remove mirrored sections after filtering
unmirror!(epochs, :both)

# Continue with analysis
```

# Use Cases
- **Before filtering**: Reduce edge artifacts
- **Any processing sensitive to edges**: FFT, wavelet analysis, etc.

# Important Notes
- Always call `unmirror!()` after processing to remove mirrored sections
- The `side` parameter for unmirroring must match the mirroring side
- Mirrored epochs will be approximately 3× longer with `:both` mirroring
- Sample rate is preserved
"""
function mirror!(dat::EpochData, side::Symbol = :both)::Nothing

    # Validate side parameter
    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    @info "Mirroring epoched data on side: $side"

    for (epoch_idx, epoch) in enumerate(dat.data)
        _mirror_dataframe!(epoch, side)
    end

    return nothing
end


"""
    mirror(dat::EpochData, side::Symbol = :both)::EpochData

Non-mutating version of mirror!. Returns new EpochData with mirrored epochs.
"""
function mirror(dat::EpochData, side::Symbol = :both)::EpochData
    dat_copy = copy(dat)
    mirror!(dat_copy, side)
    return dat_copy
end


"""
    mirror!(dat::ErpData, side::Symbol = :both)::Nothing

Mirror ERP data in-place by appending reversed data sections.

For averaged ERP data, mirroring extends the time series with time-reversed copies.
This is useful before filtering averaged data to reduce edge artifacts.

# Examples
```julia
# Load some data
erp = load("participant_1_erp.jld2", "erp")

# Mirror, filter, then unmirror
mirror!(erp, :both)
filter!(erp, 1, filter_type = "hp")
unmirror!(erp, :both)
```
"""
function mirror!(dat::ErpData, side::Symbol = :both)::Nothing

    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    @info "Mirroring ERP data on side: $side"

    _mirror_dataframe!(dat.data, side)

    return nothing
end


"""
    mirror(dat::ErpData, side::Symbol = :both)::ErpData

Non-mutating version of mirror! for ERP data.
"""
function mirror(dat::ErpData, side::Symbol = :both)::ErpData
    dat_copy = copy(dat)
    mirror!(dat_copy, side)
    return dat_copy
end


"""
    mirror!(data_vec::Vector{EpochData}, side::Symbol = :both)::Nothing

Mutating version of mirror for vector of EpochData objects.

# Arguments
- `data_vec::Vector{EpochData}`: Vector of EpochData objects to mirror (modified in-place)
- `side::Symbol`: Which side to mirror (`:left`, `:right`, or `:both`, default: `:both`)

# Returns
- `Nothing`: All objects in the vector are modified in-place

# Examples
```julia
# Mirror multiple EpochData objects
mirror!(epochs_vector, :both)  # Mirror both sides
```
"""
function mirror!(data_vec::Vector{EpochData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        mirror!(dat, side)
    end
    return nothing
end


"""
    mirror(data_vec::Vector{EpochData}, side::Symbol = :both)::Vector{EpochData}

Non-mutating version of mirror for vector of EpochData objects.

# Arguments
- `data_vec::Vector{EpochData}`: Vector of EpochData objects to mirror (NOT modified)
- `side::Symbol`: Which side to mirror (`:left`, `:right`, or `:both`, default: `:both`)

# Returns
- `Vector{EpochData}`: New vector with mirrored EpochData objects

# Examples
```julia
# Mirror multiple EpochData objects (creates new objects)
mirrored_epochs = mirror(epochs_vector, :both)
```
"""
function mirror(data_vec::Vector{EpochData}, side::Symbol = :both)::Vector{EpochData}
    return [mirror(dat, side) for dat in data_vec]
end


#=============================================================================
    UNMIRRORING FUNCTIONS
=============================================================================#

"""
    unmirror!(dat::EpochData, side::Symbol = :both)::Nothing

Remove mirrored sections from epoched data in-place.

This function removes the mirrored sections that were added by `mirror!()`,
restoring the data to its original length. Must be called with the same `side`
parameter as was used for mirroring.

# Arguments
- `dat::EpochData`: Mirrored epoched data
- `side::Symbol`: Which side(s) were mirrored - must match original mirroring

# Examples
```julia
# Mirror, process, then unmirror
mirror!(epochs, :both)
filter!(epochs, 1, filter_type = "hp")
unmirror!(epochs, :both)
```
"""
function unmirror!(dat::EpochData, side::Symbol = :both)::Nothing

    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    @info "Unmirroring epoched data on side: $side"

    # Unmirror each epoch
    for (epoch_idx, epoch) in enumerate(dat.data)
        _unmirror_dataframe!(epoch, side)
    end

    return nothing
end


"""
    unmirror(dat::EpochData, side::Symbol = :both)::EpochData

Non-mutating version of unmirror!.
"""
function unmirror(dat::EpochData, side::Symbol = :both)::EpochData
    dat_copy = copy(dat)
    unmirror!(dat_copy, side)
    return dat_copy
end


"""
    unmirror!(dat::ErpData, side::Symbol = :both)::Nothing

Remove mirrored sections from ERP data in-place.
"""
function unmirror!(dat::ErpData, side::Symbol = :both)::Nothing

    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    @info "Unmirroring ERP data on side: $side"

    _unmirror_dataframe!(dat.data, side)

    return nothing
end


"""
    unmirror(dat::ErpData, side::Symbol = :both)::ErpData

Non-mutating version of unmirror! for ERP data.
"""
function unmirror(dat::ErpData, side::Symbol = :both)::ErpData
    dat_copy = copy(dat)
    unmirror!(dat_copy, side)
    return dat_copy
end


"""
    unmirror!(data_vec::Vector{EpochData}, side::Symbol = :both)::Nothing

Mutating version of unmirror for vector of EpochData objects.

# Arguments
- `data_vec::Vector{EpochData}`: Vector of EpochData objects to unmirror (modified in-place)
- `side::Symbol`: Which side to unmirror (`:left`, `:right`, or `:both`, default: `:both`)

# Returns
- `Nothing`: All objects in the vector are modified in-place

# Examples
```julia
# Unmirror multiple EpochData objects
unmirror!(epochs_vector, :both)  # Unmirror both sides
```
"""
function unmirror!(data_vec::Vector{EpochData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        unmirror!(dat, side)
    end
    return nothing
end


"""
    unmirror(data_vec::Vector{EpochData}, side::Symbol = :both)::Vector{EpochData}

Non-mutating version of unmirror for vector of EpochData objects.

# Arguments
- `data_vec::Vector{EpochData}`: Vector of EpochData objects to unmirror (NOT modified)
- `side::Symbol`: Which side to unmirror (`:left`, `:right`, or `:both`, default: `:both`)

# Returns
- `Vector{EpochData}`: New vector with unmirrored EpochData objects

# Examples
```julia
# Unmirror multiple EpochData objects (creates new objects)
unmirrored_epochs = unmirror(epochs_vector, :both)
```
"""
function unmirror(data_vec::Vector{EpochData}, side::Symbol = :both)::Vector{EpochData}
    return [unmirror(dat, side) for dat in data_vec]
end


#=============================================================================
    INTERNAL HELPER FUNCTIONS
=============================================================================#

"""
Mirror a DataFrame in-place by reflecting the data at the boundaries.
Used for both epoch and continuous/ERP data.
"""
function _mirror_dataframe!(df::DataFrame, side::Symbol)
    n_samples = nrow(df)

    if side == :pre
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)

        # Create mirrored section (reversed, excluding first point to avoid duplication)
        mirror_section = df[end:-1:2, :]
        n_mirror = nrow(mirror_section)
        mirror_section.time = [df.time[1] - (n_mirror - i + 1) * dt for i = 1:n_mirror]

        # Update df in-place
        df_new = vcat(mirror_section, df)
        empty!(df)
        append!(df, df_new)

    elseif side == :post
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)

        # Create mirrored section (reversed, excluding last point to avoid duplication)
        mirror_section = df[(end-1):-1:1, :]
        n_mirror = nrow(mirror_section)
        mirror_section.time = [df.time[end] + i * dt for i = 1:n_mirror]

        # Update df in-place
        df_new = vcat(df, mirror_section)
        empty!(df)
        append!(df, df_new)

    else  # :both
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)

        # Create pre-mirror section (reversed, excluding first point to avoid duplication)
        pre_mirror = df[end:-1:2, :]
        n_pre = nrow(pre_mirror)
        pre_mirror.time = [df.time[1] - (n_pre - i + 1) * dt for i = 1:n_pre]

        # Create post-mirror section (reversed, excluding last point to avoid duplication)
        post_mirror = df[(end-1):-1:1, :]
        n_post = nrow(post_mirror)
        post_mirror.time = [df.time[end] + i * dt for i = 1:n_post]

        # Update df in-place
        df_new = vcat(pre_mirror, df, post_mirror)
        empty!(df)
        append!(df, df_new)
    end
end


"""
Unmirror a DataFrame in-place by removing the mirrored sections.
Used for both epoch and continuous/ERP data.
"""
function _unmirror_dataframe!(df::DataFrame, side::Symbol)
    n_samples = nrow(df)

    if side == :pre
        # After mirroring :pre: total = (original-1) + original = 2*original - 1
        original_length = div(n_samples + 1, 2)
        mirror_length = original_length - 1
        start_idx = mirror_length + 1

        df_new = df[start_idx:end, :]
        empty!(df)
        append!(df, df_new)

    elseif side == :post
        # After mirroring :post: total = original + (original-1) = 2*original - 1
        original_length = div(n_samples + 1, 2)

        df_new = df[1:original_length, :]
        empty!(df)
        append!(df, df_new)

    else  # :both
        # After mirroring :both: total = (original-1) + original + (original-1) = 3*original - 2
        original_length = div(n_samples + 2, 3)
        pre_mirror_length = original_length - 1
        start_idx = pre_mirror_length + 1
        end_idx = start_idx + original_length - 1

        df_new = df[start_idx:end_idx, :]
        empty!(df)
        append!(df, df_new)
    end
end
