"""
Data mirroring for reducing edge artifacts.
Provides functions to mirror (reverse and append) EEG data before and/or after epochs. 
"""

# === CORE MIRRORING FUNCTIONS ===

"""
    mirror!(dat::EpochData, side::Symbol = :both) -> Nothing
    mirror!(dat::ErpData, side::Symbol = :both) -> Nothing
    mirror!(data_vec::Vector{EpochData}, side::Symbol = :both) -> Nothing
    mirror!(data_vec::Vector{ErpData}, side::Symbol = :both) -> Nothing

Mirror data in-place by appending time-reversed copies before and/or after the data.
Reduces edge artifacts when filtering by creating smooth continuity at boundaries.
Always call [`unmirror!`](@ref) after processing to restore the original length.

# Arguments
- `dat` / `data_vec`: `EpochData`, `ErpData`, or a `Vector` of either
- `side`: `:pre`, `:post`, or `:both` (default: `:both`)

# Notes
- Mirrored data will be approximately 3× longer with `:both`
- The `side` argument to `unmirror!` must match

# Examples
```julia
mirror!(epochs, :both)
filter!(epochs, 1.0, filter_type = :hp)
unmirror!(epochs, :both)
```
"""
function mirror!(dat::EpochData, side::Symbol = :both)::Nothing

    # Validate side parameter
    if side ∉ [:pre, :post, :both]
        @minimal_error("side must be :pre, :post, or :both, got :$side")
    end
    @info "Mirroring epoched data on side: $side"

    for epoch_idx in eachindex(dat.data)
        dat.data[epoch_idx] = _mirror_dataframe!(dat.data[epoch_idx], side)
    end

    return nothing
end


"""
    mirror(dat::EpochData, side::Symbol = :both) -> EpochData
    mirror(dat::ErpData, side::Symbol = :both) -> ErpData
    mirror(data_vec::Vector{EpochData}, side::Symbol = :both) -> Vector{EpochData}
    mirror(data_vec::Vector{ErpData}, side::Symbol = :both) -> Vector{ErpData}

Non-mutating version of [`mirror!`](@ref). Returns a new object with mirrored data,
leaving the original unchanged.

# Arguments
- `dat` / `data_vec`: `EpochData`, `ErpData`, or a `Vector` of either
- `side`: `:pre`, `:post`, or `:both` (default: `:both`)
"""
function mirror(dat::EpochData, side::Symbol = :both)::EpochData
    dat_copy = copy(dat)
    mirror!(dat_copy, side)
    return dat_copy
end


function mirror!(dat::ErpData, side::Symbol = :both)::Nothing

    if side ∉ [:pre, :post, :both]
        @minimal_error("side must be :pre, :post, or :both, got :$side")
    end
    @info "Mirroring ERP data on side: $side"

    dat.data = _mirror_dataframe!(dat.data, side)

    return nothing
end


function mirror(dat::ErpData, side::Symbol = :both)::ErpData
    dat_copy = copy(dat)
    mirror!(dat_copy, side)
    return dat_copy
end


function mirror!(data_vec::Vector{EpochData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        mirror!(dat, side)
    end
    return nothing
end


function mirror(data_vec::Vector{EpochData}, side::Symbol = :both)::Vector{EpochData}
    return [mirror(dat, side) for dat in data_vec]
end


function mirror!(data_vec::Vector{ErpData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        mirror!(dat, side)
    end
    return nothing
end


function mirror(data_vec::Vector{ErpData}, side::Symbol = :both)::Vector{ErpData}
    return [mirror(dat, side) for dat in data_vec]
end


# === UNMIRRORING FUNCTIONS ===

"""
    unmirror!(dat::EpochData, side::Symbol = :both) -> Nothing
    unmirror!(dat::ErpData, side::Symbol = :both) -> Nothing
    unmirror!(data_vec::Vector{EpochData}, side::Symbol = :both) -> Nothing
    unmirror!(data_vec::Vector{ErpData}, side::Symbol = :both) -> Nothing

Remove mirrored sections from data in-place, restoring the original length.
Must be called with the same `side` as the preceding `mirror!` call.

# Arguments
- `dat` / `data_vec`: `EpochData`, `ErpData`, or a `Vector` of either
- `side`: `:pre`, `:post`, or `:both` (default: `:both`)

# Examples
```julia
mirror!(epochs, :both)
filter!(epochs, 1.0, filter_type = :hp)
unmirror!(epochs, :both)  # restores original length
```
"""
function unmirror!(dat::EpochData, side::Symbol = :both)::Nothing

    if side ∉ [:pre, :post, :both]
        @minimal_error("side must be :pre, :post, or :both, got :$side")
    end
    @info "Unmirroring epoched data on side: $side"

    for epoch_idx in eachindex(dat.data)
        dat.data[epoch_idx] = _unmirror_dataframe!(dat.data[epoch_idx], side)
    end

    return nothing
end


function unmirror(dat::EpochData, side::Symbol = :both)::EpochData
    dat_copy = copy(dat)
    unmirror!(dat_copy, side)
    return dat_copy
end


function unmirror!(dat::ErpData, side::Symbol = :both)::Nothing

    if side ∉ [:pre, :post, :both]
        @minimal_error("side must be :pre, :post, or :both, got :$side")
    end
    @info "Unmirroring ERP data on side: $side"

    dat.data = _unmirror_dataframe!(dat.data, side)

    return nothing
end


function unmirror(dat::ErpData, side::Symbol = :both)::ErpData
    dat_copy = copy(dat)
    unmirror!(dat_copy, side)
    return dat_copy
end


function unmirror!(data_vec::Vector{EpochData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        unmirror!(dat, side)
    end
    return nothing
end


"""
    unmirror(dat::EpochData, side::Symbol = :both) -> EpochData
    unmirror(dat::ErpData, side::Symbol = :both) -> ErpData
    unmirror(data_vec::Vector{EpochData}, side::Symbol = :both) -> Vector{EpochData}
    unmirror(data_vec::Vector{ErpData}, side::Symbol = :both) -> Vector{ErpData}

Non-mutating version of [`unmirror!`](@ref). Returns a new object with mirrored
sections removed, leaving the original unchanged.

# Arguments
- `dat` / `data_vec`: `EpochData`, `ErpData`, or a `Vector` of either
- `side`: `:pre`, `:post`, or `:both` (default: `:both`) — must match the original `mirror!` call
"""
function unmirror(data_vec::Vector{EpochData}, side::Symbol = :both)::Vector{EpochData}
    return [unmirror(dat, side) for dat in data_vec]
end


function unmirror!(data_vec::Vector{ErpData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        unmirror!(dat, side)
    end
    return nothing
end


function unmirror(data_vec::Vector{ErpData}, side::Symbol = :both)::Vector{ErpData}
    return [unmirror(dat, side) for dat in data_vec]
end


# === INTERNAL HELPER FUNCTIONS ===

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

        # Update df in-place by returning new and assigning at call site
        df_new = vcat(mirror_section, df)
        return df_new

    elseif side == :post
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)

        # Create mirrored section (reversed, excluding last point to avoid duplication)
        mirror_section = df[(end-1):-1:1, :]
        n_mirror = nrow(mirror_section)
        mirror_section.time = [df.time[end] + i * dt for i = 1:n_mirror]

        # Update df in-place by returning new and assigning at call site
        df_new = vcat(df, mirror_section)
        return df_new

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

        # Update df in-place by returning new and assigning at call site
        df_new = vcat(pre_mirror, df, post_mirror)
        return df_new
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
        return df_new

    elseif side == :post
        # After mirroring :post: total = original + (original-1) = 2*original - 1
        original_length = div(n_samples + 1, 2)

        df_new = df[1:original_length, :]
        return df_new

    else  # :both
        # After mirroring :both: total = (original-1) + original + (original-1) = 3*original - 2
        original_length = div(n_samples + 2, 3)
        pre_mirror_length = original_length - 1
        start_idx = pre_mirror_length + 1
        end_idx = start_idx + original_length - 1

        df_new = df[start_idx:end_idx, :]
        return df_new
    end
end

# === CONTINUOUS DATA SUPPORT ===

function mirror!(dat::ContinuousData, side::Symbol = :both)::Nothing
    if side ∉ [:pre, :post, :both]
        @minimal_error("side must be :pre, :post, or :both, got :$side")
    end
    @info "Mirroring continuous data on side: $side"
    dat.data = _mirror_dataframe!(dat.data, side)
    return nothing
end

function mirror!(data_vec::Vector{ContinuousData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        mirror!(dat, side)
    end
    return nothing
end

function unmirror!(dat::ContinuousData, side::Symbol = :both)::Nothing
    if side ∉ [:pre, :post, :both]
        @minimal_error("side must be :pre, :post, or :both, got :$side")
    end
    @info "Unmirroring continuous data on side: $side"
    dat.data = _unmirror_dataframe!(dat.data, side)
    return nothing
end

function unmirror!(data_vec::Vector{ContinuousData}, side::Symbol = :both)::Nothing
    for dat in data_vec
        unmirror!(dat, side)
    end
    return nothing
end

function mirror(dat::ContinuousData, side::Symbol = :both)::ContinuousData
    dat_copy = copy(dat)
    mirror!(dat_copy, side)
    return dat_copy
end

function mirror(data_vec::Vector{ContinuousData}, side::Symbol = :both)::Vector{ContinuousData}
    return [mirror(dat, side) for dat in data_vec]
end

function unmirror(dat::ContinuousData, side::Symbol = :both)::ContinuousData
    dat_copy = copy(dat)
    unmirror!(dat_copy, side)
    return dat_copy
end

function unmirror(data_vec::Vector{ContinuousData}, side::Symbol = :both)::Vector{ContinuousData}
    return [unmirror(dat, side) for dat in data_vec]
end
