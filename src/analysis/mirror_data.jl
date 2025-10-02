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
    mirror_data!(dat::EpochData, side::Symbol = :both)::Nothing

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

# How it works
For each epoch with data [a, b, c, d, e]:
- `:pre` → [e, d, c, b, a, b, c, d] (mirror prepended, excluding first point)
- `:post` → [b, c, d, e, d, c, b, a] (mirror appended, excluding last point)
- `:both` → [e, d, c, b, a, b, c, d, e, d, c, b, a] (both mirrors)

# Examples
```julia
using eegfun, JLD2

# Load epochs
epochs = load("participant_1_epochs.jld2", "epochs")

# Mirror on both sides (recommended for filtering)
mirror_data!(epochs, :both)

# Now filter the data (edges are protected)
filter!(epochs, 0.1, 30.0)

# Remove mirrored sections after filtering
unmirror_data!(epochs, :both)

# Continue with analysis
```

# Use Cases
- **Before filtering**: Reduce edge artifacts
- **Before baseline correction**: Ensure smooth baseline region
- **Any processing sensitive to edges**: FFT, wavelet analysis, etc.

# Important Notes
- Always call `unmirror_data!()` after processing to remove mirrored sections
- The `side` parameter for unmirroring must match the mirroring side
- Epochs will be approximately 3× longer with `:both` mirroring
- Sample rate is preserved
"""
function mirror_data!(dat::EpochData, side::Symbol = :both)::Nothing
    
    @info "Mirroring epoched data on side: $side"
    
    # Validate side parameter
    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    
    # Mirror each epoch
    for (epoch_idx, epoch) in enumerate(dat.data)
        _mirror_epoch!(epoch, side)
    end
    
    @info "Mirroring complete. $(length(dat.data)) epochs mirrored."
    @info "Remember to call unmirror_data!() after processing to remove mirrored sections."
    
    return nothing
end


"""
    mirror_data(dat::EpochData, side::Symbol = :both)::EpochData

Non-mutating version of mirror_data!. Returns new EpochData with mirrored epochs.
"""
function mirror_data(dat::EpochData, side::Symbol = :both)::EpochData
    # Create deep copy
    dat_copy = EpochData(
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info)
    )
    
    mirror_data!(dat_copy, side)
    
    return dat_copy
end


"""
    mirror_data!(dat::ErpData, side::Symbol = :both)::Nothing

Mirror ERP data in-place by appending reversed data sections.

For averaged ERP data, mirroring extends the time series with time-reversed copies.
This is useful before filtering averaged data to reduce edge artifacts.

# Examples
```julia
# Load averaged ERP
erp = load("participant_1_erp.jld2", "erp")

# Mirror, filter, then unmirror
mirror_data!(erp, :both)
filter!(erp, 0.1, 30.0)
unmirror_data!(erp, :both)
```
"""
function mirror_data!(dat::ErpData, side::Symbol = :both)::Nothing
    
    @info "Mirroring ERP data on side: $side"
    
    # Validate side parameter
    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    
    _mirror_dataframe!(dat.data, side)
    
    @info "Mirroring complete."
    @info "Remember to call unmirror_data!() after processing to remove mirrored sections."
    
    return nothing
end


"""
    mirror_data(dat::ErpData, side::Symbol = :both)::ErpData

Non-mutating version of mirror_data! for ERP data.
"""
function mirror_data(dat::ErpData, side::Symbol = :both)::ErpData
    # Create copy
    dat_copy = ErpData(
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
        dat.n_epochs
    )
    
    mirror_data!(dat_copy, side)
    
    return dat_copy
end


#=============================================================================
    UNMIRRORING FUNCTIONS
=============================================================================#

"""
    unmirror_data!(dat::EpochData, side::Symbol = :both)::Nothing

Remove mirrored sections from epoched data in-place.

This function removes the mirrored sections that were added by `mirror_data!()`,
restoring the data to its original length. Must be called with the same `side`
parameter as was used for mirroring.

# Arguments
- `dat::EpochData`: Mirrored epoched data
- `side::Symbol`: Which side(s) were mirrored - must match original mirroring

# Examples
```julia
# Mirror, process, then unmirror
mirror_data!(epochs, :both)
filter!(epochs, 0.1, 30.0)
unmirror_data!(epochs, :both)
```
"""
function unmirror_data!(dat::EpochData, side::Symbol = :both)::Nothing
    
    @info "Unmirroring epoched data on side: $side"
    
    # Validate side parameter
    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    
    # Unmirror each epoch
    for (epoch_idx, epoch) in enumerate(dat.data)
        _unmirror_epoch!(epoch, side)
    end
    
    @info "Unmirroring complete. Data restored to original length."
    
    return nothing
end


"""
    unmirror_data(dat::EpochData, side::Symbol = :both)::EpochData

Non-mutating version of unmirror_data!.
"""
function unmirror_data(dat::EpochData, side::Symbol = :both)::EpochData
    # Create deep copy
    dat_copy = EpochData(
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info)
    )
    
    unmirror_data!(dat_copy, side)
    
    return dat_copy
end


"""
    unmirror_data!(dat::ErpData, side::Symbol = :both)::Nothing

Remove mirrored sections from ERP data in-place.
"""
function unmirror_data!(dat::ErpData, side::Symbol = :both)::Nothing
    
    @info "Unmirroring ERP data on side: $side"
    
    # Validate side parameter
    if side ∉ [:pre, :post, :both]
        @minimal_error_throw("side must be :pre, :post, or :both, got :$side")
    end
    
    _unmirror_dataframe!(dat.data, side)
    
    @info "Unmirroring complete. Data restored to original length."
    
    return nothing
end


"""
    unmirror_data(dat::ErpData, side::Symbol = :both)::ErpData

Non-mutating version of unmirror_data! for ERP data.
"""
function unmirror_data(dat::ErpData, side::Symbol = :both)::ErpData
    # Create copy
    dat_copy = ErpData(
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
        dat.n_epochs
    )
    
    unmirror_data!(dat_copy, side)
    
    return dat_copy
end


#=============================================================================
    INTERNAL HELPER FUNCTIONS
=============================================================================#

"""
Mirror a single epoch DataFrame in-place.
"""
function _mirror_epoch!(epoch::DataFrame, side::Symbol)
    n_samples = nrow(epoch)
    epoch_duration = epoch.time[end] - epoch.time[1]
    
    if side == :pre
        # Create mirrored section (reversed, excluding last point to avoid duplication)
        mirror_section = epoch[end-1:-1:1, :]
        
        # Calculate time step
        dt = (epoch.time[end] - epoch.time[1]) / (n_samples - 1)
        
        # Set mirror times to extend backwards from first time point
        n_mirror = nrow(mirror_section)
        mirror_section.time = [epoch.time[1] - (n_mirror - i + 1) * dt for i in 1:n_mirror]
        
        # Concatenate: mirror + original (all points)
        epoch_new = vcat(mirror_section, epoch)
        empty!(epoch)
        append!(epoch, epoch_new)
        
    elseif side == :post
        # Create mirrored section (reversed, excluding first point to avoid duplication)
        mirror_section = epoch[end:-1:2, :]
        
        # Calculate time step
        dt = (epoch.time[end] - epoch.time[1]) / (n_samples - 1)
        
        # Set mirror times to extend forwards from last time point
        n_mirror = nrow(mirror_section)
        mirror_section.time = [epoch.time[end] + i * dt for i in 1:n_mirror]
        
        # Concatenate: original (all points) + mirror
        epoch_new = vcat(epoch, mirror_section)
        empty!(epoch)
        append!(epoch, epoch_new)
        
    else  # :both
        # Calculate time step
        dt = (epoch.time[end] - epoch.time[1]) / (n_samples - 1)
        
        # Create pre-mirror section (reversed, excluding last point to avoid duplication)
        pre_mirror = epoch[end-1:-1:1, :]
        n_pre = nrow(pre_mirror)
        pre_mirror.time = [epoch.time[1] - (n_pre - i + 1) * dt for i in 1:n_pre]
        
        # Create post-mirror section (reversed, excluding first point to avoid duplication)
        post_mirror = epoch[end:-1:2, :]
        n_post = nrow(post_mirror)
        post_mirror.time = [epoch.time[end] + i * dt for i in 1:n_post]
        
        # Concatenate: pre_mirror + original (all points) + post_mirror
        epoch_new = vcat(pre_mirror, epoch, post_mirror)
        empty!(epoch)
        append!(epoch, epoch_new)
    end
end


"""
Mirror a DataFrame (for ERP data) in-place.
"""
function _mirror_dataframe!(df::DataFrame, side::Symbol)
    n_samples = nrow(df)
    
    if side == :pre
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)
        
        # Create mirrored section (reversed, excluding last point to avoid duplication)
        mirror_section = df[end-1:-1:1, :]
        n_mirror = nrow(mirror_section)
        mirror_section.time = [df.time[1] - (n_mirror - i + 1) * dt for i in 1:n_mirror]
        
        # Update df in-place
        df_new = vcat(mirror_section, df)
        empty!(df)
        append!(df, df_new)
        
    elseif side == :post
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)
        
        # Create mirrored section (reversed, excluding first point to avoid duplication)
        mirror_section = df[end:-1:2, :]
        n_mirror = nrow(mirror_section)
        mirror_section.time = [df.time[end] + i * dt for i in 1:n_mirror]
        
        # Update df in-place
        df_new = vcat(df, mirror_section)
        empty!(df)
        append!(df, df_new)
        
    else  # :both
        # Calculate time step
        dt = (df.time[end] - df.time[1]) / (n_samples - 1)
        
        # Create pre-mirror section (reversed, excluding last point to avoid duplication)
        pre_mirror = df[end-1:-1:1, :]
        n_pre = nrow(pre_mirror)
        pre_mirror.time = [df.time[1] - (n_pre - i + 1) * dt for i in 1:n_pre]
        
        # Create post-mirror section (reversed, excluding first point to avoid duplication)
        post_mirror = df[end:-1:2, :]
        n_post = nrow(post_mirror)
        post_mirror.time = [df.time[end] + i * dt for i in 1:n_post]
        
        # Update df in-place
        df_new = vcat(pre_mirror, df, post_mirror)
        empty!(df)
        append!(df, df_new)
    end
end


"""
Unmirror a single epoch DataFrame in-place.
"""
function _unmirror_epoch!(epoch::DataFrame, side::Symbol)
    n_samples = nrow(epoch)
    
    if side == :pre
        # After mirroring :pre: total = (original-1) + original = 2*original - 1
        # original = (n_samples + 1) / 2
        original_length = div(n_samples + 1, 2)
        mirror_length = original_length - 1
        
        # Keep only the last original_length rows (the original data)
        start_idx = mirror_length + 1
        epoch_new = epoch[start_idx:end, :]
        empty!(epoch)
        append!(epoch, epoch_new)
        
    elseif side == :post
        # After mirroring :post: total = original + (original-1) = 2*original - 1
        # original = (n_samples + 1) / 2
        original_length = div(n_samples + 1, 2)
        
        # Keep only the first original_length rows (the original data)
        epoch_new = epoch[1:original_length, :]
        empty!(epoch)
        append!(epoch, epoch_new)
        
    else  # :both
        # After mirroring :both: total = (original-1) + original + (original-1) = 3*original - 2
        # original = (n_samples + 2) / 3
        original_length = div(n_samples + 2, 3)
        pre_mirror_length = original_length - 1
        
        # Keep middle section (the original data)
        start_idx = pre_mirror_length + 1
        end_idx = start_idx + original_length - 1
        epoch_new = epoch[start_idx:end_idx, :]
        empty!(epoch)
        append!(epoch, epoch_new)
    end
end


"""
Unmirror a DataFrame (for ERP data) in-place.
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

