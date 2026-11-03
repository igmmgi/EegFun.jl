const _IcaDataTarget = Union{DataFrame,Vector{DataFrame},ContinuousData,EpochData,AbstractVector{<:EpochData}}

_copy_ica_target(dat::Union{DataFrame,ContinuousData,EpochData}) = copy(dat)
_copy_ica_target(dats::AbstractVector) = [copy(d) for d in dats]

@inline function _pack_col_segment!(col_out::AbstractVector{Float64}, offset::Int, src::AbstractVector)
    n_rows = length(src)
    @inbounds @simd for s = 1:n_rows
        col_out[offset+s] = Float64(src[s])
    end
    return offset + n_rows
end

@inline function _unpack_col_segment!(dst::AbstractVector, col_in::AbstractVector{Float64}, offset::Int)
    n_rows = length(dst)
    @inbounds @simd for s = 1:n_rows
        dst[s] = col_in[offset+s]
    end
    return offset + n_rows
end

function _pack_df_to_matrix(df::DataFrame, labels::Vector{Symbol})
    X = Matrix{Float64}(undef, nrow(df), length(labels))
    for (ch_idx, ch) in enumerate(labels)
        _pack_col_segment!(view(X, :, ch_idx), 0, df[!, ch])
    end
    return X
end

function _unpack_matrix_to_df!(df::DataFrame, X::Matrix{Float64}, labels::Vector{Symbol})
    for (ch_idx, ch) in enumerate(labels)
        _unpack_col_segment!(df[!, ch], view(X, :, ch_idx), 0)
    end
    return nothing
end

function _pack_dfs_to_matrix(dfs::Vector{DataFrame}, labels::Vector{Symbol}, n_total_samples::Int)
    X = Matrix{Float64}(undef, n_total_samples, length(labels))
    for (ch_idx, ch) in enumerate(labels)
        col_out = view(X, :, ch_idx)
        offset = 0
        for df in dfs
            offset = _pack_col_segment!(col_out, offset, df[!, ch])
        end
    end
    return X
end

function _unpack_matrix_to_dfs!(dfs::Vector{DataFrame}, X::Matrix{Float64}, labels::Vector{Symbol})
    for (ch_idx, ch) in enumerate(labels)
        col_in = view(X, :, ch_idx)
        offset = 0
        for df in dfs
            offset = _unpack_col_segment!(df[!, ch], col_in, offset)
        end
    end
    return nothing
end

function _collect_validated_epoch_dfs(dats::AbstractVector{<:EpochData})
    all_dfs = DataFrame[]
    for d in dats
        if !isempty(d.data)
            n_samples_cond = nrow(d.data[1])
            all(nrow(df) == n_samples_cond for df in d.data) ||
                @minimal_error("All epochs within an EpochData condition must have the same number of samples ($n_samples_cond).")
            append!(all_dfs, d.data)
        end
    end
    return all_dfs
end

function _subtract_core!(X::Matrix{Float64}, ica::InfoIca, comps::Vector{Int})
    Uc = ica.unmixing[comps, :]
    Mc = ica.mixing[:, comps]

    acts = Uc * transpose(X)
    mean_proj = Uc * ica.mean
    inv_scale = 1.0 / ica.scale
    n_k, n_samples = size(acts)
    @inbounds for j = 1:n_samples
        @simd for i = 1:n_k
            acts[i, j] = (acts[i, j] - mean_proj[i]) * inv_scale
        end
    end

    mul!(X, transpose(acts), transpose(Mc), -ica.scale, 1.0)
    return acts
end

function _add_core!(X::Matrix{Float64}, ica::InfoIca, comps::Vector{Int})
    n_samples = size(X, 1)
    n_k = length(comps)
    acts = Matrix{Float64}(undef, n_k, n_samples)
    for (i, c) in enumerate(comps)
        stored = ica.removed_activations[c]
        size(stored, 2) == n_samples ||
            @minimal_error("Stored activations have $(size(stored, 2)) samples, but data has $n_samples samples.")
        @inbounds @simd for s = 1:n_samples
            acts[i, s] = stored[1, s]
        end
    end

    Mc = ica.mixing[:, comps]
    mul!(X, transpose(acts), transpose(Mc), ica.scale, 1.0)
    return nothing
end

function _validate_and_filter_subtract_components(ica::InfoIca, component_selection::Union{Function,Nothing})
    isnothing(component_selection) && @minimal_error(
        "subtract_ica_components! requires an explicit 'component_selection' (e.g., component_selection = components([1, 2]) or components(artifacts)). To subtract all components, pass component_selection = components().",
    )

    comps = sort!(unique(get_selected_components(ica, component_selection)))
    n_components = size(ica.unmixing, 1)
    all(1 .<= comps .<= n_components) || @minimal_error("Selected components must be between 1 and $n_components")

    already_removed = filter(c -> haskey(ica.removed_activations, c), comps)
    if !isempty(already_removed)
        @minimal_warning "Components $(already_removed) have already been removed from this ICA model and are in removed_activations. Skipping to avoid subtracting zeros or corrupting stored activations. Use a fresh copy(ica) or clear ica.removed_activations if applying to a new dataset."
        comps = filter(c -> !haskey(ica.removed_activations, c), comps)
    end
    return comps
end

function _validate_and_filter_add_components(ica::InfoIca, component_selection::Union{Function,Nothing})
    n_components = size(ica.unmixing, 1)
    comps = if isnothing(component_selection)
        sort(collect(keys(ica.removed_activations)))
    else
        sort!(unique(get_selected_components(ica, component_selection)))
    end

    all(1 .<= comps .<= n_components) || @minimal_error("Components must be between 1 and $n_components")

    if length(comps) == n_components && comps == collect(1:n_components)
        comps = sort(collect(keys(ica.removed_activations)))
    else
        for comp in comps
            haskey(ica.removed_activations, comp) || @minimal_error("Component $comp has no stored activations to restore.")
        end
    end

    if isempty(comps)
        @minimal_warning "No stored removed activations found to restore."
        return Int[]
    end
    return comps
end

function _store_removed_activations!(ica::InfoIca, comps::Vector{Int}, acts::Matrix{Float64})
    for (i, comp_idx) in enumerate(comps)
        ica.removed_activations[comp_idx] = acts[i:i, :]
    end
    return nothing
end

function _subtract_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    comps = _validate_and_filter_subtract_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    X = _pack_df_to_matrix(dat, labels)
    acts = _subtract_core!(X, ica, comps)
    _store_removed_activations!(ica, comps, acts)
    _unpack_matrix_to_df!(dat, X, labels)
    return comps
end

function _subtract_ica_components_dfs!(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function,Nothing} = nothing,
    require_equal_rows::Bool = false,
)
    isempty(dfs) && return Int[]
    if require_equal_rows
        n_rows_1 = nrow(dfs[1])
        all(nrow(df) == n_rows_1 for df in dfs) ||
            @minimal_error("All epochs in Vector{DataFrame} must have the same number of samples ($n_rows_1).")
    end

    comps = _validate_and_filter_subtract_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    X = _pack_dfs_to_matrix(dfs, labels, sum(nrow, dfs))
    acts = _subtract_core!(X, ica, comps)
    _store_removed_activations!(ica, comps, acts)
    _unpack_matrix_to_dfs!(dfs, X, labels)
    return comps
end

_subtract_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing) =
    _subtract_ica_components_dfs!(dfs, ica; component_selection = component_selection, require_equal_rows = true)

function _subtract_ica_components_multi!(
    dats::AbstractVector{<:EpochData},
    ica::InfoIca;
    component_selection::Union{Function,Nothing} = nothing,
)
    all_dfs = _collect_validated_epoch_dfs(dats)
    return _subtract_ica_components_dfs!(all_dfs, ica; component_selection = component_selection, require_equal_rows = false)
end

function _add_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    comps = _validate_and_filter_add_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    X = _pack_df_to_matrix(dat, labels)
    _add_core!(X, ica, comps)
    _unpack_matrix_to_df!(dat, X, labels)
    for comp in comps
        delete!(ica.removed_activations, comp)
    end
    return comps
end

function _add_ica_components_dfs!(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function,Nothing} = nothing,
    require_equal_rows::Bool = false,
)
    isempty(dfs) && return Int[]
    if require_equal_rows
        n_rows_1 = nrow(dfs[1])
        all(nrow(df) == n_rows_1 for df in dfs) ||
            @minimal_error("All epochs in Vector{DataFrame} must have the same number of samples ($n_rows_1).")
    end

    comps = _validate_and_filter_add_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    X = _pack_dfs_to_matrix(dfs, labels, sum(nrow, dfs))
    _add_core!(X, ica, comps)
    _unpack_matrix_to_dfs!(dfs, X, labels)
    for comp in comps
        delete!(ica.removed_activations, comp)
    end
    return comps
end

_add_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing) =
    _add_ica_components_dfs!(dfs, ica; component_selection = component_selection, require_equal_rows = true)

function _add_ica_components_multi!(dats::AbstractVector{<:EpochData}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    all_dfs = _collect_validated_epoch_dfs(dats)
    return _add_ica_components_dfs!(all_dfs, ica; component_selection = component_selection, require_equal_rows = false)
end

"""
    subtract_ica_components!(dat, ica::InfoIca; component_selection = nothing)

Remove selected ICA components in-place from a `DataFrame`, `Vector{DataFrame}`, `ContinuousData`,
`EpochData`, or `AbstractVector{<:EpochData}`, storing removed activations in `ica.removed_activations`
for exact restoration via [`add_ica_components!`](@ref).

An explicit `component_selection` (e.g. `components([1, 3])` or `components()`) is required.
"""
function subtract_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    _subtract_ica_components!(dat, ica; component_selection = component_selection)
    return dat
end

function subtract_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    _subtract_ica_components!(dfs, ica; component_selection = component_selection)
    return dfs
end

function subtract_ica_components!(
    dat::Union{ContinuousData,EpochData},
    ica::InfoIca;
    component_selection::Union{Function,Nothing} = nothing,
)
    actually_removed = _subtract_ica_components!(dat.data, ica; component_selection = component_selection)
    dat.analysis_info.n_ica_components_removed += length(actually_removed)
    return dat
end

function subtract_ica_components!(dats::AbstractVector{<:EpochData}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    isempty(dats) && return dats
    actually_removed = _subtract_ica_components_multi!(dats, ica; component_selection = component_selection)
    for d in dats
        d.analysis_info.n_ica_components_removed += length(actually_removed)
    end
    return dats
end

"""
    subtract_ica_components(dat, ica::InfoIca; component_selection = nothing) -> (cleaned_data, ica_updated)

Non-mutating ICA component removal; copies both `dat` and `ica` and returns `(cleaned_data, ica_updated)`.
"""
function subtract_ica_components(dat::_IcaDataTarget, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    dat_out = _copy_ica_target(dat)
    ica_out = copy(ica)
    subtract_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end

"""
    add_ica_components!(dat, ica::InfoIca; component_selection = nothing)

Restore previously subtracted ICA components in-place using activations stored in
`ica.removed_activations` (restores all removed components when `component_selection === nothing`).
"""
function add_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    _add_ica_components!(dat, ica; component_selection = component_selection)
    return dat
end

function add_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    _add_ica_components!(dfs, ica; component_selection = component_selection)
    return dfs
end

function add_ica_components!(dat::Union{ContinuousData,EpochData}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    actually_added = _add_ica_components!(dat.data, ica; component_selection = component_selection)
    dat.analysis_info.n_ica_components_removed -= length(actually_added)
    return dat
end

function add_ica_components!(dats::AbstractVector{<:EpochData}, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    isempty(dats) && return dats
    actually_added = _add_ica_components_multi!(dats, ica; component_selection = component_selection)
    for d in dats
        d.analysis_info.n_ica_components_removed -= length(actually_added)
    end
    return dats
end

"""
    add_ica_components(dat, ica::InfoIca; component_selection = nothing) -> (restored_data, ica_updated)

Non-mutating restoration of previously removed ICA components; returns `(restored_data, ica_updated)`.
"""
function add_ica_components(dat::_IcaDataTarget, ica::InfoIca; component_selection::Union{Function,Nothing} = nothing)
    dat_out = _copy_ica_target(dat)
    ica_out = copy(ica)
    add_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end

function _copy_col_to_matrix!(dst::Matrix{Float64}, col_idx::Int, src::AbstractVector, selected_samples::Vector{Int})
    isempty(selected_samples) && return dst
    min_idx, max_idx = extrema(selected_samples)
    if min_idx < 1 || max_idx > length(src)
        @minimal_error "Sample indices out of bounds: range [$min_idx, $max_idx] exceeds column length $(length(src))"
    end
    @inbounds for (j, sample_idx) in enumerate(selected_samples)
        dst[col_idx, j] = Float64(src[sample_idx])
    end
    return dst
end

function _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca, selected_samples::Vector{Int})
    relevant_cols = ica.layout.data.label
    n_samples = length(selected_samples)
    n_channels = length(relevant_cols)
    n_data_samples = nrow(dat.data)

    if !isempty(selected_samples)
        min_idx, max_idx = extrema(selected_samples)
        if min_idx < 1 || max_idx > n_data_samples
            @minimal_error "Selected sample indices out of bounds: [$min_idx, $max_idx] vs data length $n_data_samples"
        end
    end

    dat_matrix = Matrix{Float64}(undef, n_channels, n_samples)
    for (i, col) in enumerate(relevant_cols)
        _copy_col_to_matrix!(dat_matrix, i, dat.data[!, col], selected_samples)
    end

    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica.scale

    components = ica.unmixing * dat_matrix
    return components, size(components, 1)
end
