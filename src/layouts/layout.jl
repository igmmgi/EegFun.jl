"""
    channel_labels(layout::Layout) -> Vector{Symbol}

Return electrode labels from `layout`.
"""
channel_labels(layout::Layout)::Vector{Symbol} = layout.data.label

"""
    positions_polar(layout::Layout) -> DataFrame

Return polar coordinate columns (`:inc`, `:azi`) from `layout`.
"""
positions_polar(layout::Layout)::DataFrame = layout.data[:, [:inc, :azi]]

"""
    positions_2D(layout::Layout) -> DataFrame

Return 2D Cartesian coordinate columns (`:x2`, `:y2`) from `layout`.
"""
positions_2D(layout::Layout)::DataFrame = layout.data[:, [:x2, :y2]]

"""
    positions_3D(layout::Layout) -> DataFrame

Return 3D Cartesian coordinate columns (`:x3`, `:y3`, `:z3`) from `layout`.
"""
positions_3D(layout::Layout)::DataFrame = layout.data[:, [:x3, :y3, :z3]]

"""
    validate_layout(layout::Layout) -> Layout

Validate that `layout` contains required `:label`, `:inc`, and `:azi` columns.
"""
function validate_layout(layout::Layout)
    missing_cols = [c for c in (:label, :inc, :azi) if !hasproperty(layout.data, c)]
    !isempty(missing_cols) && @minimal_error "Layout missing required columns: $missing_cols"
    return layout
end

"""
    has_2d_coords(layout::Layout) -> Bool

Return `true` if `layout` has `:x2` and `:y2` columns.
"""
has_2d_coords(layout::Layout) = hasproperty(layout.data, :x2) && hasproperty(layout.data, :y2)

"""
    has_3d_coords(layout::Layout) -> Bool

Return `true` if `layout` has `:x3`, `:y3`, and `:z3` columns.
"""
has_3d_coords(layout::Layout) = hasproperty(layout.data, :x3) && hasproperty(layout.data, :y3) && hasproperty(layout.data, :z3)

@inline _any_nonzero(v::AbstractVector) = any(!iszero, v)

"""
    has_valid_coordinates(layout::Layout) -> Bool

Return `true` if `layout` has non-zero `:inc` or `:azi` coordinates suitable for topographic plotting.
"""
function has_valid_coordinates(layout::Layout)::Bool
    df = layout.data
    (hasproperty(df, :inc) && hasproperty(df, :azi)) || return false
    return _any_nonzero(df.inc) || _any_nonzero(df.azi)
end

"""
    read_layout(file::String) -> Layout

Read an electrode layout CSV file and return a `Layout` object.
"""
function read_layout(file::String)
    !isfile(file) && @minimal_error "Cannot open file: $file"
    @info "Reading layout from $file"
    df = DataFrame(CSV.File(file; types = Dict(:label => Symbol)))
    rename!(df, Symbol.(names(df)))
    return Layout(df, nothing, nothing, nothing)
end

function _ensure_coordinates_2d!(layout::Layout)
    if !has_2d_coords(layout)
        @info "Converting polar coordinates to 2D Cartesian coordinates"
        polar_to_cartesian_xy!(layout)
    end
    return layout
end

function _ensure_coordinates_3d!(layout::Layout)
    if !has_3d_coords(layout)
        @info "Converting polar coordinates to 3D Cartesian coordinates"
        polar_to_cartesian_xyz!(layout)
    end
    return layout
end

function _check_polar_cols(df::DataFrame)
    (hasproperty(df, :inc) && hasproperty(df, :azi)) || @minimal_error "Layout must contain :inc and :azi columns"
    (eltype(df.inc) <: Number && eltype(df.azi) <: Number) || @minimal_error ":inc and :azi columns must contain numeric values"
end

function _compute_xy(inc_col::AbstractVector, azi_col::AbstractVector, norm_radius::Float64, preserve_radial::Bool)
    n = length(inc_col)
    x2 = Vector{Float64}(undef, n)
    y2 = Vector{Float64}(undef, n)
    rad = Float64(pi / 180)
    if preserve_radial
        inv_90 = 1.0 / 90.0
        @inbounds for i = 1:n
            r = Float64(inc_col[i]) * inv_90
            s, c = sincos(Float64(azi_col[i]) * rad)
            x2[i] = r * c
            y2[i] = r * s
        end
    else
        max_r_sq = 0.0
        @inbounds for i = 1:n
            r = Float64(inc_col[i]) * rad
            s, c = sincos(Float64(azi_col[i]) * rad)
            xi = r * c
            yi = r * s
            x2[i] = xi
            y2[i] = yi
            r_sq = xi * xi + yi * yi
            if !isnan(r_sq) && r_sq > max_r_sq
                max_r_sq = r_sq
            end
        end
        if max_r_sq > 0.0
            scale = norm_radius / sqrt(max_r_sq)
            @inbounds @simd for i = 1:n
                x2[i] *= scale
                y2[i] *= scale
            end
        end
    end
    return x2, y2
end

"""
    polar_to_cartesian_xy!(layout::Layout; normalization_radius::Real = 1.0, preserve_radial_distance::Bool = true)

Convert polar (`:inc`, `:azi`) coordinates into 2D Cartesian (`:x2`, `:y2`) coordinates in place and clear cached neighbors.
When `preserve_radial_distance = true`, `90°` maps to radius `1.0` without shifting the vertex `(0, 0)`.
"""
function polar_to_cartesian_xy!(layout::Layout; normalization_radius::Real = 1.0, preserve_radial_distance::Bool = true)
    df = layout.data
    _check_polar_cols(df)
    x2, y2 = _compute_xy(df.inc, df.azi, Float64(normalization_radius), preserve_radial_distance)
    df[!, :x2] = x2
    df[!, :y2] = y2
    clear_neighbors!(layout)
    return layout
end

polar_to_cartesian_xy!(dat::EegFunData; kwargs...) = (polar_to_cartesian_xy!(dat.layout; kwargs...); dat)

function _compute_xyz(inc_col::AbstractVector, azi_col::AbstractVector)
    n = length(inc_col)
    x3 = Vector{Float64}(undef, n)
    y3 = Vector{Float64}(undef, n)
    z3 = Vector{Float64}(undef, n)
    n == 0 && return x3, y3, z3

    rad = Float64(pi / 180)
    min_x = min_y = min_z = Inf
    max_x = max_y = max_z = -Inf
    @inbounds for i = 1:n
        si, ci = sincos(Float64(inc_col[i]) * rad)
        sa, ca = sincos(Float64(azi_col[i]) * rad)
        xi = si * ca
        yi = si * sa
        zi = ci
        x3[i] = xi
        y3[i] = yi
        z3[i] = zi
        xi < min_x && (min_x = xi)
        xi > max_x && (max_x = xi)
        yi < min_y && (min_y = yi)
        yi > max_y && (max_y = yi)
        zi < min_z && (min_z = zi)
        zi > max_z && (max_z = zi)
    end

    max_range = max(max_x - min_x, max_y - min_y, max_z - min_z)
    if max_range > 0.0
        inv_half = 2.0 / max_range
        cx = 0.5 * (max_x + min_x)
        cy = 0.5 * (max_y + min_y)
        cz = 0.5 * (max_z + min_z)
        @inbounds @simd for i = 1:n
            x3[i] = (x3[i] - cx) * inv_half
            y3[i] = (y3[i] - cy) * inv_half
            z3[i] = (z3[i] - cz) * inv_half
        end
    end
    return x3, y3, z3
end

"""
    polar_to_cartesian_xyz!(layout::Layout)

Convert polar (`:inc`, `:azi`) coordinates into centered, normalized 3D Cartesian (`:x3`, `:y3`, `:z3`) coordinates in place and clear cached neighbors.
"""
function polar_to_cartesian_xyz!(layout::Layout)
    df = layout.data
    _check_polar_cols(df)
    x3, y3, z3 = _compute_xyz(df.inc, df.azi)
    df[!, :x3] = x3
    df[!, :y3] = y3
    df[!, :z3] = z3
    clear_neighbors!(layout)
    return layout
end

polar_to_cartesian_xyz!(dat::EegFunData) = (polar_to_cartesian_xyz!(dat.layout); dat)

@inline _distance_xy(x1::Real, y1::Real, x2::Real, y2::Real)::Float64 = sqrt((x1 - x2)^2 + (y1 - y2)^2)
@inline _squared_distance_xy(x1, y1, x2, y2) = (x1 - x2)^2 + (y1 - y2)^2
@inline _distance_xyz(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)::Float64 = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
@inline _squared_distance_xyz(x1, y1, z1, x2, y2, z2) = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2

function _idw_weights(dists::Vector{Float64})
    k = length(dists)
    weights = Vector{Float64}(undef, k)
    total_inv = 0.0
    @inbounds for idx = 1:k
        w = 1.0 / max(dists[idx], 1e-6)
        weights[idx] = w
        total_inv += w
    end
    inv_sum = 1.0 / total_inv
    @inbounds @simd for idx = 1:k
        weights[idx] *= inv_sum
    end
    return weights
end

function _finalize_neighbors_dict(syms::Vector{Symbol}, nb_ch::Vector{Vector{Symbol}}, nb_d::Vector{Vector{Float64}})
    n = length(syms)
    neighbor_dict = Dict{Symbol,Neighbors}()
    sizehint!(neighbor_dict, n)
    @inbounds for i = 1:n
        chs = nb_ch[i]
        ds = nb_d[i]
        neighbor_dict[syms[i]] = isempty(chs) ? Neighbors(chs, ds, Float64[]) : Neighbors(chs, ds, _idw_weights(ds))
    end
    return neighbor_dict
end

function _find_neighbors_2d(labels::AbstractVector, xs::AbstractVector, ys::AbstractVector, distance_criterion::Real)
    n = length(labels)
    crit_sq = Float64(distance_criterion)^2
    syms = labels isa Vector{Symbol} ? labels : Symbol.(labels)
    nb_ch = [Symbol[] for _ = 1:n]
    nb_d = [Float64[] for _ = 1:n]

    @inbounds for i = 1:n
        xi = Float64(xs[i])
        yi = Float64(ys[i])
        si = syms[i]
        for j = (i+1):n
            dx = xi - Float64(xs[j])
            dy = yi - Float64(ys[j])
            d_sq = dx * dx + dy * dy
            if d_sq <= crit_sq
                d = sqrt(d_sq)
                push!(nb_ch[i], syms[j])
                push!(nb_d[i], d)
                push!(nb_ch[j], si)
                push!(nb_d[j], d)
            end
        end
    end
    return _finalize_neighbors_dict(syms, nb_ch, nb_d)
end

function _find_neighbors_3d(labels::AbstractVector, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector, distance_criterion::Real)
    n = length(labels)
    crit_sq = Float64(distance_criterion)^2
    syms = labels isa Vector{Symbol} ? labels : Symbol.(labels)
    nb_ch = [Symbol[] for _ = 1:n]
    nb_d = [Float64[] for _ = 1:n]

    @inbounds for i = 1:n
        xi = Float64(xs[i])
        yi = Float64(ys[i])
        zi = Float64(zs[i])
        si = syms[i]
        for j = (i+1):n
            dx = xi - Float64(xs[j])
            dy = yi - Float64(ys[j])
            dz = zi - Float64(zs[j])
            d_sq = dx * dx + dy * dy + dz * dz
            if d_sq <= crit_sq
                d = sqrt(d_sq)
                push!(nb_ch[i], syms[j])
                push!(nb_d[i], d)
                push!(nb_ch[j], si)
                push!(nb_d[j], d)
            end
        end
    end
    return _finalize_neighbors_dict(syms, nb_ch, nb_d)
end

function _format_neighbors_toml(neighbors_dict::Dict{Symbol,Neighbors}, nneighbors::Real)
    toml_dict = OrderedDict{String,Any}()
    toml_dict["metadata"] = OrderedDict(
        "average_neighbors_per_electrode" => round(nneighbors; digits = 2),
        "total_electrodes" => length(neighbors_dict),
        "generated_at" => string(now()),
    )

    electrodes = OrderedDict{String,Any}()
    for electrode in sort!(collect(keys(neighbors_dict)))
        nb = neighbors_dict[electrode]
        electrodes[string(electrode)] = OrderedDict(
            "neighbors" => [string(n) for n in nb.channels],
            "distances" => [round(d; digits = 4) for d in nb.distances],
            "weights" => [round(w; digits = 6) for w in nb.weights],
            "neighbor_count" => length(nb.channels),
        )
    end
    toml_dict["electrodes"] = electrodes
    return toml_dict
end

"""
    print_layout_neighbors(layout_or_dict[, filename::String])

Write electrode neighbors, distances, and inverse-distance interpolation weights as TOML to `filename` (or `stdout` if `filename` is omitted).
"""
function print_layout_neighbors(neighbors_dict::Dict{Symbol,Neighbors}, filename::String)
    toml_data = _format_neighbors_toml(neighbors_dict, average_number_of_neighbors(neighbors_dict))
    @info "Printing neighbors to $filename"
    open(filename, "w") do io
        TOML.print(io, toml_data)
    end
end

function print_layout_neighbors(layout::Layout, filename::String)
    isnothing(layout.neighbors) && @minimal_error "No neighbors to print"
    print_layout_neighbors(layout.neighbors, filename)
end

function print_layout_neighbors(layout::Layout)
    isnothing(layout.neighbors) && @minimal_error "No neighbors to print"
    toml_data = _format_neighbors_toml(layout.neighbors, average_number_of_neighbors(layout.neighbors))
    TOML.print(x -> x isa Symbol ? string(x) : x, stdout, toml_data; sorted = true)
end

"""
    average_number_of_neighbors(neighbors_dict::Dict{Symbol,Neighbors}) -> Float64

Calculate the mean number of neighbors per electrode.
"""
function average_number_of_neighbors(neighbors_dict::Dict{Symbol,Neighbors})
    isempty(neighbors_dict) && return 0.0
    return sum(nb -> length(nb.channels), values(neighbors_dict); init = 0) / length(neighbors_dict)
end

"""
    has_neighbors(layout::Layout) -> Bool

Return `true` if neighbor relationships have been calculated on `layout`.
"""
has_neighbors(layout::Layout) = !isnothing(layout.neighbors)

"""
    neighbor_criterion(layout::Layout) -> Union{Nothing,Float64}

Return the distance criterion used for neighbor calculation on `layout`, or `nothing`.
"""
neighbor_criterion(layout::Layout) = layout.criterion

"""
    neighbors(layout::Layout) -> Union{Nothing,Dict{Symbol,Neighbors}}

Return the dictionary of electrode `Neighbors` on `layout`, or `nothing`.
"""
neighbors(layout::Layout) = layout.neighbors

"""
    clear_neighbors!(layout::Layout) -> Layout

Clear cached neighbor data and distance criterion from `layout`.
"""
function clear_neighbors!(layout::Layout)
    layout.neighbors = nothing
    layout.criterion = nothing
    layout.criterion_type = nothing
    return layout
end

function _subset_neighbors(nb::Neighbors, keep_set::Set{Symbol})
    chs = Symbol[]
    ds = Float64[]
    @inbounds for i in eachindex(nb.channels)
        ch = nb.channels[i]
        if ch in keep_set
            push!(chs, ch)
            push!(ds, nb.distances[i])
        end
    end
    isempty(chs) && return Neighbors(chs, ds, Float64[])
    return Neighbors(chs, ds, _idw_weights(ds))
end

function _format_electrode(io, electrode, nb)
    n_neighbors = length(nb.channels)
    if n_neighbors == 0
        println(io, "$(rpad(string(electrode), 6)): 0 neighbors")
    else
        avg_distance = round(mean(nb.distances); digits = 2)
        println(io, "$(rpad(string(electrode), 6)): $(n_neighbors) neighbors (avg dist: $(avg_distance))")
        print(io, "    Neighbors: ")
        for j = 1:n_neighbors
            j > 1 && print(io, ", ")
            print(io, nb.channels[j], " (", round(nb.distances[j]; digits = 2), ",", round(nb.weights[j]; digits = 3), ")")
        end
        println(io)
    end
    println(io)
end

"""
    get_neighbors_xy!(layout::Layout, distance_criterion::Real) -> Layout

Compute and cache 2D electrode neighbors within `distance_criterion` normalized head-radius units.
"""
function get_neighbors_xy!(layout::Layout, distance_criterion::Real)
    if !has_neighbors(layout) || layout.criterion != distance_criterion || layout.criterion_type != :xy
        distance_criterion <= 0 && @minimal_error "Distance criterion must be positive"
        _ensure_coordinates_2d!(layout)
        @info "Calculating 2D neighbors with distance criterion $distance_criterion"
        df = layout.data
        layout.neighbors = _find_neighbors_2d(df.label, df.x2, df.y2, distance_criterion)
        layout.criterion = Float64(distance_criterion)
        layout.criterion_type = :xy
    end
    return layout
end
get_neighbors_xy!(dat::EegFunData, distance_criterion::Real) = (get_neighbors_xy!(dat.layout, distance_criterion); dat)

"""
    get_neighbors_xyz!(layout::Layout, distance_criterion::Real) -> Layout

Compute and cache 3D electrode neighbors within `distance_criterion` normalized coordinate units.
"""
function get_neighbors_xyz!(layout::Layout, distance_criterion::Real)
    if !has_neighbors(layout) || layout.criterion != distance_criterion || layout.criterion_type != :xyz
        distance_criterion <= 0 && @minimal_error "Distance criterion must be positive"
        _ensure_coordinates_3d!(layout)
        @info "Calculating 3D neighbors with distance criterion $distance_criterion"
        df = layout.data
        layout.neighbors = _find_neighbors_3d(df.label, df.x3, df.y3, df.z3, distance_criterion)
        layout.criterion = Float64(distance_criterion)
        layout.criterion_type = :xyz
    end
    return layout
end
get_neighbors_xyz!(dat::EegFunData, distance_criterion::Real) = (get_neighbors_xyz!(dat.layout, distance_criterion); dat)

"""
    rename_channel!(layout::Layout, rename_dict::Dict{Symbol,Symbol}) -> Layout

Rename channels in `layout` in place according to `rename_dict` (`old => new`), clearing cached neighbors.
"""
function rename_channel!(layout::Layout, rename_dict::Dict{Symbol,Symbol})
    isempty(rename_dict) && return layout
    labels = layout.data.label

    matched_rows = Int[]
    target_names = Symbol[]
    @inbounds for i in eachindex(labels)
        lbl = labels[i]
        if haskey(rename_dict, lbl)
            push!(matched_rows, i)
            push!(target_names, rename_dict[lbl])
        end
    end

    if isempty(matched_rows)
        @info "rename!: No channels found to rename"
        return layout
    end

    if length(target_names) != length(unique(target_names))
        duplicate_names = filter(x -> count(==(x), target_names) > 1, unique(target_names))
        @minimal_error "Cannot rename channels to duplicate names: $(join(duplicate_names, ", "))"
    end

    @inbounds for k in eachindex(matched_rows)
        labels[matched_rows[k]] = target_names[k]
    end

    if has_neighbors(layout)
        @info "rename!: Clearing neighbors since channel names have changed"
        clear_neighbors!(layout)
    end

    @info "rename!: Renamed $(length(matched_rows)) channels"
    return layout
end

"""
    rename_channel(layout::Layout, rename_dict::Dict{Symbol,Symbol}) -> Layout

Return a copy of `layout` with channels renamed according to `rename_dict`.
"""
rename_channel(layout::Layout, rename_dict::Dict{Symbol,Symbol}) = rename_channel!(copy(layout), rename_dict)

"""
    get_selected_channels(layout::Layout, channel_selection) -> Vector{Symbol}

Return channel labels in `layout` matching `channel_selection`.
"""
get_selected_channels(layout::Layout, ::Colon) = get_selected_channels(layout, channels())
get_selected_channels(layout::Layout, channel_selection) = get_selected_channels(layout, _to_channel_func(channel_selection))
function get_selected_channels(layout::Layout, channel_selection::Function)
    all_channels = layout.data.label
    return all_channels[channel_selection(all_channels)]
end

"""
    subset_layout!(layout::Layout; channel_selection = channels()) -> Layout

Subset `layout` in place to channels matching `channel_selection`, updating cached neighbors.
"""
function subset_layout!(layout::Layout; channel_selection = channels())
    labels = layout.data.label
    n_channels_before = length(labels)
    selected_channels = get_selected_channels(layout, channel_selection)
    n_channels_after = length(selected_channels)

    if n_channels_after < n_channels_before
        remaining = Set{Symbol}(selected_channels)
        keep_rows = findall(in(remaining), labels)
        layout.data = layout.data[keep_rows, :]

        if has_neighbors(layout)
            old_nb = layout.neighbors
            new_nb = Dict{Symbol,Neighbors}()
            sizehint!(new_nb, n_channels_after)
            empty_nb = Neighbors(Symbol[], Float64[], Float64[])
            for ch in selected_channels
                nb = get(old_nb, ch, nothing)
                new_nb[ch] = isnothing(nb) ? empty_nb : _subset_neighbors(nb, remaining)
            end
            layout.neighbors = new_nb
        end
    end

    @info "subset_layout!: Subset layout from $n_channels_before to $n_channels_after channels"
    return layout
end

subset_layout!(dat::EegFunData; kwargs...) = (subset_layout!(dat.layout; kwargs...); dat)

"""
    subset_layout(layout::Layout; channel_selection = channels()) -> Layout

Return a subsetted copy of `layout` containing channels matching `channel_selection`.
"""
subset_layout(layout::Layout; channel_selection = channels()) = subset_layout!(copy(layout); channel_selection = channel_selection)
