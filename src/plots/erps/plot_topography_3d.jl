const _HEAD_CENTER_3D = GeometryBasics.Point3f(0.0f0, -0.015f0, 0.015f0)

struct _CachedHeadMesh
    mesh::GeometryBasics.Mesh
    verts::Vector{GeometryBasics.Point3f}
    v_dirs::Vector{GeometryBasics.Point3f}
    v_dists::Vector{Float32}
    grid::Matrix{Float64}
end

const _HEAD_MESH_CACHE = Ref{Union{Nothing,_CachedHeadMesh}}(nothing)
const _HEAD_MESH_LOCK = ReentrantLock()

function _parse_head_mesh()
    filename = example_path("models/head.obj")
    vertices = GeometryBasics.Point3f[]
    faces = GeometryBasics.GLTriangleFace[]

    for line in eachline(filename)
        if startswith(line, "v ")
            parts = split(line)
            push!(vertices, GeometryBasics.Point3f(parse(Float32, parts[2]), parse(Float32, parts[3]), parse(Float32, parts[4])))
        elseif startswith(line, "f ")
            parts = split(line)
            push!(faces, GeometryBasics.GLTriangleFace(parse(Int, parts[2]), parse(Int, parts[3]), parse(Int, parts[4])))
        end
    end

    mesh = GeometryBasics.Mesh(vertices, faces)
    verts = GeometryBasics.coordinates(mesh)
    n_verts = length(verts)
    v_dirs = Vector{GeometryBasics.Point3f}(undef, n_verts)
    v_dists = Vector{Float32}(undef, n_verts)
    grid = Matrix{Float64}(undef, 3, n_verts)

    @inbounds for (i, v) in enumerate(verts)
        grid[1, i] = v[1]
        grid[2, i] = v[2]
        grid[3, i] = v[3]
        v_vec = v - _HEAD_CENTER_3D
        d = norm(v_vec)
        v_dists[i] = d
        v_dirs[i] = d > 0 ? (v_vec / d) : GeometryBasics.Point3f(0.0f0, 0.0f0, 0.0f0)
    end

    return _CachedHeadMesh(mesh, verts, v_dirs, v_dists, grid)
end

function _get_cached_head_mesh()
    cached = _HEAD_MESH_CACHE[]
    !isnothing(cached) && return cached
    return lock(_HEAD_MESH_LOCK) do
        isnothing(_HEAD_MESH_CACHE[]) && (_HEAD_MESH_CACHE[] = _parse_head_mesh())
        return _HEAD_MESH_CACHE[]::_CachedHeadMesh
    end
end

function _snap_to_mesh_cached(layout_df::DataFrame, hm::_CachedHeadMesh)
    inc_col = layout_df[!, :inc]
    azi_col = layout_df[!, :azi]
    n_ch = length(inc_col)
    snapped = Vector{GeometryBasics.Point3f}(undef, n_ch)
    v_dirs = hm.v_dirs
    v_dists = hm.v_dists

    @inbounds for i = 1:n_ch
        si, ci = sincos(Float32(deg2rad(inc_col[i])))
        sa, ca = sincos(Float32(deg2rad(azi_col[i])))
        dir = GeometryBasics.Point3f(si * ca, si * sa, ci)

        best_radius = 0.100f0
        best_dot = -Inf32
        for j in eachindex(v_dirs, v_dists)
            vd = v_dists[j]
            if vd > 0
                d = dot(v_dirs[j], dir)
                if d > best_dot
                    best_dot = d
                    best_radius = vd
                end
            end
        end
        snapped[i] = _HEAD_CENTER_3D + dir * (best_radius + 0.003f0)
    end
    return snapped
end

function _sensors_to_matrix(sensors::Vector{GeometryBasics.Point3f})
    pts = Matrix{Float64}(undef, 3, length(sensors))
    @inbounds for (i, s) in enumerate(sensors)
        pts[1, i] = s[1]
        pts[2, i] = s[2]
        pts[3, i] = s[3]
    end
    return pts
end

function _interp_3d_colors(pts::Matrix{Float64}, vals::AbstractVector{Float64}, grid::Matrix{Float64})
    itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), pts, vals)
    eval_res = ScatteredInterpolation.evaluate(itp, grid)
    out = Vector{Float32}(undef, length(eval_res))
    @inbounds @simd for i in eachindex(eval_res)
        out[i] = Float32(eval_res[i])
    end
    return out
end

function _render_3d_sensors_and_labels!(ax::LScene, sensors::Vector{GeometryBasics.Point3f}, labels, plot_kwargs)
    if get(plot_kwargs, :point_plot, true)
        point_markersize = get(plot_kwargs, :point_markersize, 12) * 0.0003
        point_color = get(plot_kwargs, :point_color, :black)
        isnothing(point_color) && (point_color = _resolve_theme_textcolor(ax, point_color))
        meshscatter!(ax, sensors, color = point_color, markersize = point_markersize)
    end

    if get(plot_kwargs, :label_plot, true)
        label_fontsize = get(plot_kwargs, :label_fontsize, 20)
        label_color = get(plot_kwargs, :label_color, :black)
        isnothing(label_color) && (label_color = _resolve_theme_textcolor(ax, label_color))
        for (i, label) in enumerate(labels)
            lbl_pos = sensors[i] + (sensors[i] - _HEAD_CENTER_3D) * 0.05f0
            text!(ax, lbl_pos, text = string(label), align = (:center, :center), fontsize = label_fontsize, color = label_color)
        end
    end
    return nothing
end

function _apply_3d_camera!(ax::LScene, camera_azimuth, camera_elevation)
    (isnothing(camera_azimuth) && isnothing(camera_elevation)) && return nothing
    cam = cameracontrols(ax.scene)
    if hasproperty(cam, :eyeposition) && hasproperty(cam, :lookat)
        radius = norm(cam.eyeposition[] - cam.lookat[])
        radius < 0.01 && (radius = 3.0)
        az_rad = isnothing(camera_azimuth) ? 0.0 : deg2rad(camera_azimuth)
        el_rad = isnothing(camera_elevation) ? 0.0 : deg2rad(camera_elevation)
        eye = GeometryBasics.Point3f(radius * cos(el_rad) * sin(az_rad), radius * cos(el_rad) * cos(az_rad), radius * sin(el_rad))
        update_cam!(ax.scene, eye, cam.lookat[])
    end
    return nothing
end

const PLOT_TOPOGRAPHY_3D_KWARGS = merge(
    PLOT_TOPOGRAPHY_KWARGS,
    Dict{Symbol,Tuple{Any,String}}(
        :camera_azimuth => (nothing, "Initial camera azimuth angle (degrees) for 3D plots. (nothing for auto)"),
        :camera_elevation => (nothing, "Initial camera elevation angle (degrees) for 3D plots. (nothing for auto)"),
    ),
)

"""
    plot_topography_3d(dat::SingleDataFrameEeg; scrubber::Bool = false, kwargs...)
    plot_topography_3d(dat::MultiDataFrameEeg, epoch::Int; scrubber::Bool = false, kwargs...)

Plot a smooth 3D interpolation of EEG data mapped onto a realistically proportioned head mesh.
When `scrubber = true`, opens an interactive time scrubber with slider, wall-clock playback,
and amplitude scaling controls.

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_3D_KWARGS))
"""
function plot_topography_3d(
    dat::SingleDataFrameEeg;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    scrubber::Bool = false,
    display_plot = true,
    interactive = true,
    kwargs...,
)
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. 3D plots require GLMakie for full interactivity. You may want to run `using GLMakie; GLMakie.activate!()`."
    end

    !isnothing(baseline_interval) && (dat = baseline(dat, baseline_interval))
    dat_subset =
        subset(dat, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)

    if scrubber && nrow(dat_subset.data) > 1
        return _plot_topography_3d_scrubber(dat_subset; display_plot = display_plot, interactive = interactive, kwargs...)
    end

    layout_df = dat_subset.layout.data
    (!hasproperty(layout_df, :inc) || !hasproperty(layout_df, :azi)) &&
        @minimal_error "Cannot create 3D topographic plot: layout missing spherical coordinates (:inc, :azi)."

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_3D_KWARGS, kwargs)
    channel_data = Float64[mean(dat_subset.data[!, lbl]) for lbl in layout_df.label]

    _set_window_title(_generate_window_title(dat))
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    ax = LScene(fig[1, 1], show_axis = false)

    show_title = pop!(plot_kwargs, :show_title, true)
    if show_title && !isnothing(plot_kwargs[:plot_title]) && !isempty(plot_kwargs[:plot_title])
        Label(fig[1, 1, Top()], plot_kwargs[:plot_title], fontsize = plot_kwargs[:plot_title_fontsize], font = :bold)
    end

    hm = _get_cached_head_mesh()
    sensors = _snap_to_mesh_cached(layout_df, hm)

    colormap = pop!(plot_kwargs, :colormap)
    colorrange = something(pop!(plot_kwargs, :colorrange), _symmetric_colorrange(channel_data))

    pts = _sensors_to_matrix(sensors)
    vertex_colors = _interp_3d_colors(pts, channel_data, hm.grid)

    m = mesh!(ax, hm.mesh, color = vertex_colors, colormap = _resolve_theme_colormap(ax, colormap), colorrange = colorrange, shading = true)
    _render_3d_sensors_and_labels!(ax, sensors, layout_df.label, plot_kwargs)

    if pop!(plot_kwargs, :colorbar_plot, true)
        colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))
        cb_pos = _get_colorbar_position(pop!(plot_kwargs, :colorbar_position, :right), 1:1, 1:1)
        Colorbar(fig[cb_pos...], m; colorbar_kwargs...)
    end

    _apply_3d_camera!(ax, pop!(plot_kwargs, :camera_azimuth, nothing), pop!(plot_kwargs, :camera_elevation, nothing))

    if get(plot_kwargs, :interactive, true)
        _setup_topo_keyboard_handlers!(fig, ax; zoom_step = get(plot_kwargs, :zoom_step, 0.2))
    end

    display_plot && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax], ax = ax)
end

function plot_topography_3d(
    dat::MultiDataFrameEeg,
    epoch::Int;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    scrubber::Bool = false,
    display_plot = true,
    interactive = true,
    kwargs...,
)
    return plot_topography_3d(
        epoch_to_continuous(dat, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        scrubber = scrubber,
        display_plot = display_plot,
        interactive = interactive,
        kwargs...,
    )
end

function _plot_topography_3d_scrubber(dat::SingleDataFrameEeg; display_plot = true, interactive = true, kwargs...)
    layout_df = dat.layout.data
    if !hasproperty(layout_df, :inc) || !hasproperty(layout_df, :azi)
        @minimal_warning "Cannot create 3D topographic plot: layout missing spherical coordinates (:inc, :azi)."
        return
    end

    df = dat.data
    df_props = Set(propertynames(df))
    matched_indices = Int[]
    available_channels = Symbol[]
    for (i, ch_str) in enumerate(layout_df.label)
        ch_sym = Symbol(ch_str)
        if ch_sym in df_props
            push!(matched_indices, i)
            push!(available_channels, ch_sym)
        end
    end
    if isempty(available_channels)
        @minimal_warning "No layout-matched channels found in data."
        return
    end

    time_vec = df[!, :time]
    n_samples = length(time_vec)
    t_min, t_max = extrema(time_vec)

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_3D_KWARGS, kwargs)
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    colormap = get(plot_kwargs, :colormap, nothing)
    colorbar_plot = get(plot_kwargs, :colorbar_plot, true)
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))

    hm = _get_cached_head_mesh()
    layout_subset = layout_df[matched_indices, :]
    sensors = _snap_to_mesh_cached(layout_subset, hm)
    pts = _sensors_to_matrix(sensors)

    X = Matrix{Float64}(undef, length(available_channels), n_samples)
    for (i, ch) in enumerate(available_channels)
        X[i, :] .= df[!, ch]
    end

    user_colorrange = get(plot_kwargs, :colorrange, nothing)
    ylim_obs = if !isnothing(user_colorrange)
        Observable{Tuple{Float64,Float64}}((Float64(user_colorrange[1]), Float64(user_colorrange[2])))
    else
        q99 = Statistics.quantile!(abs.(vec(copy(X))), 0.99)
        clim = max(round(q99, digits = 1), 1.0)
        Observable{Tuple{Float64,Float64}}((-clim, clim))
    end

    vertex_colors_obs = Observable(_interp_3d_colors(pts, @view(X[:, 1]), hm.grid))

    fig = _eegfun_figure(size = (700, 750); figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    _set_window_title("Topography 3D Scrubber | $(_generate_window_title(dat))")

    current_idx = Observable{Int}(1)
    title_obs = Observable("Topography 3D | Latency: $(@sprintf("%.3f s (%d ms)", time_vec[1], round(Int, time_vec[1] * 1000)))")

    ax = LScene(fig[1, 1], show_axis = false)
    Label(fig[1, 1, Top()], title_obs, fontsize = 16, font = :bold)

    m = mesh!(
        ax,
        hm.mesh,
        color = vertex_colors_obs,
        colormap = _resolve_theme_colormap(ax, colormap),
        colorrange = ylim_obs,
        shading = true,
    )
    _render_3d_sensors_and_labels!(ax, sensors, layout_subset.label, plot_kwargs)
    colorbar_plot && Colorbar(fig[1, 2], m; label = "μV", height = Relative(0.75), colorbar_kwargs...)
    _apply_3d_camera!(ax, get(plot_kwargs, :camera_azimuth, nothing), get(plot_kwargs, :camera_elevation, nothing))

    _setup_scrubber_controls!(
        fig,
        time_vec,
        current_idx,
        ylim_obs,
        interactive,
        function (idx::Int)
            c_idx = clamp(idx, 1, n_samples)
            current_idx[] = c_idx
            t = time_vec[c_idx]
            title_obs[] = @sprintf("Topography 3D | Latency: %.3f s (%d ms)", t, round(Int, t * 1000))
            vertex_colors_obs[] = _interp_3d_colors(pts, @view(X[:, c_idx]), hm.grid)
        end,
        function ()
            mean_vals = vec(mean(X, dims = 2))
            vertex_colors_obs[] = _interp_3d_colors(pts, mean_vals, hm.grid)
            title_obs[] = @sprintf("Topography 3D | Mean over %.3f – %.3f s", t_min, t_max)
        end,
    )

    screen = display_plot ? _display_figure(fig) : nothing
    _set_window_title("Makie")
    return (fig = fig, axes = [ax], ax = ax, screen = screen)
end
