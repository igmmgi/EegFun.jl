# === 3D Topographic Plot ===

"""
    _load_head_mesh()

Internal lightweight parser for the .obj mesh without requiring FileIO/MeshIO.
Returns a `GeometryBasics.Mesh`.
"""
function _load_head_mesh()
    filename = example_path("models/head.obj")

    vertices = GeometryBasics.Point3f[]
    faces = GeometryBasics.GLTriangleFace[]

    for line in eachline(filename)
        parts = split(line)
        if isempty(parts)
            continue
        end
        if parts[1] == "v"
            push!(vertices, GeometryBasics.Point3f(parse(Float32, parts[2]), parse(Float32, parts[3]), parse(Float32, parts[4])))
        elseif parts[1] == "f"
            # OBJ indices are 1-based, Makie GLTriangleFace is 1-based
            push!(faces, GeometryBasics.GLTriangleFace(parse(Int, parts[2]), parse(Int, parts[3]), parse(Int, parts[4])))
        end
    end

    return GeometryBasics.Mesh(vertices, faces)
end

"""
    _snap_to_mesh(layout_df, head_verts)

Ray-casting algorithm to map spherical electrode coordinates onto the 3D head mesh.
"""
function _snap_to_mesh(layout_df, head_verts)
    # The mathematical center of the Fieldtrip head is roughly at y=-0.015, z=0.015
    center = GeometryBasics.Point3f(0.0, -0.015, 0.015)

    snapped_sensors = GeometryBasics.Point3f[]
    for row in eachrow(layout_df)
        inc_rad = deg2rad(row.inc)
        azi_rad = deg2rad(row.azi)

        # Theoretical direction vector (unit length)
        dir = GeometryBasics.Point3f(sin(inc_rad) * cos(azi_rad), sin(inc_rad) * sin(azi_rad), cos(inc_rad))

        # Find the mesh vertex that aligns best with this direction
        best_radius = 0.100 # fallback
        best_dot = -Inf
        for v in head_verts
            v_vec = v - center
            v_dist = norm(v_vec)
            if v_dist > 0
                v_dir = v_vec / v_dist
                d = dot(v_dir, dir)
                if d > best_dot
                    best_dot = d
                    best_radius = v_dist
                end
            end
        end

        # Snap the sensor to this radius, plus 3mm so it floats visibly on top of the skin
        push!(snapped_sensors, center + dir * (best_radius + 0.003f0))
    end
    return snapped_sensors
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

Plot a perfectly smooth 3D interpolation of EEG data mapped onto a realistically proportioned head mesh.
When `scrubber = true`, opens an interactive time scrubber with slider, real-time wall-clock playback, and amplitude scaling controls.

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
    # Issue a friendly warning if CairoMakie is active
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. 3D plots require GLMakie for full interactivity. You may want to run `using GLMakie; GLMakie.activate!()`."
    end

    # Apply baseline correction if requested
    if !isnothing(baseline_interval)
        dat = baseline(dat, baseline_interval)
    end

    if scrubber
        dat_subset =
            subset(dat, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)
        if nrow(dat_subset.data) > 1
            return _plot_topography_3d_scrubber(dat_subset; display_plot = display_plot, interactive = interactive, kwargs...)
        end
    end

    # Apply subsetting
    dat_subset =
        subset(dat, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)

    # Validate layout has spherical coordinates
    if !hasproperty(dat_subset.layout.data, :inc) || !hasproperty(dat_subset.layout.data, :azi)
        @minimal_error "Cannot create 3D topographic plot: layout missing spherical coordinates (:inc, :azi)."
    end

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_3D_KWARGS, kwargs)

    # Compute channel data averaged over time window
    channel_data = mean.(eachcol(dat_subset.data[!, dat_subset.layout.data.label]))

    # Set up Figure and 3D Axis (LScene)
    _set_window_title(_generate_window_title(dat))
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = Figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    ax = LScene(fig[1, 1], show_axis = false)

    # Add title
    if !isnothing(plot_kwargs[:plot_title]) && !isempty(plot_kwargs[:plot_title])
        if plot_kwargs[:plot_title] != ""
            title_str = plot_kwargs[:plot_title]
        else
            time_min, time_max = extrema(dat_subset.data.time)
            time_unit = get(plot_kwargs, :time_unit, :s)
            title_str = _format_time_range(time_min, time_max, time_unit)
        end
        Label(fig[1, 1, Top()], title_str, fontsize = plot_kwargs[:plot_title_fontsize], font = :bold)
    end

    # Load 3D head mesh and extract vertices
    head_mesh = _load_head_mesh()
    head_verts = GeometryBasics.coordinates(head_mesh)

    # Coregister sensors
    sensors = _snap_to_mesh(dat_subset.layout.data, head_verts)

    # Extract rendering parameters
    colormap = pop!(plot_kwargs, :colormap)
    ylim = pop!(plot_kwargs, :ylim)

    if isnothing(ylim)
        valid_data = filter(!isnan, channel_data)
        if !isempty(valid_data)
            max_abs = maximum(abs.(valid_data))
            ylim = (-max_abs, max_abs)
        else
            ylim = (-1.0, 1.0)
        end
    end

    # Interpolate using ScatteredInterpolation (ThinPlate)
    pts = hcat([[s[1], s[2], s[3]] for s in sensors]...)
    grid = hcat([[v[1], v[2], v[3]] for v in head_verts]...)
    itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), pts, Float64.(channel_data))
    vertex_colors = Float32.(ScatteredInterpolation.evaluate(itp, grid))

    actual_colormap = _resolve_theme_colormap(ax, colormap)
    # Render continuous mesh
    m = mesh!(ax, head_mesh, color = vertex_colors, colormap = actual_colormap, colorrange = ylim, shading = true)

    # Extract point and label kwargs
    point_plot = get(plot_kwargs, :point_plot, true)
    # Scale 2D markersize down for 3D world space coordinates
    point_markersize = get(plot_kwargs, :point_markersize, 12) * 0.0003
    point_color = get(plot_kwargs, :point_color, :black)
    if isnothing(point_color)
        point_color = _resolve_theme_textcolor(ax, point_color)
    end

    label_plot = get(plot_kwargs, :label_plot, true)
    label_fontsize = get(plot_kwargs, :label_fontsize, 20)
    label_color = get(plot_kwargs, :label_color, :black)
    if isnothing(label_color)
        label_color = _resolve_theme_textcolor(ax, label_color)
    end

    # Render Points
    if point_plot
        meshscatter!(ax, sensors, color = point_color, markersize = point_markersize)
    end

    # Render Labels
    if label_plot
        for (i, label) in enumerate(dat_subset.layout.data.label)
            # Offset labels slightly more outward from the head
            lbl_pos = sensors[i] + (sensors[i] - GeometryBasics.Point3f(0.0, -0.015, 0.015)) * 0.05
            text!(ax, lbl_pos, text = string(label), align = (:center, :center), fontsize = label_fontsize, color = label_color)
        end
    end

    # Add Colorbar
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot, true)
    if colorbar_plot
        colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))
        cb_pos = _get_colorbar_position(pop!(plot_kwargs, :colorbar_position, :right), 1:1, 1:1)
        Colorbar(fig[cb_pos...], m; colorbar_kwargs...)
    end

    # Configure initial camera view if requested
    camera_azimuth = pop!(plot_kwargs, :camera_azimuth, nothing)
    camera_elevation = pop!(plot_kwargs, :camera_elevation, nothing)

    if !isnothing(camera_azimuth) || !isnothing(camera_elevation)
        cam = cameracontrols(ax.scene)

        if hasproperty(cam, :eyeposition) && hasproperty(cam, :lookat)
            # Use current radius or default to a reasonable distance
            radius = norm(cam.eyeposition[] - cam.lookat[])
            if radius < 0.01
                radius = 3.0
            end

            # Default to current azimuth/elevation if only one is provided
            az_rad = isnothing(camera_azimuth) ? 0.0 : deg2rad(camera_azimuth)
            el_rad = isnothing(camera_elevation) ? 0.0 : deg2rad(camera_elevation)

            # Convert spherical back to cartesian for camera eye position
            # We use a standard spherical mapping where azimuth=0 looks from +Y (nose)
            x = radius * cos(el_rad) * sin(az_rad)
            y = radius * cos(el_rad) * cos(az_rad)
            z = radius * sin(el_rad)

            update_cam!(ax.scene, GeometryBasics.Point3f(x, y, z), cam.lookat[])
        end
    end

    if get(plot_kwargs, :interactive, true)
        _setup_topo_keyboard_handlers!(fig, ax; zoom_step = get(plot_kwargs, :zoom_step, 0.2))
    end

    if display_plot
        _display_figure(fig)
    end

    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end

"""
    plot_topography_3d(dat::MultiDataFrameEeg, epoch::Int; kwargs...)

Create a 3D topographic plot for a specific epoch of `MultiDataFrameEeg` data.
"""
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
    plot_topography_3d(
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

"""
    _plot_topography_3d_scrubber(dat::SingleDataFrameEeg; display_plot = true, interactive = true, kwargs...)

Internal interactive 3D topography scrubber window with time slider, real-time playback,
and amplitude scaling controls.
"""
function _plot_topography_3d_scrubber(dat::SingleDataFrameEeg; display_plot = true, interactive = true, kwargs...)
    layout = dat.layout
    if !hasproperty(layout.data, :inc) || !hasproperty(layout.data, :azi)
        @minimal_warning "Cannot create 3D topographic plot: layout missing spherical coordinates (:inc, :azi)."
        return
    end

    df = dat.data
    layout_df = layout.data
    matched_indices = Int[]
    available_channels = Symbol[]
    for (i, ch_str) in enumerate(layout_df.label)
        ch_sym = Symbol(ch_str)
        if ch_sym in propertynames(df)
            push!(matched_indices, i)
            push!(available_channels, ch_sym)
        end
    end

    if isempty(available_channels)
        @minimal_warning "No layout-matched channels found in data."
        return
    end

    time_vec = df.time
    n_samples = length(time_vec)
    t_min, t_max = extrema(time_vec)

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_3D_KWARGS, kwargs)
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    colormap = get(plot_kwargs, :colormap, nothing)
    colorbar_plot = get(plot_kwargs, :colorbar_plot, true)
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))

    head_mesh = _load_head_mesh()
    head_verts = GeometryBasics.coordinates(head_mesh)
    layout_subset = layout_df[matched_indices, :]
    sensors = _snap_to_mesh(layout_subset, head_verts)

    pts = hcat([[s[1], s[2], s[3]] for s in sensors]...)
    grid = hcat([[v[1], v[2], v[3]] for v in head_verts]...)

    # Global symmetric amplitude scale across the selected segment
    user_ylim = get(plot_kwargs, :ylim, nothing)
    if !isnothing(user_ylim)
        ylim_obs = Observable{Tuple{Float64,Float64}}((Float64(user_ylim[1]), Float64(user_ylim[2])))
    else
        all_vals = Float64[]
        for ch in available_channels
            append!(all_vals, df[!, ch])
        end
        q99 = Statistics.quantile(abs.(all_vals), 0.99)
        clim = max(round(q99, digits = 1), 1.0)
        ylim_obs = Observable{Tuple{Float64,Float64}}((-clim, clim))
    end

    init_vals = Float64[df[1, ch] for ch in available_channels]
    itp_init = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), pts, init_vals)
    init_colors = Float32.(ScatteredInterpolation.evaluate(itp_init, grid))
    vertex_colors_obs = Observable(init_colors)

    fig = Figure(size = (700, 750); figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    _set_window_title("Topography 3D Scrubber | $(_generate_window_title(dat))")

    current_idx = Observable{Int}(1)
    title_obs = Observable("Topography 3D | Latency: $(@sprintf("%.3f s (%d ms)", time_vec[1], round(Int, time_vec[1] * 1000)))")

    ax = LScene(fig[1, 1], show_axis = false)
    Label(fig[1, 1, Top()], title_obs, fontsize = 16, font = :bold)

    actual_colormap = _resolve_theme_colormap(ax, colormap)
    m = mesh!(ax, head_mesh, color = vertex_colors_obs, colormap = actual_colormap, colorrange = ylim_obs, shading = true)

    point_plot = get(plot_kwargs, :point_plot, true)
    point_markersize = get(plot_kwargs, :point_markersize, 12) * 0.0003
    point_color = get(plot_kwargs, :point_color, :black)
    if isnothing(point_color)
        point_color = _resolve_theme_textcolor(ax, point_color)
    end

    label_plot = get(plot_kwargs, :label_plot, true)
    label_fontsize = get(plot_kwargs, :label_fontsize, 20)
    label_color = get(plot_kwargs, :label_color, :black)
    if isnothing(label_color)
        label_color = _resolve_theme_textcolor(ax, label_color)
    end

    if point_plot
        meshscatter!(ax, sensors, color = point_color, markersize = point_markersize)
    end

    if label_plot
        for (i, label) in enumerate(layout_subset.label)
            lbl_pos = sensors[i] + (sensors[i] - GeometryBasics.Point3f(0.0, -0.015, 0.015)) * 0.05
            text!(ax, lbl_pos, text = string(label), align = (:center, :center), fontsize = label_fontsize, color = label_color)
        end
    end

    if colorbar_plot
        Colorbar(fig[1, 2], m; label = "μV", height = Relative(0.75), colorbar_kwargs...)
    end

    # Controls Panel (Bottom)
    controls = fig[2, 1:2] = GridLayout(tellwidth = false)
    Label(controls[1, 1], "Time:", font = :bold)
    slider = Slider(controls[1, 2], range = 1:n_samples, startvalue = 1)
    time_label = Label(controls[1, 3], lift(idx -> @sprintf("%.3f s", time_vec[idx]), current_idx), font = :bold, width = 75)
    btn_prev = Button(controls[1, 4], label = "<", width = 36)
    btn_play = Button(controls[1, 5], label = "Play", width = 54)
    btn_next = Button(controls[1, 6], label = ">", width = 36)
    btn_mean = Button(controls[1, 7], label = "Mean", width = 54)

    function _update_sample(idx::Int)
        c_idx = clamp(idx, 1, n_samples)
        current_idx[] = c_idx
        t = time_vec[c_idx]
        title_obs[] = @sprintf("Topography 3D | Latency: %.3f s (%d ms)", t, round(Int, t * 1000))
        vals = Float64[df[c_idx, ch] for ch in available_channels]
        itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), pts, vals)
        vertex_colors_obs[] = Float32.(ScatteredInterpolation.evaluate(itp, grid))
    end

    on(slider.value) do idx
        if idx != current_idx[]
            _update_sample(idx)
        end
    end

    on(btn_prev.clicks) do _
        _stop_play()
        set_close_to!(slider, max(1, current_idx[] - 1))
    end

    on(btn_next.clicks) do _
        _stop_play()
        set_close_to!(slider, min(n_samples, current_idx[] + 1))
    end

    is_playing = Observable(false)
    play_task = Ref{Union{Nothing,Task}}(nothing)

    function _start_play()
        is_playing[] = true
        btn_play.label[] = "Pause"
        play_task[] = @async begin
            t_wall_start = time()
            t_data_start = time_vec[current_idx[]]
            if t_data_start >= t_max - 0.002
                t_data_start = t_min
                set_close_to!(slider, 1)
            end
            while is_playing[] && to_value(events(fig).window_open)
                elapsed_wall = time() - t_wall_start
                target_data_time = t_data_start + elapsed_wall
                if target_data_time > t_max
                    t_wall_start = time()
                    t_data_start = t_min
                    target_data_time = t_min
                end
                target_idx = clamp(searchsortedfirst(time_vec, target_data_time), 1, n_samples)
                if target_idx != current_idx[]
                    set_close_to!(slider, target_idx)
                end
                sleep(0.016)
            end
        end
    end

    function _stop_play()
        is_playing[] = false
        btn_play.label[] = "Play"
    end

    on(events(fig).window_open) do is_open
        if !is_open
            _stop_play()
        end
    end

    function _toggle_play()
        if is_playing[]
            _stop_play()
        else
            _start_play()
        end
    end

    on(btn_play.clicks) do _
        _toggle_play()
    end

    on(btn_mean.clicks) do _
        _stop_play()
        mean_vals = Float64[Statistics.mean(df[!, ch]) for ch in available_channels]
        itp_mean = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), pts, mean_vals)
        vertex_colors_obs[] = Float32.(ScatteredInterpolation.evaluate(itp_mean, grid))
        title_obs[] = @sprintf("Topography 3D | Mean over %.3f – %.3f s", t_min, t_max)
    end

    # Initial camera view if requested
    camera_azimuth = get(plot_kwargs, :camera_azimuth, nothing)
    camera_elevation = get(plot_kwargs, :camera_elevation, nothing)

    if !isnothing(camera_azimuth) || !isnothing(camera_elevation)
        cam = cameracontrols(ax.scene)
        if hasproperty(cam, :eyeposition) && hasproperty(cam, :lookat)
            radius = norm(cam.eyeposition[] - cam.lookat[])
            if radius < 0.01
                radius = 3.0
            end
            az_rad = isnothing(camera_azimuth) ? 0.0 : deg2rad(camera_azimuth)
            el_rad = isnothing(camera_elevation) ? 0.0 : deg2rad(camera_elevation)
            x = radius * cos(el_rad) * sin(az_rad)
            y = radius * cos(el_rad) * cos(az_rad)
            z = radius * sin(el_rad)
            update_cam!(ax.scene, GeometryBasics.Point3f(x, y, z), cam.lookat[])
        end
    end

    if interactive
        on(events(fig).keyboardbutton) do event
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if event.key == Keyboard.left
                    _stop_play()
                    set_close_to!(slider, max(1, current_idx[] - 1))
                elseif event.key == Keyboard.right
                    _stop_play()
                    set_close_to!(slider, min(n_samples, current_idx[] + 1))
                elseif event.key == Keyboard.space
                    _toggle_play()
                elseif event.key == Keyboard.up
                    curr_clim = ylim_obs[][2]
                    new_clim = max(1.0, curr_clim * 0.8)
                    ylim_obs[] = (-new_clim, new_clim)
                elseif event.key == Keyboard.down
                    curr_clim = ylim_obs[][2]
                    new_clim = curr_clim * 1.25
                    ylim_obs[] = (-new_clim, new_clim)
                end
            end
        end
    end

    screen = nothing
    if display_plot
        screen = _display_figure(fig)
    end
    _set_window_title("Makie")
    return (fig = fig, axes = [ax], screen = screen)
end
