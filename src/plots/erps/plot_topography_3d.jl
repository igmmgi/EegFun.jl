# ==========================================
# 3D Topographic Plot
# ==========================================

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
    plot_topography_3d(dat::SingleDataFrameEeg; kwargs...)

Plot a perfectly smooth 3D interpolation of EEG data mapped onto a realistically proportioned head mesh.

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_3D_KWARGS))
"""
function plot_topography_3d(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = times(),
    display_plot = true,
    kwargs...
)
    # Issue a friendly warning if CairoMakie is active
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. 3D plots require GLMakie for full interactivity. You may want to run `using GLMakie; GLMakie.activate!()`."
    end
    
    # Apply baseline correction if requested
    if !isnothing(baseline_interval)
        dat = baseline(dat, baseline_interval)
    end
    
    # Apply subsetting
    dat_subset = subset(dat, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)
    
    # Validate layout has spherical coordinates
    if !hasproperty(dat_subset.layout.data, :inc) || !hasproperty(dat_subset.layout.data, :azi)
        @minimal_error "Cannot create 3D topographic plot: layout missing spherical coordinates (:inc, :azi)."
    end
    
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_3D_KWARGS, kwargs)
    
    # Compute channel data averaged over time window
    channel_data = mean.(eachcol(dat_subset.data[!, dat_subset.layout.data.label]))
    
    # Set up Figure and 3D Axis (LScene)
    _set_window_title(_generate_window_title(dat))
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])
    fig = Figure(figure_padding = plot_kwargs[:figure_padding])
    ax = LScene(fig[1, 1], show_axis = false)
    
    # Add title
    if plot_kwargs[:show_title]
        if plot_kwargs[:title] != ""
            title_str = plot_kwargs[:title]
        else
            time_min, time_max = extrema(dat_subset.data.time)
            time_unit = get(plot_kwargs, :time_unit, :s)
            title_str = _format_time_range(time_min, time_max, time_unit)
        end
        Label(fig[1, 1, Top()], title_str, fontsize = plot_kwargs[:title_fontsize], font = :bold)
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
    
    # Render continuous mesh
    m = mesh!(ax, head_mesh, color = vertex_colors, colormap = colormap, colorrange = ylim, shading = true)
    
    # Extract point and label kwargs
    point_plot = get(plot_kwargs, :point_plot, true)
    # Scale 2D markersize down for 3D world space coordinates
    point_markersize = get(plot_kwargs, :point_markersize, 12) * 0.0003 
    point_color = get(plot_kwargs, :point_color, :black)
    
    label_plot = get(plot_kwargs, :label_plot, true)
    label_fontsize = get(plot_kwargs, :label_fontsize, 20)
    label_color = get(plot_kwargs, :label_color, :black)
    
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
        colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
        colorbar_position = pop!(plot_kwargs, :colorbar_position, (1, 2))
        Colorbar(fig[colorbar_position...], m; colorbar_kwargs...)
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
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = times(),
    display_plot = true,
    kwargs...
)
    plot_topography_3d(
        epoch_to_continuous(dat, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        baseline_interval = baseline_interval,
        display_plot = display_plot,
        kwargs...
    )
end
