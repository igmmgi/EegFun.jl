"""
    plot_erp_measurement_gui(filepath::String; kwargs...)
    plot_erp_measurement_gui(erp::ErpData; kwargs...)
    plot_erp_measurement_gui(erps::Vector{ErpData}; kwargs...)

Launch interactive GUI for exploring ERP measurements with live visual feedback.

Useful for:
- **Teaching**: Show students how measurements are extracted
- **Exploration**: Try different intervals/parameters interactively
- **Visual Validation**: Check measurement intervals before batch processing

# Arguments
- `filepath::String`: Path to a `.jld2` file containing ERP data
- `erp::ErpData` or `erps::Vector{ErpData}`: Single ERP or multiple conditions to overlay

# Keyword Arguments
- `channel::Union{Symbol,Nothing}`: Initial channel to display (default: first channel)
- `analysis_type::String`: Initial measurement type (default: "mean_amplitude")
- `analysis_interval::Union{Tuple{Real,Real},Nothing}`: Initial measurement interval (default: auto)
- `baseline_interval::Union{Tuple{Real,Real},Nothing}`: Initial baseline interval (default: full range)

# Interactive Controls
- **Channel Menu**: Switch between channels
- **Analysis Type Menu**: Select from all 13 measurement types
- **Analysis Interval Sliders**: Adjust start/end times (thin blue band around y=0)
- **Baseline Interval Sliders**: Adjust baseline interval (thin gray band around y=0)
- **Show Result Markers**: Toggle visualization of measurement results

# Examples
```julia
# Load from file
plot_erp_measurement_gui("participant_1_erps_good.jld2")

# Single condition
erp = read_data("participant_1.jld2")
plot_erp_measurement_gui(erp)

# Multiple conditions overlaid
erps = read_data("grand_average.jld2")
plot_erp_measurement_gui(erps, channel = :Cz)

# With initial settings for teaching
plot_erp_measurement_gui(erp,
                    analysis_type = "max_peak_latency",
                    analysis_interval = (0.3, 0.5),
                    baseline_interval = (-0.2, 0.0))
```

# Notes
- For batch measurement extraction, use `erp_measurements()` function
- This GUI is designed for visual exploration and teaching, not batch processing
"""
# String filepath - load data and dispatch
function plot_erp_measurement_gui(filepath::String; kwargs...)
    data = read_data(filepath)
    if isnothing(data)
        @minimal_error "No data found in file: $filepath"
    end
    return plot_erp_measurement_gui(data; kwargs...)
end
# Single ErpData - dispatch to vector version
function plot_erp_measurement_gui(
    erp::ErpData;
    channel::Union{Symbol,Vector{Symbol},Nothing} = nothing,
    analysis_type::String = "mean_amplitude",
    analysis_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    display_plot::Bool = true,
)
    return plot_erp_measurement_gui([erp]; channel, analysis_type, analysis_interval, baseline_interval, display_plot)
end

# Vector of ErpData - main implementation
function plot_erp_measurement_gui(
    erp_vec::Vector{ErpData};
    channel::Union{Symbol,Vector{Symbol},Nothing} = nothing,
    analysis_type::String = "mean_amplitude",
    analysis_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    display_plot::Bool = true,
)

    # Get available channels from first ERP
    first_erp = erp_vec[1]
    metadata_cols = meta_labels(first_erp)
    all_channels = setdiff(propertynames(first_erp.data), metadata_cols)

    if isempty(all_channels)
        @minimal_error "No channels found in data!"
    end

    # Set initial channel
    initial_channels = isnothing(channel) ? [all_channels[1]] : (channel isa Symbol ? [channel] : channel)
    for ch in initial_channels
        if !(ch in all_channels)
            @minimal_error "Channel $ch not found in data. Available: $(all_channels)"
        end
    end

    # Get time range from first ERP
    time_data = time_vector(first_erp)
    time_min, time_max = extrema(time_data)

    # Set default measurement interval if not provided
    if isnothing(analysis_interval)
        analysis_interval = (time_min, time_max)
    end

    # Set default baseline to full time interval if not provided
    if isnothing(baseline_interval)
        baseline_interval = (time_min, time_max)
    end

    # ===== OBSERVABLES FOR REACTIVE UPDATES =====
    selected_channels = Observable{Vector{Symbol}}(initial_channels)
    selected_type = Observable(analysis_type)
    analysis_interval_obs = Observable(analysis_interval)
    baseline_interval_obs = Observable(baseline_interval)
    show_markers_obs = Observable(false)  # Off by default

    # Use plot_erp to create the ERP plot
    # We need to handle baseline reactively, so we'll recreate the plot when baseline changes
    fig = Figure(size = (1280, 960))
    _set_window_title("ERP Measurement Tool")

    # Create main grid: controls on left (25%), plot area on right (75%)
    main_grid = fig[1, 1] = GridLayout()
    controls_grid = main_grid[1, 1] = GridLayout(valign = :top, tellheight = false)  # Don't constrain row height
    plot_area_grid = main_grid[1, 2] = GridLayout()

    # Set column widths
    colsize!(main_grid, 1, Relative(0.2))
    colsize!(main_grid, 2, Relative(0.8))

    # Let row expand to fill available space
    rowsize!(main_grid, 1, Auto())

    # ===== CONTROLS PANEL =====

    # Title
    Label(controls_grid[1, 1], "Measurement Controls", fontsize = 16, font = :bold, halign = :left)

    # Channel selection
    Label(controls_grid[2, 1], "Channel:", halign = :left)
    channel_menu = Menu(controls_grid[3, 1], options = zip(string.(all_channels), all_channels), default = string(initial_channels[1]))

    # Measurement type selection
    Label(controls_grid[4, 1], "Measurement Type:", halign = :left)
    # Create menu with display names (labels) that map to full pairs (values)
    type_menu = Menu(
        controls_grid[5, 1],
        options = zip([p.first for p in MEASUREMENT_TYPE_LABELS], MEASUREMENT_TYPE_LABELS),
        default = "Mean Amplitude",
    )

    # Measurement interval sliders
    Label(controls_grid[6, 1], "Measurement Interval:", halign = :left)
    meas_window_slider = IntervalSlider(controls_grid[7, 1], range = time_min:0.005:time_max, startvalues = analysis_interval)
    meas_window_label = Label(controls_grid[8, 1], @sprintf("%.3f s - %.3f s", analysis_interval...), halign = :left)

    # Baseline interval sliders
    Label(controls_grid[9, 1], "Baseline Interval:", halign = :left)
    baseline_interval_slider = IntervalSlider(controls_grid[10, 1], range = time_min:0.005:time_max, startvalues = baseline_interval)
    baseline_interval_label = Label(controls_grid[11, 1], @sprintf("%.3f s - %.3f s", baseline_interval...), halign = :left)

    peak_lbl1 = Label(controls_grid[12, 1], "Peak Local Window:", halign = :left)
    local_interval_slider = Slider(controls_grid[13, 1], range = 1:15, startvalue = 3)
    peak_lbl2 = Label(controls_grid[14, 1], lift(v -> "$v samples", local_interval_slider.value), halign = :left)

    frac_lbl1 = Label(controls_grid[15, 1], "Fraction:", halign = :left)
    fraction_slider = Slider(controls_grid[16, 1], range = 0.05:0.05:0.95, startvalue = 0.5)
    frac_lbl2 = Label(controls_grid[17, 1], lift(v -> @sprintf("%.2f", v), fraction_slider.value), halign = :left)

    function update_param_sliders_visibility!(atype)
        is_peak = atype in [
            "max_peak_amplitude",
            "min_peak_amplitude",
            "max_peak_latency",
            "min_peak_latency",
            "peak_to_peak_amplitude",
            "peak_to_peak_latency",
            "fractional_peak_latency",
        ]
        is_frac = atype in ["fractional_area_latency", "fractional_peak_latency"]

        peak_lbl1.blockscene.visible[] = is_peak
        local_interval_slider.blockscene.visible[] = is_peak
        peak_lbl2.blockscene.visible[] = is_peak

        frac_lbl1.blockscene.visible[] = is_frac
        fraction_slider.blockscene.visible[] = is_frac
        frac_lbl2.blockscene.visible[] = is_frac
    end
    update_param_sliders_visibility!("mean_amplitude")

    Label(controls_grid[18, 1], "Show Result Markers:", halign = :left)
    show_markers_toggle = Toggle(controls_grid[19, 1], active = false)

    Label(controls_grid[20, 1], "Result:", fontsize = 14, font = :bold, halign = :right, justification = :right)
    result_label = Label(controls_grid[21, 1], "—", fontsize = 16, halign = :right, justification = :right)

    ax_topo = Axis(controls_grid[22, 1], aspect = DataAspect())
    plot_layout_2d!(fig, ax_topo, first_erp.layout, point_plot = true, label_plot = false, head_linewidth = 1)

    # Hide axis borders
    hidedecorations!(ax_topo)
    hidespines!(ax_topo)
    deregister_interaction!(ax_topo, :rectanglezoom)
    deregister_interaction!(ax_topo, :scrollzoom)

    # Drag state for channel selection
    ch_selection_active = Ref(false)
    ch_selection_start = Ref(Point2f(0, 0))
    ch_selection_rect = poly!(ax_topo, Point2f[(0, 0), (0, 0), (0, 0), (0, 0)], color = (:blue, 0.3), strokecolor = :blue, visible = false)

    # Interactive channel drag / click
    on(events(ax_topo).mousebutton, priority = 50) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                # Check if mouse is actually inside the topography axis
                ax_pos = ax_topo.scene.viewport[]
                mouse_pos = events(fig).mouseposition[]
                is_inside =
                    mouse_pos[1] >= ax_pos.origin[1] &&
                    mouse_pos[1] <= ax_pos.origin[1] + ax_pos.widths[1] &&
                    mouse_pos[2] >= ax_pos.origin[2] &&
                    mouse_pos[2] <= ax_pos.origin[2] + ax_pos.widths[2]

                if !is_inside
                    return Consume(false)
                end

                # Start selection box
                ch_selection_active[] = true
                ch_selection_start[] = mouseposition(ax_topo)
                ch_selection_rect[1] = [ch_selection_start[], ch_selection_start[], ch_selection_start[], ch_selection_start[]]
                ch_selection_rect.visible[] = true
                return Consume(true)
            elseif event.action == Mouse.release && ch_selection_active[]
                # End selection box
                ch_selection_active[] = false
                current_pos = mouseposition(ax_topo)
                ch_selection_rect.visible[] = false

                # compute bounding box and select channels
                x1, x2 = minmax(ch_selection_start[][1], current_pos[1])
                y1, y2 = minmax(ch_selection_start[][2], current_pos[2])

                layout_df = first_erp.layout.data
                selected_this_time = Symbol[]
                for row in eachrow(layout_df)
                    if x1 <= row.x2 <= x2 && y1 <= row.y2 <= y2
                        closest_ch = Symbol(row.label)
                        if closest_ch in all_channels
                            push!(selected_this_time, closest_ch)
                        end
                    end
                end

                # Fallback to nearest point if area was too small (simple click)
                if isempty(selected_this_time) && abs(x2 - x1) < 0.05 && abs(y2 - y1) < 0.05
                    dists = (layout_df.x2 .- current_pos[1]) .^ 2 .+ (layout_df.y2 .- current_pos[2]) .^ 2
                    min_dist, min_idx = findmin(dists)
                    if min_dist < 0.2^2
                        closest_ch = Symbol(layout_df.label[min_idx])
                        if closest_ch in all_channels
                            push!(selected_this_time, closest_ch)
                        end
                    end
                end

                if !isempty(selected_this_time)
                    is_shift_held =
                        Makie.Keyboard.left_shift in events(fig).keyboardstate || Makie.Keyboard.right_shift in events(fig).keyboardstate
                    is_ctrl_held =
                        Makie.Keyboard.left_control in events(fig).keyboardstate ||
                        Makie.Keyboard.right_control in events(fig).keyboardstate

                    if is_shift_held || is_ctrl_held
                        # Toggle behavior
                        new_sel = copy(selected_channels[])
                        for c in selected_this_time
                            if c in new_sel
                                filter!(x -> x != c, new_sel)
                            else
                                push!(new_sel, c)
                            end
                        end

                        if isempty(new_sel)
                            new_sel = [all_channels[1]]
                        end
                        selected_channels[] = new_sel
                    else
                        selected_channels[] = selected_this_time
                    end

                    if length(selected_channels[]) == 1
                        channel_menu.selection[] = string(selected_channels[][1])
                    else
                        channel_menu.selection[] = nothing
                    end
                end

                return Consume(true)
            end
        end
        return Consume(false)
    end

    on(events(ax_topo).mouseposition) do pos
        if ch_selection_active[]
            current_pos = mouseposition(ax_topo)
            start_pos = ch_selection_start[]
            ch_selection_rect[1] = Point2f[
                Point2f(start_pos[1], start_pos[2]),
                Point2f(current_pos[1], start_pos[2]),
                Point2f(current_pos[1], current_pos[2]),
                Point2f(start_pos[1], current_pos[2]),
            ]
        end
    end

    # Highlight selected channels with red dots
    active_pos = lift(selected_channels) do chs
        pts = Point2f[]
        for ch in chs
            idx = findfirst(==(Symbol(ch)), first_erp.layout.data.label)
            if !isnothing(idx)
                push!(pts, Point2f(first_erp.layout.data.x2[idx], first_erp.layout.data.y2[idx]))
            end
        end
        return pts
    end
    scatter!(ax_topo, active_pos, color = :red, markersize = 14)

    # Set row gaps
    rowgap!(controls_grid, 10)

    # ===== CREATE ERP PLOT USING plot_erp =====

    # Create axis in plot area
    ax = Axis(
        plot_area_grid[1, 1],
        xlabel = "Time (s)",
        ylabel = "μV",
        title = "ERP: $(_print_vector(initial_channels))",
        xgridvisible = false,
        ygridvisible = false,
        xminorgridvisible = false,
        yminorgridvisible = false,
    )

    # ===== INTERVAL BAND VISUALIZATION (thin band near y=0 + edge vlines for dragging) =====

    # Observables for band plot objects (recreated on each update_erp_plot!)
    analysis_band_plots = AbstractPlot[]
    baseline_band_plots = AbstractPlot[]

    """Draw a colored interval band (poly + dashed vlines) on the axis."""
    function _draw_interval_band!(plots_list, interval, band_color, edge_color)
        for p in plots_list
            delete!(ax, p)
        end
        empty!(plots_list)
        # Thin band near y=0 (original style)
        lims = ax.finallimits[]
        yrange = lims.widths[2]
        band_height = 0.01 * yrange
        p_band = poly!(
            ax,
            Point2f[(interval[1], -band_height), (interval[2], -band_height), (interval[2], band_height), (interval[1], band_height)],
            color = (band_color, 0.3),
        )
        # Subtle edge lines for drag interaction
        p_left = vlines!(ax, interval[1], color = (edge_color, 0.2), linewidth = 1, linestyle = :dash)
        p_right = vlines!(ax, interval[2], color = (edge_color, 0.2), linewidth = 1, linestyle = :dash)
        push!(plots_list, p_band, p_left, p_right)
    end

    """Redraw the analysis (blue) interval band."""
    update_analysis_band!(mw) = _draw_interval_band!(analysis_band_plots, mw, :blue, :blue)
    """Redraw the baseline (gray) interval band."""
    update_baseline_band!(bw) = _draw_interval_band!(baseline_band_plots, bw, :gray, :gray)

    # Result markers visualization (vertical for latencies, horizontal for amplitudes)
    marker_plots = AbstractPlot[]
    """Add / clear marker lines (latency → vlines, amplitude → hlines) on the axis."""
    function update_result_markers!(results, show)
        # Clear existing markers
        for p in marker_plots
            delete!(ax, p)
        end
        empty!(marker_plots)

        if !show || isempty(results)
            return
        end

        colors = length(erp_vec) > 1 ? Makie.cgrad(:jet, length(erp_vec), categorical = true) : [:black]

        for (idx, result) in enumerate(results)
            if haskey(result, :error) || !haskey(result, :value) || isnothing(result.value) || isnan(result.value)
                continue
            end

            c = (colors[idx], 0.8)

            if haskey(result, :point)
                if _is_latency_measurement(selected_type[])
                    p = vlines!(ax, result.point[1]; color = c, linewidth = 2)
                else
                    p = hlines!(ax, result.point[2]; color = c, linewidth = 2)
                end
                push!(marker_plots, p)
            elseif haskey(result, :points)
                p = lines!(ax, [Point2f(pt) for pt in result.points]; color = c, linestyle = :dash, linewidth = 2)
                push!(marker_plots, p)
            elseif haskey(result, :span)
                p = lines!(ax, [Point2f(result.span[1], result.value), Point2f(result.span[2], result.value)]; color = c, linewidth = 4)
                push!(marker_plots, p)
            elseif haskey(result, :curve)
                y_band = result.curve
                y_baseline = zeros(length(y_band))
                y_upper = y_band
                y_lower = y_baseline
                if selected_type[] == "positive_area"
                    y_lower = y_baseline
                    y_upper = max.(y_band, 0.0)
                elseif selected_type[] == "negative_area"
                    y_upper = y_baseline
                    y_lower = min.(y_band, 0.0)
                end

                p = band!(ax, result.x, y_lower, y_upper; color = (colors[idx], 0.3))
                push!(marker_plots, p)
            else
                if _is_latency_measurement(selected_type[])
                    p = vlines!(ax, result.value, color = c, linestyle = :dot, linewidth = 2)
                else
                    p = hlines!(ax, result.value, color = c, linestyle = :dot, linewidth = 2)
                end
                push!(marker_plots, p)
            end
        end
    end

    # State for stable y-limits
    cached_ylims = Ref{Tuple{Float64,Float64}}((0.0, 1.0))
    cached_ch = Ref{Vector{Symbol}}(Symbol[])

    # Plot ERPs using plot_erp! to add to existing axis
    function update_erp_plot!()
        empty!(ax)  # Clear existing plot

        chs = selected_channels[]
        bi = baseline_interval_obs[]

        # Compute mathematically guaranteed limits the first time a channel is selected
        if chs != cached_ch[]
            max_range = 0.0
            for erp in erp_vec
                channels_df = erp.data[!, chs]
                raw_data = length(chs) == 1 ? channels_df[!, 1] : vec(sum(Matrix(channels_df), dims = 2) ./ length(chs))
                r_min = minimum(raw_data)
                r_max = maximum(raw_data)
                max_range = max(max_range, r_max - r_min)
            end

            # The baselined signal will mathematically ALWAYS be within ±(max_raw - min_raw)
            # no matter what baseline interval is dragged! Add a 5% margin.
            if max_range == 0.0
                max_range = 1.0
            end

            cached_ylims[] = (-1.05 * max_range, 1.05 * max_range)
            cached_ch[] = chs
        end

        if length(chs) > 1
            lbl_sym = Symbol("ROI_$(length(chs))_channels")
            plot_vec = channel_average(erp_vec, channel_selections = [channels(chs)], output_labels = [lbl_sym], reduce = true)
            plot_ch_selection = channels([lbl_sym])
        else
            plot_vec = erp_vec
            plot_ch_selection = channels(chs)
        end

        # Plot using plot_erp!
        plot_erp!(
            fig,
            ax,
            plot_vec,
            channel_selection = plot_ch_selection,
            baseline_interval = bi,
            legend = false,  # Disable auto-legend to prevent redrawing
            legend_position = :rt,
        )

        # Apply stable limits
        ylims!(ax, cached_ylims[][1], cached_ylims[][2])

        # Redraw band overlays (empty! cleared them)
        update_analysis_band!(analysis_interval_obs[])
        update_baseline_band!(baseline_interval_obs[])
        _set_window_title("ERP Measurement Tool")

    end

    # Initial plot
    update_erp_plot!()

    # TODO: interactive legend not working here!?
    # Add manual legend for multi-condition plots (outside axis, won't be cleared by empty!)
    if length(erp_vec) > 1
        # Get the line elements from the axis
        line_elements = filter(p -> p isa Lines, ax.scene.plots)
        labels = [erp.condition_name for erp in erp_vec]
        Legend(plot_area_grid[1, 2], line_elements[1:length(erp_vec)], labels, halign = :right, valign = :top)
    end

    # ===== DRAG STATE FOR INTERACTIVE EDGE/BAND DRAGGING =====
    # Targets: :none, :analysis_left, :analysis_right, :analysis_middle, :baseline_left, :baseline_right, :baseline_middle
    drag_target = Observable{Symbol}(:none)
    drag_suppress_slider = Ref(false)  # prevent feedback loops during drag
    drag_offset = Ref(0.0)  # offset from mouse to interval start (for middle drag)

    # Hit detection: find which drag target (if any) is at the mouse x position
    """Hit-test the mouse x position against interval edges and bands."""
    function _find_drag_target(mouse_x)
        time_range = time_max - time_min
        tolerance = 0.02 * time_range  # 2% of time range

        ai = analysis_interval_obs[]
        bi = baseline_interval_obs[]

        # Edge detection first (higher priority than middle)
        if abs(mouse_x - ai[1]) < tolerance
            return :analysis_left
        elseif abs(mouse_x - ai[2]) < tolerance
            return :analysis_right
        elseif abs(mouse_x - bi[1]) < tolerance
            return :baseline_left
        elseif abs(mouse_x - bi[2]) < tolerance
            return :baseline_right
        end

        # Middle detection: click inside the band (analysis has priority)
        if ai[1] < mouse_x < ai[2]
            return :analysis_middle
        elseif bi[1] < mouse_x < bi[2]
            return :baseline_middle
        end

        return :none
    end

    # Mouse press: start dragging if on an edge or inside a band
    on(events(ax).mousebutton, priority = 100) do event
        if event.button == Mouse.left && event.action == Mouse.press
            # Check if mouse is over the axis
            mouse_pos = events(fig).mouseposition[]
            area = ax.scene.viewport[]
            if !(
                area.origin[1] <= mouse_pos[1] <= area.origin[1] + area.widths[1] &&
                area.origin[2] <= mouse_pos[2] <= area.origin[2] + area.widths[2]
            )
                return Consume(false)
            end

            world_x = mouseposition(ax)[1]
            target = _find_drag_target(world_x)
            if target != :none
                drag_target[] = target
                # For middle drag, store offset from mouse to interval left edge
                if target == :analysis_middle
                    drag_offset[] = world_x - analysis_interval_obs[][1]
                elseif target == :baseline_middle
                    drag_offset[] = world_x - baseline_interval_obs[][1]
                end
                return Consume(true)  # consume to prevent axis pan
            end
        elseif event.button == Mouse.left && event.action == Mouse.release
            if drag_target[] != :none
                drag_target[] = :none
                return Consume(true)
            end
        end
        return Consume(false)
    end

    # Mouse move: update position during drag
    on(events(fig).mouseposition) do _
        if drag_target[] == :none
            return
        end

        world_x = mouseposition(ax)[1]
        world_x = clamp(world_x, time_min, time_max)

        target = drag_target[]
        drag_suppress_slider[] = true

        if target == :analysis_left
            ai = analysis_interval_obs[]
            new_left = min(world_x, ai[2] - 0.005)
            analysis_interval_obs[] = (new_left, ai[2])
            set_close_to!(meas_window_slider, new_left, ai[2])
            meas_window_label.text = @sprintf("%.3f s - %.3f s", new_left, ai[2])
        elseif target == :analysis_right
            ai = analysis_interval_obs[]
            new_right = max(world_x, ai[1] + 0.005)
            analysis_interval_obs[] = (ai[1], new_right)
            set_close_to!(meas_window_slider, ai[1], new_right)
            meas_window_label.text = @sprintf("%.3f s - %.3f s", ai[1], new_right)
        elseif target == :analysis_middle
            ai = analysis_interval_obs[]
            width = ai[2] - ai[1]
            new_left = world_x - drag_offset[]
            new_left = clamp(new_left, time_min, time_max - width)
            new_right = new_left + width
            analysis_interval_obs[] = (new_left, new_right)
            set_close_to!(meas_window_slider, new_left, new_right)
            meas_window_label.text = @sprintf("%.3f s - %.3f s", new_left, new_right)
        elseif target == :baseline_left
            bi = baseline_interval_obs[]
            new_left = min(world_x, bi[2] - 0.005)
            baseline_interval_obs[] = (new_left, bi[2])
            set_close_to!(baseline_interval_slider, new_left, bi[2])
            baseline_interval_label.text = @sprintf("%.3f s - %.3f s", new_left, bi[2])
            update_erp_plot!()
        elseif target == :baseline_right
            bi = baseline_interval_obs[]
            new_right = max(world_x, bi[1] + 0.005)
            baseline_interval_obs[] = (bi[1], new_right)
            set_close_to!(baseline_interval_slider, bi[1], new_right)
            baseline_interval_label.text = @sprintf("%.3f s - %.3f s", bi[1], new_right)
            update_erp_plot!()
        elseif target == :baseline_middle
            bi = baseline_interval_obs[]
            width = bi[2] - bi[1]
            new_left = world_x - drag_offset[]
            new_left = clamp(new_left, time_min, time_max - width)
            new_right = new_left + width
            baseline_interval_obs[] = (new_left, new_right)
            set_close_to!(baseline_interval_slider, new_left, new_right)
            baseline_interval_label.text = @sprintf("%.3f s - %.3f s", new_left, new_right)
            update_erp_plot!()
        end

        drag_suppress_slider[] = false
    end

    # Connect menu/slider changes to observables
    on(channel_menu.selection) do ch
        if !isnothing(ch)
            sym_ch = Symbol(ch)
            if length(selected_channels[]) > 1 || selected_channels[][1] != sym_ch
                selected_channels[] = [sym_ch]
            end
        end
    end

    # Plot update triggered by channel selection change
    on(selected_channels) do chs
        update_erp_plot!()
        ax.title = "$(_print_vector(chs)): $(selected_type[])"
        # Redraw markers after plot update
        yield()
        update_result_markers!(measurement_results[], show_markers_obs[])
    end

    on(type_menu.selection) do type_pair
        # Menu returns a Pair (display_name => analysis_type), extract the VALUE
        type_str = type_pair isa Pair ? type_pair[2] : type_pair
        selected_type[] = type_str
        ax.title = "$(_print_vector(selected_channels[])): $(type_str)"
        update_param_sliders_visibility!(type_str)
    end

    on(meas_window_slider.interval) do interval
        drag_suppress_slider[] && return  # skip if drag is driving this
        analysis_interval_obs[] = (interval[1], interval[2])
        meas_window_label.text = @sprintf("%.3f s - %.3f s", interval[1], interval[2])
    end

    on(baseline_interval_slider.interval) do interval
        drag_suppress_slider[] && return  # skip if drag is driving this
        baseline_interval_obs[] = (interval[1], interval[2])
        baseline_interval_label.text = @sprintf("%.3f s - %.3f s", interval[1], interval[2])
        update_erp_plot!()  # Baseline always enabled, so always update
        # Manually redraw markers after plot clears them
        yield()  # Small delay to let measurement_results update
        update_result_markers!(measurement_results[], show_markers_obs[])
    end

    # Connect observables to band visuals
    on(analysis_interval_obs) do mw
        update_analysis_band!(mw)
    end

    on(baseline_interval_obs) do bw
        update_baseline_band!(bw)
    end

    # Computed observable for measurement values (one per condition)
    measurement_results = lift(
        selected_channels,
        selected_type,
        analysis_interval_obs,
        baseline_interval_obs,
        local_interval_slider.value,
        fraction_slider.value,
    ) do chs, type_str, mw, bw, li, frac
        # Compute for all conditions (baseline always enabled)
        results = NamedTuple[]

        measurement_kwargs = Dict{Symbol,Any}(k => v[1] for (k, v) in ERP_MEASUREMENTS_KWARGS)
        measurement_kwargs[:local_interval] = li
        measurement_kwargs[:fractional_area_fraction] = frac
        measurement_kwargs[:fractional_peak_fraction] = frac

        for erp in erp_vec
            try
                result = _compute_gui_measurement(erp, chs, type_str, mw, bw, measurement_kwargs)
                push!(results, (condition = erp.condition_name, result...))
            catch e
                push!(results, (condition = erp.condition_name, value = NaN, error = string(e)))
            end
        end
        return results
    end

    # Update result label to show all condition values
    on(measurement_results) do results
        if isempty(results)
            result_label.text = "—"
            return
        end

        # Format based on measurement type
        function format_value(val)
            if isnothing(val) || isnan(val)
                return "   N/A  "
            elseif _is_latency_measurement(selected_type[])
                return @sprintf("% .3f s", val)
            elseif selected_type[] in ["rectified_area", "integral", "positive_area", "negative_area"]
                return @sprintf("% .3f μV·s", val)
            else  # Amplitudes
                return @sprintf("% .3f μV", val)
            end
        end

        if length(results) == 1
            result = results[1]
            if haskey(result, :error)
                result_label.text = "Error: $(result.error)"
            else
                result_label.text = format_value(result.value)
            end
        else
            lines = String[]
            for result in results
                if haskey(result, :error)
                    push!(lines, "$(result.condition): $(result.error)")
                else
                    push!(lines, "$(result.condition): $(format_value(result.value))")
                end
            end
            result_label.text = join(lines, "\n")
        end

        # Update result markers
        update_result_markers!(results, show_markers_obs[])
    end

    # Connect toggle to observable
    on(show_markers_toggle.active) do enabled
        show_markers_obs[] = enabled
        update_result_markers!(measurement_results[], enabled)
    end

    # Trigger initial result display 
    notify(measurement_results)

    if display_plot
        _display_figure(fig)
    end
    return (fig = fig)

end


"""
Helper function to compute measurements for GUI.
Returns NamedTuple with :value field and geometric objects for marker plotting.
"""
function _compute_gui_measurement(
    erp::ErpData,
    channels::Vector{Symbol},
    analysis_type,
    analysis_interval,
    baseline_interval,
    measurement_kwargs,
)
    # Apply baseline correction using existing infrastructure
    erp_data = baseline(erp, baseline_interval)
    time_data = time_vector(erp_data)

    channels_df = erp_data.data[!, channels]
    channel_data = length(channels) == 1 ? channels_df[!, 1] : vec(sum(Matrix(channels_df), dims = 2) ./ length(channels))

    time_mask = (time_data .>= analysis_interval[1]) .& (time_data .<= analysis_interval[2])

    if !any(time_mask)
        return (value = NaN,)
    end

    selected_data = channel_data[time_mask]
    selected_times = time_data[time_mask]

    dummy_channel_sym = length(channels) == 1 ? channels[1] : Symbol("ROI_$(length(channels))_channels")
    value = _compute_measurement(selected_data, selected_times, analysis_type, measurement_kwargs, dummy_channel_sym)
    if isnothing(value)
        return (value = NaN,)
    end

    result = Dict{Symbol,Any}(:value => value)

    if analysis_type == "mean_amplitude"
        result[:span] = (analysis_interval[1], analysis_interval[2])
    elseif analysis_type in ["max_peak_amplitude", "min_peak_amplitude", "max_peak_latency", "min_peak_latency"]
        peak_type = startswith(analysis_type, "max") ? :max : :min
        local_interval = measurement_kwargs[:local_interval]
        peak_val, peak_idx = _compute_peak_measurement(selected_data, peak_type, local_interval, dummy_channel_sym)
        if !isnothing(peak_val) && !isnothing(peak_idx)
            result[:point] = (selected_times[peak_idx], peak_val)
        end
    elseif analysis_type in ["peak_to_peak_amplitude", "peak_to_peak_latency"]
        local_interval = measurement_kwargs[:local_interval]
        max_val, max_idx = _compute_peak_measurement(selected_data, :max, local_interval, dummy_channel_sym)
        min_val, min_idx = _compute_peak_measurement(selected_data, :min, local_interval, dummy_channel_sym)
        if !isnothing(max_val) && !isnothing(min_val)
            result[:points] = [(selected_times[max_idx], max_val), (selected_times[min_idx], min_val)]
        end
    elseif analysis_type in ["rectified_area", "integral", "positive_area", "negative_area"]
        result[:curve] = selected_data
        result[:x] = selected_times
    elseif analysis_type in ["fractional_area_latency", "fractional_peak_latency"]
        idx = argmin(abs.(time_data .- value))
        result[:point] = (value, channel_data[idx])
    end

    return (; result...)
end
