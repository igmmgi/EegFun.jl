"""
    plot_erp_measurement_gui(erp::ErpData; kwargs...)
    plot_erp_measurement_gui(erps::Vector{ErpData}; kwargs...)

Launch interactive GUI for exploring ERP measurements with live visual feedback.

Perfect for:
- **Teaching**: Show students how measurements are extracted
- **Exploration**: Try different windows/parameters interactively
- **Visual Validation**: Check measurement windows before batch processing

# Arguments
- `erp::ErpData` or `erps::Vector{ErpData}`: Single ERP or multiple conditions to overlay

# Keyword Arguments
- `channel::Union{Symbol,Nothing}`: Initial channel to display (default: first channel)
- `measurement_type::String`: Initial measurement type (default: "mean_amplitude")
- `measurement_window::Union{Tuple{Real,Real},Nothing}`: Initial measurement window (default: auto)
- `baseline_window::Union{Tuple{Real,Real},Nothing}`: Initial baseline window (default: full range)

# Interactive Controls
- **Channel Menu**: Switch between channels
- **Measurement Type Menu**: Select from all 13 measurement types
- **Measurement Window Sliders**: Adjust start/end times (thin blue band around y=0)
- **Baseline Window Sliders**: Adjust baseline interval (thin gray band around y=0)
- **Show Result Markers**: Toggle visualization of measurement results

# Examples
```julia
# Single condition
erp = load_erp("participant_1.jld2")[1]
plot_erp_measurement_gui(erp)

# Multiple conditions overlaid
erps = load_erp("grand_average.jld2")
plot_erp_measurement_gui(erps, channel = :Cz)

# With initial settings for teaching
plot_erp_measurement_gui(erp,
                    measurement_type = "max_peak_latency",
                    measurement_window = (0.3, 0.5),
                    baseline_window = (-0.2, 0.0))
```

# Notes
- For batch measurement extraction, use `erp_measurements()` function
- This GUI is designed for visual exploration and teaching, not batch processing
"""
# Single ErpData - dispatch to vector version
function plot_erp_measurement_gui(
    erp::ErpData;
    channel::Union{Symbol,Nothing} = nothing,
    measurement_type::String = "mean_amplitude",
    measurement_window::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_window::Union{Tuple{Real,Real},Nothing} = nothing,
)
    return plot_erp_measurement_gui([erp]; channel, measurement_type, measurement_window, baseline_window)
end

# Vector of ErpData - main implementation
function plot_erp_measurement_gui(
    erp_vec::Vector{ErpData};
    channel::Union{Symbol,Nothing} = nothing,
    measurement_type::String = "mean_amplitude",
    measurement_window::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_window::Union{Tuple{Real,Real},Nothing} = nothing,
)

    # Get available channels from first ERP
    first_erp = erp_vec[1]
    metadata_cols = meta_labels(first_erp)
    all_channels = setdiff(propertynames(first_erp.data), metadata_cols)

    if isempty(all_channels)
        @minimal_error_throw "No channels found in ERP data"
    end

    # Set initial channel
    initial_channel = isnothing(channel) ? all_channels[1] : channel
    if !(initial_channel in all_channels)
        @minimal_error_throw "Channel $initial_channel not found in data. Available: $(all_channels)"
    end

    # Get time range from first ERP
    time_data = time(first_erp)
    time_min, time_max = extrema(time_data)

    # Set default measurement window if not provided
    if isnothing(measurement_window)
        measurement_window = (time_min, time_max)
    end

    # Set default baseline to full time window if not provided
    if isnothing(baseline_window)
        baseline_window = (time_min, time_max)
    end

    # ===== OBSERVABLES FOR REACTIVE UPDATES =====
    selected_channel = Observable(initial_channel)
    selected_type = Observable(measurement_type)
    meas_window_obs = Observable(measurement_window)
    baseline_window_obs = Observable(baseline_window)
    show_markers_obs = Observable(false)  # Off by default

    # Use plot_erp to create the ERP plot
    # We need to handle baseline reactively, so we'll recreate the plot when baseline changes
    fig = Figure(size = (1200, 700), title = "ERP Measurement Tool")

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
    channel_menu = Menu(controls_grid[3, 1], options = zip(string.(all_channels), all_channels), default = string(initial_channel))

    # Measurement type selection
    Label(controls_grid[4, 1], "Measurement Type:", halign = :left)
    measurement_types = [
        "Mean Amplitude" => "mean_amplitude",
        "Max Peak Amplitude" => "max_peak_amplitude",
        "Min Peak Amplitude" => "min_peak_amplitude",
        "Max Peak Latency" => "max_peak_latency",
        "Min Peak Latency" => "min_peak_latency",
        "Peak-to-Peak Amplitude" => "peak_to_peak_amplitude",
        "Peak-to-Peak Latency" => "peak_to_peak_latency",
        "Rectified Area" => "rectified_area",
        "Integral" => "integral",
        "Positive Area" => "positive_area",
        "Negative Area" => "negative_area",
        "Fractional Area Latency" => "fractional_area_latency",
        "Fractional Peak Latency" => "fractional_peak_latency",
    ]
    # Create menu with display names (labels) that map to full pairs (values)
    type_menu =
        Menu(controls_grid[5, 1], options = zip([p.first for p in measurement_types], measurement_types), default = "Mean Amplitude")

    # Measurement window sliders
    Label(controls_grid[6, 1], "Measurement Window:", halign = :left)
    meas_window_slider = IntervalSlider(controls_grid[7, 1], range = time_min:0.005:time_max, startvalues = measurement_window)
    meas_window_label = Label(controls_grid[8, 1], @sprintf("%.3f s - %.3f s", measurement_window...), halign = :left)

    # Baseline window sliders
    Label(controls_grid[9, 1], "Baseline Window:", halign = :left)
    baseline_window_slider = IntervalSlider(controls_grid[10, 1], range = time_min:0.005:time_max, startvalues = baseline_window)
    baseline_window_label = Label(controls_grid[11, 1], @sprintf("%.3f s - %.3f s", baseline_window...), halign = :left)

    # Show markers toggle
    Label(controls_grid[12, 1], "Show Result Markers:", halign = :left)
    show_markers_toggle = Toggle(controls_grid[13, 1], active = false)  # Off by default

    # Result display
    Label(controls_grid[14, 1], "Result:", fontsize = 14, font = :bold, halign = :left)
    result_label = Label(controls_grid[15, 1], "—", fontsize = 16, halign = :left)

    # Set row gaps
    rowgap!(controls_grid, 10)

    # ===== CREATE ERP PLOT USING plot_erp =====

    # Create axis in plot area
    ax = Axis(plot_area_grid[1, 1], xlabel = "Time (s)", ylabel = "μV", title = "ERP: $(initial_channel)")

    # Define visualization helpers (need ax to be defined first)
    meas_vspan_plots = []
    function update_meas_window!(mw)
        for p in meas_vspan_plots
            delete!(ax, p)
        end
        empty!(meas_vspan_plots)
        # Draw thin band around y=0, ±5% of y-axis range
        lims = ax.finallimits[]
        ymin = lims.origin[2]
        ymax = lims.origin[2] + lims.widths[2]
        yrange = ymax - ymin
        band_height = 0.01 * yrange
        p = poly!(
            ax,
            Point2f[(mw[1], -band_height), (mw[2], -band_height), (mw[2], band_height), (mw[1], band_height)],
            color = (:blue, 0.3),
        )
        push!(meas_vspan_plots, p)
    end

    baseline_vspan_plots = []
    function update_baseline_window!(bw, enabled)
        for p in baseline_vspan_plots
            delete!(ax, p)
        end
        empty!(baseline_vspan_plots)
        if enabled
            # Draw thin band around y=0
            lims = ax.finallimits[]
            ymin = lims.origin[2]
            ymax = lims.origin[2] + lims.widths[2]
            yrange = ymax - ymin
            band_height = 0.01 * yrange
            p = poly!(
                ax,
                Point2f[(bw[1], -band_height), (bw[2], -band_height), (bw[2], band_height), (bw[1], band_height)],
                color = (:gray, 0.2),
            )
            push!(baseline_vspan_plots, p)
        end
    end

    # Result markers visualization (vertical for latencies, horizontal for amplitudes)
    marker_plots = []
    function update_result_markers!(results, show)
        # Clear existing markers
        for p in marker_plots
            delete!(ax, p)
        end
        empty!(marker_plots)

        if !show || isempty(results)
            return
        end

        # Determine if this is a latency or amplitude measurement
        is_latency = selected_type[] in [
            "max_peak_latency",
            "min_peak_latency",
            "peak_to_peak_latency",
            "fractional_area_latency",
            "fractional_peak_latency",
            "onset_latency",
            "offset_latency",
        ]

        # Draw markers for each condition
        colors = length(erp_vec) > 1 ? Makie.cgrad(:jet, length(erp_vec), categorical = true) : [:black]

        for (idx, result) in enumerate(results)
            if haskey(result, :error) || isnothing(result.value) || isnan(result.value)
                continue
            end

            if is_latency
                # Vertical line at latency time point
                p = vlines!(ax, result.value, color = (colors[idx], 0.8), linestyle = :dot, linewidth = 2)
            else
                # Horizontal line at amplitude value
                p = hlines!(ax, result.value, color = (colors[idx], 0.8), linestyle = :dot, linewidth = 2)
            end
            push!(marker_plots, p)
        end
    end

    # Plot ERPs using plot_erp! to add to existing axis
    function update_erp_plot!()
        empty!(ax)  # Clear existing plot

        # Baseline is always enabled
        baseline_int = baseline_window_obs[]

        # Plot using plot_erp!
        plot_erp!(
            fig,
            ax,
            erp_vec,
            channel_selection = channels(selected_channel[]),
            baseline_interval = baseline_int,
            legend = false,  # Disable auto-legend to prevent redrawing
            legend_position = :rt,
        )

        # Redraw vspan overlays (empty! cleared them)
        update_meas_window!(meas_window_obs[])
        update_baseline_window!(baseline_window_obs[], true)  # Always enabled

    end

    # Initial plot
    update_erp_plot!()

    # Add manual legend for multi-condition plots (outside axis, won't be cleared by empty!)
    if length(erp_vec) > 1
        # Get the line elements from the axis
        line_elements = filter(p -> p isa Lines, ax.scene.plots)
        labels = [erp.condition_name for erp in erp_vec]
        Legend(plot_area_grid[1, 2], line_elements[1:length(erp_vec)], labels, halign = :right, valign = :top)
    end

    # Connect menu/slider changes to observables
    on(channel_menu.selection) do ch
        selected_channel[] = ch
        update_erp_plot!()
        ax.title = "$(print_vector([ch])): $(selected_type[])"
        # Redraw markers after plot update
        sleep(0.01)
        update_result_markers!(measurement_results[], show_markers_obs[])
    end

    on(type_menu.selection) do type_pair
        # Menu returns a Pair (display_name => measurement_type), extract the VALUE
        type_str = type_pair isa Pair ? type_pair[2] : type_pair
        selected_type[] = type_str
        ax.title = "$(print_vector([selected_channel[]])): $(type_str)"
    end

    on(meas_window_slider.interval) do interval
        meas_window_obs[] = (interval[1], interval[2])
        meas_window_label.text = @sprintf("%.3f s - %.3f s", interval[1], interval[2])
    end

    on(baseline_window_slider.interval) do interval
        baseline_window_obs[] = (interval[1], interval[2])
        baseline_window_label.text = @sprintf("%.3f s - %.3f s", interval[1], interval[2])
        update_erp_plot!()  # Baseline always enabled, so always update
        # Manually redraw markers after plot clears them
        sleep(0.01)  # Small delay to let measurement_results update
        update_result_markers!(measurement_results[], show_markers_obs[])
    end

    # Initialize visualizations
    update_meas_window!(measurement_window)
    update_baseline_window!(baseline_window, true)

    # Connect to observables
    on(meas_window_obs) do mw
        update_meas_window!(mw)
    end

    on(baseline_window_obs) do bw
        update_baseline_window!(bw, true)  # Always enabled
    end

    # Computed observable for measurement values (one per condition)
    measurement_results = lift(selected_channel, selected_type, meas_window_obs, baseline_window_obs) do ch, type_str, mw, bw
        # Compute for all conditions (baseline always enabled)
        results = []
        for erp in erp_vec
            try
                result = _compute_gui_measurement(erp, ch, type_str, mw, bw)
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
        format_value = function (val)
            if isnothing(val) || isnan(val)
                return "N/A"
            elseif selected_type[] in [
                "max_peak_latency",
                "min_peak_latency",
                "peak_to_peak_latency",
                "fractional_area_latency",
                "fractional_peak_latency",
            ]
                return @sprintf("%.4f s", val)
            elseif selected_type[] in ["rectified_area", "integral", "positive_area", "negative_area"]
                return @sprintf("%.4f μV·s", val)
            else  # Amplitudes
                return @sprintf("%.4f μV", val)
            end
        end

        # Build multi-line text for multiple conditions
        if length(results) == 1
            # Single condition - simple display
            result = results[1]
            if haskey(result, :error)
                result_label.text = "Error: $(result.error)"
            else
                result_label.text = format_value(result.value)
            end
        else
            # Multiple conditions - show each with label
            lines = String[]
            for result in results
                if haskey(result, :error)
                    # Show actual error for debugging
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

    display(fig)
    return (fig = fig)
end


"""
Helper function to compute measurements for GUI.
Returns NamedTuple with :value field (and optionally :error field if failed).
"""
function _compute_gui_measurement(
    erp::ErpData,
    channel::Symbol,
    measurement_type::String,
    measurement_window::Tuple{Real,Real},
    baseline_window::Union{Tuple{Real,Real},Nothing},
)
    # Get data
    time_data = time(erp)
    channel_data = copy(erp.data[!, channel])

    # Apply baseline correction if requested
    if !isnothing(baseline_window)
        time_mask = (time_data .>= baseline_window[1]) .& (time_data .<= baseline_window[2])
        if any(time_mask)
            baseline_mean = mean(channel_data[time_mask])
            channel_data .-= baseline_mean
        end
    end

    # Get measurement window mask
    time_mask = (time_data .>= measurement_window[1]) .& (time_data .<= measurement_window[2])

    if !any(time_mask)
        return (value = NaN,)
    end

    selected_data = channel_data[time_mask]
    selected_times = time_data[time_mask]

    # Compute measurement using existing internal function logic
    local_window = 3  # Fixed neighborhood size for peak detection
    measurement_kwargs = Dict{Symbol,Any}(
        :local_window => local_window,
        :fractional_area_fraction => 0.5,  # Default for fractional area latency
        :fractional_peak_fraction => 0.5,  # Default for fractional peak latency  
        :fractional_peak_direction => :onset,  # Default direction (before peak)
    )

    try
        value = _compute_measurement(selected_data, selected_times, measurement_type, measurement_kwargs, channel)
        if isnothing(value)
            return (value = NaN, error = "Measurement returned nothing")
        end
        return (value = value,)
    catch e
        return (value = NaN, error = string(e))
    end
end
