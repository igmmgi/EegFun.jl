"""
    plot_erp_filter_gui(erp::ErpData; kwargs...)
    plot_erp_filter_gui(erps::Vector{ErpData}; kwargs...)

Launch interactive GUI for exploring filter effects on ERP waveforms.

Perfect for:
- **Teaching**: Show students how filters affect ERP signals
- **Exploration**: Try different filter parameters interactively
- **Visual Validation**: Compare original vs filtered waveforms before batch processing

# Arguments
- `erp::ErpData` or `erps::Vector{ErpData}`: Single ERP or multiple conditions

# Keyword Arguments
- `channel::Union{Symbol,Nothing}`: Initial channel to display (default: first channel)

# Interactive Controls
- **Channel Menu**: Switch between channels
- **Filter Method Menu**: butterworth (IIR), fir (Hamming)
- **Filter Function Menu**: filtfilt (zero-phase), filt (one-pass)
- **Lowpass Filter**: Toggle, cutoff, order sliders
- **Highpass Filter**: Toggle, cutoff, order sliders

# Notes
- Both filters disabled by default
- For batch filtering, use `lowpass_filter` or `highpass_filter` functions
"""
# Single ErpData - dispatch to vector version
function plot_erp_filter_gui(erp::ErpData; channel::Union{Symbol,Nothing} = nothing)
    return plot_erp_filter_gui([erp]; channel)
end

# Vector of ErpData - main implementation
function plot_erp_filter_gui(erp_vec::Vector{ErpData}; channel::Union{Symbol,Nothing} = nothing)

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

    # Get sample rate
    sample_rate = first_erp.sample_rate
    nyquist = sample_rate / 2

    # ===== OBSERVABLES FOR REACTIVE UPDATES =====
    selected_channel = Observable(initial_channel)

    # Lowpass controls (disabled by default)
    lowpass_enabled = Observable(false)
    lowpass_cutoff = Observable(30.0)
    lowpass_order = Observable(2)

    # Highpass controls (disabled by default)
    highpass_enabled = Observable(false)
    highpass_cutoff = Observable(0.1)
    highpass_order = Observable(1)

    # Common controls
    filter_method_obs = Observable("butterworth")
    filter_func_obs = Observable("filtfilt")
    mirror_enabled = Observable(false)

    # Create Figure
    fig = Figure(size = (1200, 700), title = "ERP Filter Tool")

    # Create main grid: controls on left (20%), plot area on right (80%)
    main_grid = fig[1, 1] = GridLayout()
    controls_grid = main_grid[1, 1] = GridLayout(valign = :top, tellheight = false)
    plot_area_grid = main_grid[1, 2] = GridLayout()

    # Set column widths
    colsize!(main_grid, 1, Relative(0.2))
    colsize!(main_grid, 2, Relative(0.8))
    rowsize!(main_grid, 1, Auto())

    # ===== CONTROLS PANEL =====
    row = 1

    Label(controls_grid[row, 1], "Filter Controls", fontsize = 16, font = :bold, halign = :left)
    row += 1

    # Channel selection
    Label(controls_grid[row, 1], "Channel:", halign = :left)
    row += 1
    channel_menu = Menu(controls_grid[row, 1], options = zip(string.(all_channels), all_channels), default = string(initial_channel))
    row += 1

    # Filter method selection
    Label(controls_grid[row, 1], "Filter Method:", halign = :left)
    row += 1
    filter_methods = ["Butterworth (IIR)" => "butterworth", "FIR (Hamming)" => "fir"]
    method_menu =
        Menu(controls_grid[row, 1], options = zip([p.first for p in filter_methods], filter_methods), default = "Butterworth (IIR)")
    row += 1

    # Filter function selection
    Label(controls_grid[row, 1], "Filter Function:", halign = :left)
    row += 1
    filter_funcs = ["filtfilt" => "filtfilt", "filt" => "filt"]
    func_menu = Menu(controls_grid[row, 1], options = zip([p.first for p in filter_funcs], filter_funcs), default = "filtfilt")
    row += 1

    # Mirror option
    Label(controls_grid[row, 1], "Mirror (for edge artifacts):", halign = :left)
    row += 1

    hbox_mirror = GridLayout(controls_grid[row, 1])
    mirror_toggle = Toggle(hbox_mirror[1, 1], active = false, halign = :left)
    Label(hbox_mirror[1, 2], "Enable", halign = :left)
    colsize!(hbox_mirror, 1, Auto())
    colsize!(hbox_mirror, 2, Auto())
    row += 1

    # === LOWPASS SECTION ===
    Label(controls_grid[row, 1], "Lowpass Filter:", fontsize = 14, font = :bold, halign = :left)
    row += 1

    hbox = GridLayout(controls_grid[row, 1])
    lowpass_toggle = Toggle(hbox[1, 1], active = false, halign = :left)
    Label(hbox[1, 2], "Enable", halign = :left)
    colsize!(hbox, 1, Auto())
    colsize!(hbox, 2, Auto())
    row += 1

    Label(controls_grid[row, 1], "Cutoff Frequency:", halign = :left)
    row += 1
    lowpass_cutoff_slider = Slider(controls_grid[row, 1], range = 1.0:1.0:(nyquist*0.9), startvalue = 30.0)
    row += 1
    lowpass_cutoff_label = Label(controls_grid[row, 1], "30.0 Hz", halign = :left)
    row += 1

    Label(controls_grid[row, 1], "Order:", halign = :left)
    row += 1
    lowpass_order_slider = Slider(controls_grid[row, 1], range = 1:10, startvalue = 2)
    row += 1
    lowpass_order_label = Label(controls_grid[row, 1], "2", halign = :left)
    row += 1

    # === HIGHPASS SECTION ===
    Label(controls_grid[row, 1], "Highpass Filter:", fontsize = 14, font = :bold, halign = :left)
    row += 1

    hbox2 = GridLayout(controls_grid[row, 1])
    highpass_toggle = Toggle(hbox2[1, 1], active = false, halign = :left)
    Label(hbox2[1, 2], "Enable", halign = :left)
    colsize!(hbox2, 1, Auto())
    colsize!(hbox2, 2, Auto())
    row += 1

    Label(controls_grid[row, 1], "Cutoff Frequency:", halign = :left)
    row += 1
    highpass_cutoff_slider = Slider(controls_grid[row, 1], range = 0.1:0.1:10.0, startvalue = 0.1)
    row += 1
    highpass_cutoff_label = Label(controls_grid[row, 1], "0.1 Hz", halign = :left)
    row += 1

    Label(controls_grid[row, 1], "Order:", halign = :left)
    row += 1
    highpass_order_slider = Slider(controls_grid[row, 1], range = 1:10, startvalue = 1)
    row += 1
    highpass_order_label = Label(controls_grid[row, 1], "1", halign = :left)
    row += 1

    # Set row gaps
    rowgap!(controls_grid, 8)

    # ===== CREATE SUBPLOTS FOR EACH CONDITION =====
    n_conditions = length(erp_vec)
    nrows, ncols = best_rect(n_conditions)

    axes = []
    for i = 1:n_conditions
        row_idx = div(i - 1, ncols) + 1
        col_idx = mod(i - 1, ncols) + 1
        ax = Axis(
            plot_area_grid[row_idx, col_idx],
            xlabel = row_idx == nrows ? "Time (s)" : "",
            ylabel = col_idx == 1 ? "μV" : "",
            title = erp_vec[i].condition_name,
        )
        push!(axes, ax)
    end

    # Link all x-axes
    if length(axes) > 1
        for i = 2:lastindex(axes)
            linkxaxes!(axes[1], axes[i])
        end
    end

    # Function to update all subplots
    function update_plot!()
        for (idx, ax) in enumerate(axes)
            empty!(ax)

            # Plot original ERP in black
            plot_erp!(fig, ax, [erp_vec[idx]], channel_selection = channels(selected_channel[]), legend = false)

            # Apply filters if enabled
            if lowpass_enabled[] || highpass_enabled[]
                filtered_erp = copy(erp_vec[idx])

                # Mirror if enabled
                if mirror_enabled[]
                    mirror!(filtered_erp, :both)
                end

                filt_method = filter_method_obs[] == "butterworth" ? "iir" : "fir"
                filt_func = filter_func_obs[]

                try
                    if lowpass_enabled[]
                        lowpass_filter!(
                            filtered_erp,
                            lowpass_cutoff[];
                            order = lowpass_order[],
                            filter_method = filt_method,
                            filter_func = filt_func,
                        )
                    end
                    if highpass_enabled[]
                        highpass_filter!(
                            filtered_erp,
                            highpass_cutoff[];
                            order = highpass_order[],
                            filter_method = filt_method,
                            filter_func = filt_func,
                        )
                    end

                    # Unmirror if enabled
                    if mirror_enabled[]
                        unmirror!(filtered_erp, :both)
                    end

                    # Plot filtered ERP in red (solid line)
                    plot_erp!(fig, ax, [filtered_erp], channel_selection = channels(selected_channel[]), legend = false, color = :red)

                    # Add legend to this subplot
                    original_line = filter(plot_obj -> plot_obj isa Lines, ax.scene.plots)[1]
                    filtered_line = filter(plot_obj -> plot_obj isa Lines, ax.scene.plots)[2]
                    axislegend(ax, [original_line, filtered_line], ["Original", "Filtered"], position = :rt, framevisible = false)
                catch e
                    if e isa BoundsError
                        # FIR filter too long for data
                        text!(
                            ax,
                            0.5,
                            0.5,
                            text = "Filter error:\nFIR too long for data.\nUse IIR instead.",
                            align = (:center, :center),
                            fontsize = 10,
                            color = :red,
                        )
                    else
                        rethrow(e)
                    end
                end
            end
        end
    end

    # Initial plot
    update_plot!()

    # ===== CONNECT OBSERVABLES =====
    on(channel_menu.selection) do ch
        selected_channel[] = ch
        update_plot!()
    end

    on(method_menu.selection) do method_pair
        method_str = method_pair isa Pair ? method_pair[2] : method_pair
        filter_method_obs[] = method_str
        update_plot!()
    end

    on(func_menu.selection) do func_pair
        func_str = func_pair isa Pair ? func_pair[2] : func_pair
        filter_func_obs[] = func_str
        update_plot!()
    end

    on(mirror_toggle.active) do val
        mirror_enabled[] = val
        update_plot!()
    end

    # Lowpass controls
    on(lowpass_toggle.active) do val
        lowpass_enabled[] = val
        update_plot!()
    end

    on(lowpass_cutoff_slider.value) do val
        lowpass_cutoff[] = val
        lowpass_cutoff_label.text = @sprintf("%.1f Hz", val)
        update_plot!()
    end

    on(lowpass_order_slider.value) do val
        lowpass_order[] = val
        lowpass_order_label.text = "$(val)"
        update_plot!()
    end

    # Highpass controls
    on(highpass_toggle.active) do val
        highpass_enabled[] = val
        update_plot!()
    end

    on(highpass_cutoff_slider.value) do val
        highpass_cutoff[] = val
        highpass_cutoff_label.text = @sprintf("%.1f Hz", val)
        update_plot!()
    end

    on(highpass_order_slider.value) do val
        highpass_order[] = val
        highpass_order_label.text = "$(val)"
        update_plot!()
    end

    display(fig)
    return fig
end
