
# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ERP_IMAGE_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("ERP Image Plot", "Title for the plot window"),
    :interactive => (true, "Enable interactive features (true/false)"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits for ERP waveform as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("Epoch", "Label for y-axis"),
    :yreversed => (false, "Whether to reverse the y-axis"),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Show title (true/false)"),

    # Image styling
    :colormap => (:jet, "Colormap for the image"),
    :colorrange => (nothing, "Color range for the image. If nothing, automatically determined"),

    # ERP overlay
    :plot_erp => (true, "Whether to plot ERP average overlay"),
    :boxcar_average => (1, "Boxcar average window size for smoothing the ERP image (1 = no smoothing)"),

    # Colorbar
    :colorbar_plot => (true, "Whether to show the colorbar"),
    :colorbar_position => ((1, 2), "Position of the colorbar as (row, col) tuple"),
    :colorbar_width => (30, "Width of the colorbar in pixels"),
    :colorbar_label => ("μV", "Label for the colorbar"),
    :colorbar_plot_numbers =>
        ([], "Plot indices for which to show colorbars. Empty list shows colorbars for all plots."),

    # Grid
    :xgrid => (false, "Show x-axis grid (true/false)"),
    :ygrid => (false, "Show y-axis grid (true/false)"),
    :xminorgrid => (false, "Show x-axis minor grid (true/false)"),
    :yminorgrid => (false, "Show y-axis minor grid (true/false)"),

    # Origin lines
    :add_xy_origin => (true, "Add origin lines at x=0 and y=0 (true/false)"),

    # Layout parameters (for topo and other layouts)
    :layout_topo_plot_width => (0.07, "Width of individual plots (fraction of figure width)"),
    :layout_topo_plot_height => (0.07, "Height of individual plots (fraction of figure height)"),
    :layout_topo_scale_pos => ((0.95, 0.05), "Fallback position for scale plot in topo layout as (x, y) tuple"),

    # Grid layout parameters
    :layout_grid_rowgap => (10, "Gap between rows (in pixels)"),
    :layout_grid_colgap => (10, "Gap between columns (in pixels)"),
    :layout_grid_dims =>
        (nothing, "Grid dimensions as (rows, cols) tuple for grid layouts. If nothing, automatically determined"),
    :layout_grid_skip_positions =>
        (nothing, "Positions to skip in grid layout as vector of (row, col) tuples, e.g., [(2,1), (2,3)]"),

    # General layout parameters
    :figure_padding =>
        ((10, 10, 10, 10), "Padding around entire figure as (left, right, top, bottom) tuple (in pixels)"),
)

"""
    plot_erp_image(dat::EpochData; 
                   layout::Union{Symbol, PlotLayout} = :single,
                   channel_selection::Function = channels(),
                   sample_selection::Function = samples(),
                   kwargs...)

Plot ERP image for specified channels and samples with flexible layout options.

# Arguments
- `dat::EpochData`: Epoch data structure
- `layout`: Layout specification:
  - `:single` (default): Single ERP image plot
  - `:grid`: Auto-calculated grid layout for multiple channels
  - `:topo`: Topographic layout based on channel positions
  - `PlotLayout`: Custom layout object
  - `Vector{Int}`: Custom grid dimensions [rows, cols]
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)

$(generate_kwargs_doc(PLOT_ERP_IMAGE_KWARGS))

# Returns
- `Figure`: The Makie Figure object
- `Vector{Axis}`: Vector of axes for the ERP images

# Examples
```julia
# Plot all channels and samples
plot_erp_image(dat)

# Plot specific channels
plot_erp_image(dat, channel_selection = channels([:Fp1, :Fp2]))

# Exclude reference channels
plot_erp_image(dat, channel_selection = channels_not([:M1, :M2]))

# Plot frontal channels only
plot_erp_image(dat, channel_selection = channels(1:10))

# Custom predicates
plot_erp_image(dat, 
    channel_selection = x -> startswith.(string.(x), "F"),
    sample_selection = x -> x .> 0.0  # Only positive time points
)
```
"""
function plot_erp_image(
    dat::EpochData;
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    kwargs...,
)
    # Merge default kwargs with provided kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_IMAGE_KWARGS, kwargs)

    # Use subset to get the data we want to plot (same pattern as other functions)
    dat_subset = subset(dat; channel_selection = channel_selection, sample_selection = sample_selection)

    # Generate window title from dataset
    title_str = _generate_window_title(dat_subset)
    set_window_title(title_str)

    # Get all available channels
    selected_channels = channel_labels(dat_subset)  # Gets EEG channels from layout
    extra_channels = extra_labels(dat_subset)       # Gets extra channels (EOG, etc.)
    all_plot_channels = vcat(selected_channels, extra_channels)

    # Validate we have channels to plot
    if isempty(all_plot_channels)
        @minimal_error("No channels available for plotting after filtering")
    end

    # Set default plot title only for single layouts (same as plot_erp)
    # For grid/topo layouts, we want individual channel names, not a global title
    if plot_kwargs[:show_title] && plot_kwargs[:title] == "" && layout == :single
        plot_kwargs[:title] =
            length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"
    end

    # Handle colorbar positioning for grid layouts
    colorbar_enabled = plot_kwargs[:colorbar_plot]
    user_colorbar_position = plot_kwargs[:colorbar_position]
    colorbar_plot_numbers = plot_kwargs[:colorbar_plot_numbers]

    # Extract layout_* parameters, remove prefix, and pass to create_layout
    layout_kwargs = _extract_layout_kwargs(plot_kwargs)

    # Create figure and apply layout system
    fig = Figure()

    # For grid layouts with colorbars, we need to expand the grid
    if layout == :grid && colorbar_enabled
        # Get the grid dimensions that would be created
        temp_layout = create_layout(layout, all_plot_channels, dat_subset.layout; layout_kwargs...)
        rows, cols = temp_layout.dims

        # Expand grid to accommodate colorbars (default: to the right)
        if user_colorbar_position !== nothing && user_colorbar_position isa Tuple
            cb_row_offset, cb_col_offset = user_colorbar_position
            if cb_row_offset > 1
                # Colorbars below: double rows
                total_rows = rows * 2
                total_cols = cols
            else
                # Colorbars to the right: double columns
                total_rows = rows
                total_cols = cols * 2
            end
        else
            # Default: colorbars to the right
            total_rows = rows
            total_cols = cols * 2
        end

        # Create a modified layout with expanded dimensions
        # We'll manually create axes in the expanded grid
        plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout; layout_kwargs...)
        axes = Axis[]
        channels = Symbol[]

        # Create axes in the expanded grid
        for (idx, channel) in enumerate(plot_layout.channels)
            base_row = div(idx - 1, cols) + 1
            base_col = mod1(idx, cols)

            if user_colorbar_position !== nothing && user_colorbar_position isa Tuple
                cb_row_offset, cb_col_offset = user_colorbar_position
                if cb_row_offset > 1
                    # Colorbars below
                    plot_row = (base_row - 1) * 2 + 1
                    plot_col = base_col
                else
                    # Colorbars to the right
                    plot_row = base_row
                    plot_col = (base_col - 1) * 2 + 1
                end
            else
                # Default: colorbars to the right
                plot_row = base_row
                plot_col = (base_col - 1) * 2 + 1
            end

            ax = Axis(fig[plot_row, plot_col])
            push!(axes, ax)
            push!(channels, channel)
        end
    else
        # For single or topo layouts, use normal layout system
        plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout; layout_kwargs...)
        axes, channels = _apply_layout!(fig, plot_layout; plot_kwargs...)
    end

    # Disable default interactions that conflict with our custom selection
    if plot_kwargs[:interactive]
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end
    end

    # Plot ERP images for each channel
    heatmaps = []
    for (ax, channel) in zip(axes, channels)
        if plot_layout.type == :single
            # For single layout: average across all selected channels
            data = zeros(length(dat_subset.data), nrow(dat_subset.data[1]))
            for epoch in eachindex(dat_subset.data)
                data[epoch, :] = colmeans(dat_subset.data[epoch], all_plot_channels)
            end
        else
            # For grid/topo layouts: show individual channel
            data = zeros(length(dat_subset.data), nrow(dat_subset.data[1]))
            for epoch in eachindex(dat_subset.data)
                data[epoch, :] = dat_subset.data[epoch][!, channel]
            end
        end

        # Apply boxcar average if requested
        if plot_kwargs[:boxcar_average] > 1
            data = _apply_boxcar_average(data, plot_kwargs[:boxcar_average])
        end

        # Set color range if not provided
        if isnothing(plot_kwargs[:colorrange])
            plot_kwargs[:colorrange] = extrema(data)
        end

        # Create heatmap
        hm = heatmap!(
            ax,
            dat_subset.data[1].time,
            1:length(dat_subset.data),
            transpose(data),
            colormap = plot_kwargs[:colormap],
            colorrange = plot_kwargs[:colorrange],
        )
        push!(heatmaps, hm)

        # Set axis properties (only for single layout)
        if plot_layout.type == :single
            # Don't show x-axis labels on heatmap if ERP trace will be shown below
            ax.xlabel = plot_kwargs[:plot_erp] ? "" : plot_kwargs[:xlabel]
            ax.ylabel = plot_kwargs[:ylabel]
            # Hide x-axis tick labels on heatmap if ERP trace will be shown below
            if plot_kwargs[:plot_erp]
                ax.xticklabelsvisible = false
            end
        elseif plot_layout.type == :topo
            # Set channel name as title for topo layout
            ax.title = string(channel)
        end

        ax.yreversed = plot_kwargs[:yreversed]

        # Set x-axis limits if provided, otherwise use data limits
        if !isnothing(plot_kwargs[:xlim])
            xlims!(ax, plot_kwargs[:xlim])
        else
            xlims!(ax, extrema(dat_subset.data[1].time))
        end

        # Set y-axis limits for ERP image (always 1 to number of epochs)
        ylims!(ax, (1, length(dat_subset.data)))

        # Add colorbar if requested
        if plot_kwargs[:colorbar_plot]
            # Determine if this plot should have a colorbar
            # For single layout, always show if colorbar_plot is true
            # For grid layout, check colorbar_plot_numbers
            should_show_colorbar = if plot_layout.type == :single
                true
            elseif plot_layout.type == :grid
                idx = findfirst(==(channel), channels)
                isempty(colorbar_plot_numbers) || idx in colorbar_plot_numbers
            else
                false
            end

            if should_show_colorbar
                if plot_layout.type == :single
                    # Single layout: use provided position
                    colorbar_position = plot_kwargs[:colorbar_position]
                    Colorbar(
                        fig[colorbar_position...],
                        hm,
                        width = plot_kwargs[:colorbar_width],
                        label = plot_kwargs[:colorbar_label],
                    )
                elseif plot_layout.type == :grid
                    # Grid layout: calculate position based on plot position
                    rows, cols = plot_layout.dims
                    idx = findfirst(==(channel), channels)
                    base_row = div(idx - 1, cols) + 1
                    base_col = mod1(idx, cols)

                    if user_colorbar_position !== nothing && user_colorbar_position isa Tuple
                        cb_row_offset, cb_col_offset = user_colorbar_position
                        if cb_row_offset > 1
                            # Colorbars below
                            plot_row = (base_row - 1) * 2 + 1
                            plot_col = base_col
                            colorbar_row = plot_row + (cb_row_offset - 1)
                            colorbar_col = plot_col + (cb_col_offset - 1)
                        else
                            # Colorbars to the right
                            plot_row = base_row
                            plot_col = (base_col - 1) * 2 + 1
                            colorbar_row = plot_row + (cb_row_offset - 1)
                            colorbar_col = plot_col + (cb_col_offset - 1)
                        end
                    else
                        # Default: colorbar to the right
                        plot_row = base_row
                        plot_col = (base_col - 1) * 2 + 1
                        colorbar_row = plot_row
                        colorbar_col = plot_col + 1
                    end

                    Colorbar(
                        fig[colorbar_row, colorbar_col],
                        hm,
                        width = plot_kwargs[:colorbar_width],
                        label = plot_kwargs[:colorbar_label],
                    )
                end
            end
        end
    end

    # Apply axis properties to all axes (sets labels on all axes first)
    _apply_axis_properties!.(axes; plot_kwargs...)

    # Then apply layout-specific properties (clears labels on inner axes for grid layouts)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...)

    # For plot_erp_image, we want to show tick marks (but not labels) on all axes (including topo layouts)
    # Override the decoration hiding that _apply_layout_axis_properties! does for topo layouts
    if plot_layout.type == :topo
        for ax in axes
            # Show ticks but hide tick labels (keep spines hidden for topo aesthetic)
            ax.xticksvisible = true
            ax.yticksvisible = true
            ax.xticklabelsvisible = false
            ax.yticklabelsvisible = false
        end
    end

    # Add origin lines to all axes
    for ax in axes
        _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
    end

    # Add scale plot for topo layout (similar to plot_epochs)
    # For topo layouts, all plots are positioned in fig[1, 1] with halign/valign
    # We add the scale axis the same way - positioned absolutely in fig[1, 1]
    # This matches how plot_epochs does it and won't interfere with the topo plot layout
    scale_pos = nothing
    if plot_layout.type == :topo
        # Get scale position from kwargs
        scale_pos = plot_kwargs[:layout_topo_scale_pos]

        # Use the same width/height as the channel plots (from layout metadata)
        # This ensures consistency - channel plots use plot_layout.metadata[:topo_plot_width]
        scale_width = plot_layout.metadata[:topo_plot_width]
        scale_height = plot_layout.metadata[:topo_plot_height]

        # Create scale axis positioned at the specified location (axis only, no data)
        # This is positioned absolutely in fig[1, 1] using halign/valign, just like topo plots
        scale_ax = Axis(
            fig[1, 1],
            width = Relative(scale_width),
            height = Relative(scale_height),
            halign = scale_pos[1],
            valign = scale_pos[2],
        )
        push!(axes, scale_ax)

        # Make sure the scale axis is visible - don't hide decorations like other topo axes
        # The scale axis should show all labels and decorations

        # Set up scale axis properties (show time and epoch labels)
        # Use the same approach as plot_epochs
        scale_ax.title = ""
        tmin, tmax = extrema(dat_subset.data[1].time)

        # Set axis limits first
        if !isnothing(plot_kwargs[:xlim])
            xlims!(scale_ax, plot_kwargs[:xlim])
        else
            xlims!(scale_ax, (tmin, tmax))
        end
        ylims!(scale_ax, (1, length(dat_subset.data)))

        # Use _set_axis_properties! like plot_epochs does
        _set_axis_properties!(
            scale_ax;
            xlim = isnothing(plot_kwargs[:xlim]) ? (tmin, tmax) : plot_kwargs[:xlim],
            ylim = (1, length(dat_subset.data)),
            xlabel = plot_kwargs[:xlabel],
            ylabel = plot_kwargs[:ylabel],
            yreversed = plot_kwargs[:yreversed],
        )

        # Use _set_axis_grid! like plot_epochs does
        _set_axis_grid!(
            scale_ax;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )

        # Ensure scale axis spines are visible (unlike other topo axes which have spines hidden)
        scale_ax.bottomspinevisible = true
        scale_ax.topspinevisible = true
        scale_ax.leftspinevisible = true
        scale_ax.rightspinevisible = true

        # Add origin lines like plot_epochs does (if the parameter exists)
        add_xy_origin = get(plot_kwargs, :add_xy_origin, true)
        _set_origin_lines!(scale_ax; add_xy_origin = add_xy_origin)
    end

    # Add colorbar for topo layout AFTER layout is finalized
    # Create it in fig[1, 1] with halign/valign - this is the same approach as the scale axis
    # The key is that all topo elements use absolute positioning within the same grid cell
    # and don't participate in grid size calculations
    if plot_layout.type == :topo && plot_kwargs[:colorbar_plot] && !isempty(heatmaps) && scale_pos !== nothing
        # Calculate colorbar position: place it to the right of the scale axis
        # The scale axis width is Relative(plot_kwargs[:layout_topo_plot_width])
        # We add a small fixed offset (0.02) to position the colorbar next to it
        # scale_width = plot_kwargs[:layout_topo_plot_width]
        colorbar_halign = scale_pos[1] + 0.015  # Position to the right of scale axis

        # Get scale height from layout metadata (same as channel plots and scale axis)
        scale_height = plot_layout.metadata[:topo_plot_height]

        # Create colorbar in fig[1, 1] with halign/valign, positioned to the right of scale axis
        # Use tellwidth=false and tellheight=false to prevent it from affecting grid layout
        Colorbar(
            fig[1, 1],
            heatmaps[1],
            label = plot_kwargs[:colorbar_label],
            width = plot_kwargs[:colorbar_width],
            height = Relative(scale_height),
            halign = colorbar_halign,
            valign = scale_pos[2],
            tellwidth = false,
            tellheight = false,
        )
    end

    # Set up interactivity AFTER heatmaps to ensure rectangles are drawn on top
    if plot_kwargs[:interactive]
        # Set up custom interactivity for ERP images (left/right keys only)
        _setup_erp_image_interactivity!(fig, axes, heatmaps)

        # Set up selection system (rectangles created AFTER heatmaps)
        selection_state = SharedSelectionState(axes)
        _setup_unified_selection!(fig, axes, selection_state, dat_subset, plot_layout)

        # Set up channel selection events for topo and grid layouts
        if plot_layout.type == :topo
            _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes, :topo)
        elseif plot_layout.type == :grid
            _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes, :grid)
        end
    end

    # Add ERP trace if requested (only for single layout)
    if plot_kwargs[:plot_erp] && plot_layout.type == :single
        # Convert epoch data to ERP data by averaging across epochs
        erp_data = average_epochs(dat_subset)

        # Create ERP axis in our main figure
        ax_erp = Axis(fig[2, 1])

        # Add ERP axis to the axes list so it gets keyboard navigation
        push!(axes, ax_erp)

        # For single layout, show averaged ERP trace across all selected channels
        avg_data = colmeans(erp_data.data, all_plot_channels)
        lines!(ax_erp, erp_data.data.time, avg_data, color = :black, linewidth = 1)

        # Set axis properties consistent with plot_erp
        ax_erp.xlabel = plot_kwargs[:xlabel]
        ax_erp.ylabel = "Amplitude (μV)"
        ax_erp.xgridvisible = false
        ax_erp.ygridvisible = false

        # Set y-axis limits for ERP trace
        if !isnothing(plot_kwargs[:ylim])
            ylims!(ax_erp, plot_kwargs[:ylim])
        else
            # Auto-determine y-limits based on data
            ylims!(ax_erp, extrema(avg_data))
        end
        ax_erp.xminorgridvisible = false
        ax_erp.yminorgridvisible = false

        # Add origin lines like plot_erp does
        vlines!(ax_erp, [0], color = :gray, linewidth = 0.5)
        hlines!(ax_erp, [0], color = :gray, linewidth = 0.5)

        # Set limits if provided, otherwise use the same limits as the heatmap
        if !isnothing(plot_kwargs[:xlim])
            xlims!(ax_erp, plot_kwargs[:xlim])
        else
            # Use the same x-limits as the heatmap axis
            heatmap_xlims = axes[1].xaxis.attributes.limits[]
            xlims!(ax_erp, heatmap_xlims)
        end

        # Resize rows: 2/3 for heatmap, 1/3 for ERP trace
        rowsize!(fig.layout, 1, Relative(2/3))
        rowsize!(fig.layout, 2, Relative(1/3))
    end

    # Link only x-axes (time dimension) but not y-axes (epoch dimension)
    # since each channel may have different numbers of epochs
    if length(axes) > 1
        linkxaxes!(axes...)
    end

    # Display plot if requested
    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    set_window_title("Makie")
    return fig, axes
end

"""
    _setup_erp_image_interactivity!(fig::Figure, axes::Vector{Axis}, heatmaps::Vector)

Set up custom interactivity for ERP images with context-aware key controls based on mouse position.
"""
function _setup_erp_image_interactivity!(fig::Figure, axes::Vector{Axis}, heatmaps::Vector)
    _setup_keyboard_tracking(fig)

    # Track mouse position to determine active axis
    active_axis = Ref(axes[1])  # Default to first axis (ERP image)

    # Update active axis based on mouse position
    on(events(fig).mouseposition) do mouse_pos
        for ax in axes
            viewport = ax.scene.viewport[]
            if mouse_pos[1] >= viewport.origin[1] &&
               mouse_pos[1] <= viewport.origin[1] + viewport.widths[1] &&
               mouse_pos[2] >= viewport.origin[2] &&
               mouse_pos[2] <= viewport.origin[2] + viewport.widths[2]
                active_axis[] = ax
                break
            end
        end
    end

    # Define keyboard actions for ERP images
    keyboard_actions =
        Dict(Keyboard.up => :y_more, Keyboard.down => :y_less, Keyboard.left => :x_less, Keyboard.right => :x_more)

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            # Show help for ERP image
            show_plot_help(:erp_image)
        elseif event.action in (Keyboard.press, Keyboard.repeat) && haskey(keyboard_actions, event.key)
            action = keyboard_actions[event.key]
            _handle_erp_image_navigation!(axes, heatmaps, action, active_axis[])
        end
    end
end

"""
    _handle_erp_image_navigation!(axes::Vector{Axis}, heatmaps::Vector, action::Symbol)

Handle navigation actions for ERP images.
"""
function _handle_erp_image_navigation!(axes::Vector{Axis}, heatmaps::Vector, action::Symbol, active_axis::Axis)
    if action in (:x_less, :x_more) # Adjust time axis
        func = action == :x_less ? xmore! : xless!
        func(active_axis)
    elseif action in (:y_less, :y_more) # Adjust y-axis based on active plot
        if active_axis == axes[1] # ERP image axis
            # For ERP image, adjust color scale instead of y-axis
            factor = action == :y_more ? 1.25 : 0.8
            for hm in heatmaps
                hm.colorrange[] = hm.colorrange[] .* factor
            end
        else # ERP waveform axis
            # For ERP waveform, adjust y-axis limits
            func = action == :y_less ? ymore! : yless!
            func(active_axis)
        end
    end
end


"""
    _apply_boxcar_average(data::Matrix, window_size::Int)

Apply boxcar (moving average) smoothing to the ERP image data.

# Arguments
- `data::Matrix`: The ERP image data (epochs × time points)
- `window_size::Int`: Size of the moving average window

# Returns
- `Matrix`: Smoothed data with the same dimensions

# Examples
```julia
# Apply 5-epoch boxcar average
smoothed_data = _apply_boxcar_average(data, 5)
```
"""
function _apply_boxcar_average(data::Matrix, window_size::Int)
    n_epochs, n_timepoints = size(data)
    smoothed_data = similar(data)

    # Calculate half window size for centering
    half_window = window_size ÷ 2

    for epoch = 1:n_epochs
        for timepoint = 1:n_timepoints
            # Define the window boundaries
            start_epoch = max(1, epoch - half_window)
            end_epoch = min(n_epochs, epoch + half_window)

            # Calculate the average within the window
            smoothed_data[epoch, timepoint] = mean(data[start_epoch:end_epoch, timepoint])
        end
    end

    return smoothed_data
end
