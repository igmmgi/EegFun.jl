# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_EPOCHS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Data processing
    :average_channels => (false, "Whether to average across channels"),
    :plot_avg_trials => (true, "Whether to draw ERP average overlay"),

    # Axis limits, labels, and direction
    :title => (nothing, "Plot title. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :yreversed => (false, "Whether to reverse the y-axis"),

    # Line styling
    :linewidth => ([1, 2], "Line width for epoch traces and average"),
    :color => ([:grey, :red], "Colors for epoch traces and average"),
    :alpha => ([0.3, 1.0], "Transparency for epoch traces and average"),

    # Layout configuration
    :layout => (:single, "Layout type: :single, :grid, or :topo"),
    :layout_plot_width => (0.14, "Width of individual plots in layout"),
    :layout_plot_height => (0.14, "Height of individual plots in layout"),
    :layout_margin => (0.02, "Margin between plots in layout"),
    :layout_show_scale => (true, "Whether to show scale in layout"),
    :layout_scale_position => ([0.99, 0.01], "Position of scale in layout"),
    :dims => (nothing, "Grid dimensions as (rows, cols). If nothing, automatically determined"),

    # Display options
    :theme_fontsize => (24, "Font size for theme"),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Origin lines
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),

    # Interactive features
    :interactive => (true, "Whether to enable interactive features"),
)

"""
    plot_epochs(dat::EpochData; 
                channel_selection::Function = channels(),
                sample_selection::Function = samples(), 
                epoch_selection::Function = epochs(),
                include_extra::Bool = false,
                layout = :single,
                kwargs...)

Plot epoch data with flexible layout options.

# Arguments
- `dat::EpochData`: Epoch data structure containing multiple trials
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: `channels()`)
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: `samples()`)
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering (default: `epochs()`)
- `include_extra::Bool`: Whether to include extra channels (default: `false`)
- `layout`: Layout specification:
  - `:single` (default): Single plot with all channels
  - `:grid`: Auto-calculated grid layout
  - `:topo`: Topographic layout based on channel positions
  - `Vector{Int}`: Custom grid dimensions [rows, cols]

$(generate_kwargs_doc(PLOT_EPOCHS_KWARGS))

# Returns
- `Figure`: The Makie Figure object
- `Union{Axis, Vector{Axis}}`: Single axis for single layout, or vector of axes for grid/topo layouts

# Examples
```julia
# Single plot with all channels
fig, ax = plot_epochs(dat)

# Grid layout
fig, axes = plot_epochs(dat, layout = :grid)

# Custom grid dimensions
fig, axes = plot_epochs(dat, layout = [2, 3])

# Don't display plot
fig, ax = plot_epochs(dat; display_plot = false)

# Custom styling
fig, ax = plot_epochs(dat; 
    color = [:blue, :red],
    linewidth = [1, 3],
    title = "Custom Epochs"
)
```
"""
function plot_epochs(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,  # :single (default), :grid, :topo, or [rows, cols]
    kwargs...,
)::Tuple{Figure,Union{Axis,Vector{Axis}}}

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)

    # Use subset to get the data we want to plot (same pattern as other functions)
    dat_subset = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )

    # Get all non-metadata channels from the subset (it already contains only what we want)
    selected_channels = channel_labels(dat_subset)  # Gets EEG channels from layout
    extra_channels = extra_labels(dat_subset)       # Gets extra channels (EOG, etc.)
    all_plot_channels = vcat(selected_channels, extra_channels)

    # Validate we have channels to plot
    isempty(all_plot_channels) && throw(ArgumentError("No channels selected for plotting"))

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)

    # Info about what we're plotting
    @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(length(dat_subset.data)) epochs"

    fig = Figure()
    axes = Axis[]  # Keep track of all axes created

    # Apply theme font size early to affect all elements
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])

    # Use single plot if explicitly requested via average_channels
    if plot_kwargs[:average_channels]

        # Single plot averaging across channels
        ax = Axis(fig[1, 1])
        push!(axes, ax)

        # Precompute averaged channel once via channel_average!
        dat_avg = copy(dat_subset)
        channel_average!(
            dat_avg;
            channel_selections = [channels(all_plot_channels)],
            output_labels = [:avg],
            reduce = false,
        )
        _plot_epochs!(ax, dat_avg, [:avg], plot_kwargs)

        # Optional ERP overlay (compute only when needed) from averaged data
        if plot_kwargs[:plot_avg_trials]
            erp_dat = average_epochs(dat_avg)
            _plot_erp_average!(ax, erp_dat, [:avg], plot_kwargs)
        end

        # Draw axes through origin if requested
        # _set_origin_lines!(ax; add_xy_origin = get(plot_kwargs, :add_xy_origin, false))

        # Set axis properties
        if length(all_plot_channels) == 1
            ax.title = string(all_plot_channels[1])
        else
            ax.title = "Avg: $(print_vector(all_plot_channels, max_length = 8, n_ends = 3))"
        end
        
        # Use the existing refactored helper functions
        _set_axis_properties!(ax; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim], 
                           xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], yreversed = plot_kwargs[:yreversed])
        _set_axis_grid!(ax; 
                         xgrid = plot_kwargs[:xgrid], 
                         ygrid = plot_kwargs[:ygrid],
                         xminorgrid = plot_kwargs[:xminorgrid], 
                         yminorgrid = plot_kwargs[:yminorgrid])
        _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])

    else

        # Separate subplot for each channel
        n_channels = length(all_plot_channels)
        
        # Compute ERP once for all layout types (avoids duplication)
        erp_dat = plot_kwargs[:plot_avg_trials] ? average_epochs(dat_subset) : nothing
        
        # Store grid dimensions for reuse in interactive section
        grid_dims = nothing
        
        # Handle layout parameter: :single, :grid, :topo, or custom dims
        if layout === :topo
            _plot_epochs_topo!(fig, axes, dat_subset, erp_dat, all_plot_channels, plot_kwargs)
        elseif layout === :grid || !isnothing(plot_kwargs[:dims])
            # Determine grid dimensions
            grid_dims = if !isnothing(plot_kwargs[:dims])
                dims = plot_kwargs[:dims]
                if dims[1] * dims[2] < n_channels
                    throw(ArgumentError("dims grid ($(dims[1])×$(dims[2])=$(dims[1]*dims[2])) is too small for $n_channels channels"))
                end
                dims
            else
                best_rect(n_channels)
            end
            _plot_epochs_grid!(fig, axes, dat_subset, erp_dat, all_plot_channels, grid_dims, plot_kwargs)
        elseif layout === :single
            _plot_epochs_single!(fig, axes, dat_subset, erp_dat, all_plot_channels, plot_kwargs)
        else
            throw(ArgumentError("layout must be :single, :grid, or :topo"))
        end
    end

    # Link axes (both x and y) when multiple axes are present
    if length(axes) > 1
        Makie.linkaxes!(axes...)
    end

    # Add interactive functionality (keyboard zoom + time/channel selection)
    if plot_kwargs[:interactive]
        # Disable default Makie interactions that would conflict with our custom handling
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end

        # Setup keyboard interactivity (zoom with arrow keys)
        _setup_shared_interactivity!(fig, axes, :epochs)

        # Add unified selection functionality (time + channel selection)
        # Create a selection state for both time and channel selection
        selection_state = SharedSelectionState(axes)

        # Set up unified selection (time selection with Shift, channel selection with Ctrl)
        if layout === :topo
            _setup_interactivity_topo!(fig, axes, selection_state, dat_subset, all_plot_channels)
        elseif layout === :grid || !isnothing(plot_kwargs[:dims])
            _setup_interactivity_grid!(fig, axes, selection_state, dat_subset, all_plot_channels, grid_dims)
        elseif layout === :single
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, nothing)
        end
    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    # Return fig and all axes (or single axis if only one)
    return fig, length(axes) == 1 ? first(axes) : axes

end

function _plot_epochs!(ax, dat, channels, plot_kwargs)::Nothing
    # This function expects exactly one channel; callers pass [:avg] or [channel]
    @info "plot_epochs: $(print_vector(channels))"
    @assert length(channels) == 1 "_plot_epochs! expects a single channel"
    ch = channels[1]

    # Cache time vector and styles
    time_vec = dat.data[1][!, :time]
    trial_color = plot_kwargs[:color][1]
    trial_linewidth = plot_kwargs[:linewidth][1]

    # Concatenate all trials with NaN separators into single buffers
    trials = dat.data
    m = length(trials)
    n = length(time_vec)
    total_len = m * n + (m - 1)

    time_cat = Vector{Float64}(undef, total_len)
    y_cat = Vector{Float64}(undef, total_len)

    pos = 1
    @inbounds for t = 1:m
        df = trials[t]
        y = df[!, ch]
        for i = 1:n
            time_cat[pos] = time_vec[i]
            y_cat[pos] = y[i]
            pos += 1
        end
        if t != m
            time_cat[pos] = NaN
            y_cat[pos] = NaN
            pos += 1
        end
    end

    lines!(ax, time_cat, y_cat, color = trial_color, linewidth = trial_linewidth)

    return nothing
end


function _plot_erp_average!(ax, erp_dat::ErpData, channels::Vector{Symbol}, plot_kwargs)::Nothing
    @assert length(channels) == 1 "_plot_erp_average! expects a single channel"
    ch = channels[1]

    time_vec = erp_dat.data[!, :time]
    avg_color = plot_kwargs[:color][2]
    avg_linewidth = plot_kwargs[:linewidth][2]

    lines!(ax, time_vec, erp_dat.data[!, ch], color = avg_color, linewidth = avg_linewidth)
    return nothing
end

"""
    _plot_epochs_topo!(fig, axes, dat, erp_dat, all_plot_channels, kwargs)

Create a topographic layout for plotting epochs.
"""
function _plot_epochs_topo!(
    fig::Figure,
    axes::Vector{Axis},
    dat::EpochData,
    erp_dat,
    all_plot_channels::Vector{Symbol},
    plot_kwargs::Dict,
)
    # Ensure 2D coordinates exist
    if !all(in.([:x2, :y2], Ref(propertynames(dat.layout.data))))
        polar_to_cartesian_xy!(dat.layout)
    end

    # Determine global y-lims for consistency across small axes
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(dat; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    # Normalize positions to [0,1]
    x2 = dat.layout.data.x2
    y2 = dat.layout.data.y2
    minx, maxx = extrema(x2)
    miny, maxy = extrema(y2)
    xrange = maxx - minx
    yrange = maxy - miny
    xrange = xrange == 0 ? 1.0 : xrange
    yrange = yrange == 0 ? 1.0 : yrange

    plot_w = get(plot_kwargs, :layout_plot_width, 0.12)
    plot_h = get(plot_kwargs, :layout_plot_height, 0.12)
    margin = get(plot_kwargs, :layout_margin, 0.02)

    # Map channel -> position
    pos_map = Dict{Symbol,Tuple{Float64,Float64}}()
    for (lab, x, y) in zip(dat.layout.data.label, x2, y2)
        nx = (x - minx) / xrange
        ny = (y - miny) / yrange
        pos_map[Symbol(lab)] = (nx, ny)
    end

    # Create axes at positions
    for ch in all_plot_channels
        pos = get(pos_map, ch, (0.5, 0.5))
        ax = Axis(
            fig[1, 1],
            width = Relative(plot_w),
            height = Relative(plot_h),
            halign = clamp(pos[1], margin, 1 - margin),
            valign = clamp(pos[2], margin, 1 - margin),
        )
        push!(axes, ax)
        _plot_epochs!(ax, dat, [ch], plot_kwargs)
        erp_dat !== nothing && _plot_erp_average!(ax, erp_dat, [ch], plot_kwargs)


        # Suppress axis labels on all but the final axis; set only limits and title for now
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim, :xlabel => "", :ylabel => ""))
        ax.title = string(ch)
        _set_axis_properties!(ax; xlim = axis_kwargs[:xlim], ylim = axis_kwargs[:ylim],
                           xlabel = axis_kwargs[:xlabel], ylabel = axis_kwargs[:ylabel], yreversed = axis_kwargs[:yreversed])
        _set_axis_grid!(ax; 
                         xgrid = axis_kwargs[:xgrid], 
                         ygrid = axis_kwargs[:ygrid],
                         xminorgrid = axis_kwargs[:xminorgrid], 
                         yminorgrid = axis_kwargs[:yminorgrid])
        _set_origin_lines!(ax; add_xy_origin = axis_kwargs[:add_xy_origin])
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticksvisible = false
        hidespines!(ax)  # Hide spines for clean topographic layout (like plot_erp)
    end

    # Optional extra scale axis in bottom-right
    if plot_kwargs[:layout_show_scale]
        scale_ax = Axis(fig[1, 1], 
            width = Relative(plot_kwargs[:layout_plot_width]), 
            height = Relative(plot_kwargs[:layout_plot_height]),
            halign = plot_kwargs[:layout_scale_position][1], 
            valign = plot_kwargs[:layout_scale_position][2])
        push!(axes, scale_ax)
        # No data in this axis; just show labels and limits
        tmin, tmax = (dat.data[1].time[1], dat.data[1].time[end])
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim, :xlim => (tmin, tmax)))
        scale_ax.title = ""
        _set_axis_properties!(scale_ax; xlim = axis_kwargs[:xlim], ylim = axis_kwargs[:ylim],
                           xlabel = axis_kwargs[:xlabel], ylabel = axis_kwargs[:ylabel], yreversed = axis_kwargs[:yreversed])
        _set_axis_grid!(scale_ax; 
                         xgrid = axis_kwargs[:xgrid], 
                         ygrid = axis_kwargs[:ygrid],
                         xminorgrid = axis_kwargs[:xminorgrid], 
                         yminorgrid = axis_kwargs[:yminorgrid])
        _set_origin_lines!(scale_ax; add_xy_origin = axis_kwargs[:add_xy_origin])
        # scale_ax.xticklabelsvisible = true
        # scale_ax.yticklabelsvisible = true
    end
end

"""
    _plot_epochs_grid!(fig, axes, dat, erp_dat, all_plot_channels, grid_dims, kwargs)

Create a grid layout for plotting epochs.

# Arguments
- `grid_dims::Tuple{Int, Int}`: Grid dimensions as (rows, cols)
"""
function _plot_epochs_grid!(
    fig::Figure,
    axes::Vector{Axis},
    dat::EpochData,
    erp_dat,
    all_plot_channels::Vector{Symbol},
    grid_dims::Tuple{Int, Int},
    plot_kwargs::Dict,
)
    rows, cols = grid_dims

    # Calculate y-range if not provided 
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(dat; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    for (idx, channel) in enumerate(all_plot_channels)
        row = fld(idx-1, cols) + 1
        col = mod(idx-1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)
        _plot_epochs!(ax, dat, [channel], plot_kwargs)
        if erp_dat !== nothing
            _plot_erp_average!(ax, erp_dat, [channel], plot_kwargs)
        end


        # Set axis properties with ylim
        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim))

        # Only add x and y labels to outer left column and bottom row
        if col != 1
            axis_kwargs = merge(axis_kwargs, Dict(:ylabel => ""))
            ax.yticklabelsvisible = false
        end
        if row != rows
            axis_kwargs = merge(axis_kwargs, Dict(:xlabel => ""))
            ax.xticklabelsvisible = false
        end

        # Set axis styling using shared functions
        _set_axis_properties!(ax; xlim = axis_kwargs[:xlim], ylim = axis_kwargs[:ylim],
                           xlabel = axis_kwargs[:xlabel], ylabel = axis_kwargs[:ylabel], yreversed = axis_kwargs[:yreversed])
        _set_axis_grid!(ax; 
                         xgrid = axis_kwargs[:xgrid], 
                         ygrid = axis_kwargs[:ygrid],
                         xminorgrid = axis_kwargs[:xminorgrid], 
                         yminorgrid = axis_kwargs[:yminorgrid])
        _set_origin_lines!(ax; add_xy_origin = axis_kwargs[:add_xy_origin])
        
        # Set axis title
        ax.title = isnothing(axis_kwargs[:title]) ? "$channel" : axis_kwargs[:title]
    end
end

"""
    _plot_epochs_single!(fig, axes, dat, erp_dat, all_plot_channels, kwargs)

Create a single axis layout for plotting epochs with all channels on the same plot.
"""
function _plot_epochs_single!(
    fig::Figure,
    axes::Vector{Axis},
    dat::EpochData,
    erp_dat,
    all_plot_channels::Vector{Symbol},
    plot_kwargs::Dict,
)
    # Create a single axis
    ax = Axis(fig[1, 1])
    push!(axes, ax)

    # Plot each channel on the same axis
    for channel in all_plot_channels
        _plot_epochs!(ax, dat, [channel], plot_kwargs)
        erp_dat !== nothing && _plot_erp_average!(ax, erp_dat, [channel], plot_kwargs)
    end

    # Set title based on selected channels (like plot_erp does)
    channel_title =
        length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"

    # Set axis properties
    ax.title = channel_title
    _set_axis_properties!(ax; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim],
                       xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], yreversed = plot_kwargs[:yreversed])
    _set_axis_grid!(ax; 
                     xgrid = plot_kwargs[:xgrid], 
                     ygrid = plot_kwargs[:ygrid],
                     xminorgrid = plot_kwargs[:xminorgrid], 
                     yminorgrid = plot_kwargs[:yminorgrid])
    _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
end

"""
    _setup_interactivity_topo!(fig, axes, selection_state, dat, all_plot_channels)

Set up topographic layout interactivity (unified selection + channel selection).
"""
function _setup_interactivity_topo!(fig, axes, selection_state, dat, all_plot_channels)
    topo_layout = create_topo_layout(dat.layout, all_plot_channels)
    _setup_unified_selection!(fig, axes, selection_state, dat, topo_layout)
    _setup_channel_selection_events!(fig, selection_state, topo_layout, dat, axes, :topo)
end

"""
    _setup_interactivity_grid!(fig, axes, selection_state, dat, all_plot_channels, grid_dims)

Set up grid layout interactivity (unified selection + channel selection).
"""
function _setup_interactivity_grid!(fig, axes, selection_state, dat, all_plot_channels, grid_dims)
    grid_layout = create_grid_layout(all_plot_channels, rows = grid_dims[1], cols = grid_dims[2])
    _setup_unified_selection!(fig, axes, selection_state, dat, grid_layout)
    _setup_channel_selection_events!(fig, selection_state, grid_layout, dat, axes, :grid)
end


