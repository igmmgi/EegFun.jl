# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_EPOCHS_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    
    # Data processing
    :average_channels => (false, "Whether to average across channels"),
    :plot_avg_trials => (true, "Whether to draw ERP average overlay"),
    
    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :title => (nothing, "Plot title. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    
    # Line styling
    :linewidth => ([1, 2], "Line width for epoch traces and average"),
    :color => ([:grey, :red], "Colors for epoch traces and average"),
    :alpha => ([0.3, 1.0], "Transparency for epoch traces and average"),
    :yreversed => (false, "Whether to reverse the y-axis"),
    
    # Layout configuration
    :layout => (:single, "Layout type: :single, :grid, or :topo"),
    :layout_plot_width => (0.12, "Width of individual plots in layout"),
    :layout_plot_height => (0.12, "Height of individual plots in layout"),
    :layout_margin => (0.02, "Margin between plots in layout"),
    :layout_show_scale => (true, "Whether to show scale in layout"),
    :layout_scale_position => ([0.95, 0.05], "Position of scale in layout"),
    :layout_scale_width => (0.14, "Width of scale in layout"),
    :layout_scale_height => (0.14, "Height of scale in layout"),
    :dims => (nothing, "Grid dimensions as (rows, cols). If nothing, automatically determined"),
    
    # Display options
    :hidedecorations => (false, "Whether to hide axis decorations"),
    :theme_fontsize => (24, "Font size for theme"),
    
    # Legend
    :legend => (true, "Whether to show the legend"),
    :legend_label => ("", "Custom label for the legend"),
    
    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    
    # Origin lines
    :axes_through_origin => (true, "Whether to add origin lines at x=0 and y=0"),
    
    # Minor grid
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    
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
    kwargs...)::Tuple{Figure, Union{Axis, Vector{Axis}}}
    
    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)
    
    # Use subset to get the data we want to plot (same pattern as other functions)
    dat_subset = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra
    )
    
    # Get all non-metadata channels from the subset (it already contains only what we want)
    selected_channels = channel_labels(dat_subset)  # Gets EEG channels from layout
    extra_channels = extra_labels(dat_subset)       # Gets extra channels (EOG, etc.)
    all_plot_channels = vcat(selected_channels, extra_channels)
    
    # Validate we have channels to plot
    isempty(all_plot_channels) && throw(ArgumentError("No channels selected for plotting"))
    
    # Info about what we're plotting
    @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(length(dat_subset.data)) epochs"

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)

    fig = Figure()
    axes = Axis[]  # Keep track of all axes created

    # Early branch: average across channels path
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
            _plot_epochs_from_erp!(ax, erp_dat, [:avg], plot_kwargs)
        end

        # Draw axes through origin if requested
        if get(plot_kwargs, :axes_through_origin, false)
            hlines!(ax, [0.0], color = :black, linewidth = 1)
            vlines!(ax, [0.0], color = :black, linewidth = 1)
        end

        # Set axis properties
        if length(all_plot_channels) == 1
            _set_axis_properties!(ax, plot_kwargs, "$(all_plot_channels[1])")
        else
            _set_axis_properties!(ax, plot_kwargs, "Avg: $(print_vector(all_plot_channels, max_length = 8, n_ends = 3))")
        end

    else

        # Separate subplot for each channel
        n_channels = length(all_plot_channels)
        # Handle layout parameter: :single, :grid, :topo, or custom [rows, cols]
        
        if layout === :topo
            # Topographic layout using data's layout semantics
            erp_dat = plot_kwargs[:plot_avg_trials] ? average_epochs(dat_subset) : nothing
            _plot_epochs_layout!(fig, axes, dat_subset, erp_dat, all_plot_channels, plot_kwargs)
        elseif layout === :grid
            # Auto grid layout
            rows, cols = best_rect(n_channels)
            erp_dat = plot_kwargs[:plot_avg_trials] ? average_epochs(dat_subset) : nothing
            _plot_epochs_grid!(fig, axes, dat_subset, erp_dat, all_plot_channels, rows, cols, plot_kwargs)
        elseif layout === :single
            # Single plot - create one axis with all channels
            erp_dat = plot_kwargs[:plot_avg_trials] ? average_epochs(dat_subset) : nothing
            
            # Create a single axis
            ax = Axis(fig[1, 1])
            push!(axes, ax)
            
            # Plot each channel on the same axis
            for channel in all_plot_channels
                _plot_epochs!(ax, dat_subset, [channel], plot_kwargs)
                erp_dat !== nothing && _plot_epochs_from_erp!(ax, erp_dat, [channel], plot_kwargs)
            end
            
            # Draw axes through origin if requested
            if get(plot_kwargs, :axes_through_origin, false)
                hlines!(ax, [0.0], color = :black, linewidth = 1)
                vlines!(ax, [0.0], color = :black, linewidth = 1)
            end
            
            # Set title based on selected channels (like plot_erp does)
            channel_title = length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"
            
            # Set axis properties
            _set_axis_properties!(ax, plot_kwargs, channel_title)
        elseif typeof(layout) <: Vector{<:Integer}
            # Custom grid dimensions [rows, cols]
            if length(layout) != 2 || any(x -> x <= 0, layout)
                throw(ArgumentError("layout must be a 2-element vector [rows, cols] with positive integers"))
            end
            rows, cols = layout
            if rows * cols < n_channels
                throw(ArgumentError("layout grid ($(rows)×$(cols)=$(rows*cols)) is too small for $n_channels channels"))
            end
            erp_dat = plot_kwargs[:plot_avg_trials] ? average_epochs(dat_subset) : nothing
            _plot_epochs_grid!(fig, axes, dat_subset, erp_dat, all_plot_channels, rows, cols, plot_kwargs)
        else
            throw(ArgumentError("layout must be :single, :grid, :topo, or a 2-element vector [rows, cols]"))
        end
    end

    # Link axes (both x and y) when multiple axes are present
    if length(axes) > 1
        Makie.linkaxes!(axes...)
    end
    
    # Add keyboard interactivity (zoom with arrow keys)
    if plot_kwargs[:interactive]
        # Disable default Makie interactions that would conflict with our custom handling
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end
        
        _setup_shared_interactivity!(fig, axes)
    end
    
    # Add unified selection functionality (time + channel selection)
    if plot_kwargs[:interactive]
        # Create a selection state for both time and channel selection
        selection_state = SharedSelectionState(axes)
        
        # Set up unified selection (time selection with Shift, channel selection with Ctrl)
        if layout === :topo
            # Topographic layout - use unified selection + topo channel selection
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, create_topo_layout(dat_subset.layout, all_plot_channels))
            _setup_channel_selection_events!(fig, selection_state, create_topo_layout(dat_subset.layout, all_plot_channels), dat_subset, axes, :topo)
        elseif layout === :grid || typeof(layout) <: Vector{<:Integer}
            # Grid layout - use unified selection + grid channel selection
            if layout === :grid
                rows, cols = best_rect(length(all_plot_channels))
            else
                rows, cols = layout
            end
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, create_grid_layout(all_plot_channels, rows = rows, cols = cols))
            _setup_channel_selection_events!(fig, selection_state, create_grid_layout(all_plot_channels, rows = rows, cols = cols), dat_subset, axes, :grid)
        elseif layout === :single
            # Single plot layout - only time selection, no channel selection needed
            _setup_unified_selection!(fig, axes, selection_state, dat_subset, nothing)
        end
    end

    # Theme adjustments
    fontsize_theme = Theme(fontsize = plot_kwargs[:theme_fontsize])
    update_theme!(fontsize_theme)

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
    @inbounds for t in 1:m
        df = trials[t]
        y = df[!, ch]
        for i in 1:n
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


function _plot_epochs_from_erp!(ax, erp_dat::ErpData, channels::Vector{Symbol}, plot_kwargs)::Nothing
    @assert length(channels) == 1 "_plot_epochs_from_erp! expects a single channel"
    ch = channels[1]

    time_vec = erp_dat.data[!, :time]
    avg_color = plot_kwargs[:color][2]
    avg_linewidth = plot_kwargs[:linewidth][2]

    lines!(ax, time_vec, erp_dat.data[!, ch], color = avg_color, linewidth = avg_linewidth)
    return nothing
end

function _plot_epochs_layout!(fig::Figure, axes::Vector{Axis}, dat::EpochData, erp_dat, all_plot_channels::Vector{Symbol}, kwargs::Dict)
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

    plot_w = get(kwargs, :layout_plot_width, 0.12)
    plot_h = get(kwargs, :layout_plot_height, 0.12)
    margin = get(kwargs, :layout_margin, 0.02)

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
        erp_dat !== nothing && _plot_epochs_from_erp!(ax, erp_dat, [ch], plot_kwargs)

        # Draw axes through origin if requested
        if plot_kwargs[:axes_through_origin]
            hlines!(ax, [0.0], color = :black, linewidth = 1)
            vlines!(ax, [0.0], color = :black, linewidth = 1)
        end

        # Suppress axis labels on all but the final axis; set only limits and title for now
        axis_kwargs = merge(kwargs, Dict(:ylim => ylim, :xlabel => "", :ylabel => ""))
        _set_axis_properties!(ax, axis_kwargs, string(ch))
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticksvisible = false
        hidespines!(ax)  # Hide spines for clean topographic layout (like plot_erp)
    end

    # Optional extra scale axis in bottom-right
    if get(kwargs, :layout_show_scale, true)
        sp = get(kwargs, :layout_scale_position, [0.95, 0.05])
        sw = get(kwargs, :layout_scale_width, 0.14)
        sh = get(kwargs, :layout_scale_height, 0.14)
        scale_ax = Axis(
            fig[1, 1],
            width = Relative(sw),
            height = Relative(sh),
            halign = sp[1],
            valign = sp[2],
        )
        push!(axes, scale_ax)
        # No data in this axis; just show labels and limits
        tmin, tmax = (dat.data[1].time[1], dat.data[1].time[end])
        axis_kwargs = merge(kwargs, Dict(:ylim => ylim, :xlim => (tmin, tmax)))
        _set_axis_properties!(scale_ax, axis_kwargs, "")
        scale_ax.xticklabelsvisible = true
        scale_ax.yticklabelsvisible = true
    end
end

"""
    _plot_epochs_grid!(fig, axes, dat, erp_dat, all_plot_channels, rows, cols, kwargs)

Create a grid layout for plotting epochs.
"""
function _plot_epochs_grid!(fig::Figure, axes::Vector{Axis}, dat::EpochData, erp_dat, all_plot_channels::Vector{Symbol}, rows::Int, cols::Int, kwargs::Dict)
    n_channels = length(all_plot_channels)
    
    # Calculate y-range if not provided (using shared yrange helper)
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
        erp_dat !== nothing && _plot_epochs_from_erp!(ax, erp_dat, [channel], plot_kwargs)

        # Draw axes through origin if requested
        if plot_kwargs[:axes_through_origin]
            hlines!(ax, [0.0], color = :black, linewidth = 1)
            vlines!(ax, [0.0], color = :black, linewidth = 1)
        end

        # Set axis properties with ylim
        axis_kwargs = merge(kwargs, Dict(:ylim => ylim))

        # Only add x and y labels to outer left column and bottom row
        if col != 1
            axis_kwargs = merge(axis_kwargs, Dict(:ylabel => ""))
            ax.yticklabelsvisible = false
        end
        if row != rows
            axis_kwargs = merge(axis_kwargs, Dict(:xlabel => ""))
            ax.xticklabelsvisible = false
        end

        _set_axis_properties!(ax, axis_kwargs, "$channel")
    end
end

function _set_axis_properties!(ax::Axis, kwargs::Dict, default_title::String)::Nothing
    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    ax.title = isnothing(kwargs[:title]) ? default_title : kwargs[:title]
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]
    
    # Apply grid settings
    ax.xgridvisible = kwargs[:xgrid]
    ax.ygridvisible = kwargs[:ygrid]
    ax.xminorgridvisible = kwargs[:xminorgrid]
    ax.yminorgridvisible = kwargs[:yminorgrid]
    
    return nothing
end
