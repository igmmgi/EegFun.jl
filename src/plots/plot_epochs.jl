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
    :linewidth => (1, "Line width for epoch traces"),
    :avg_linewidth_multiplier => (2.0, "Multiplier for average line width (average linewidth = linewidth * avg_linewidth_multiplier)"),
    :trial_alpha => (0.25, "Alpha (transparency) for individual trial traces"),
    :color => (:black, "Color for epoch traces (can be a single color or a vector of colors, one per condition)"),
    :colormap => (:jet, "Colormap for multi-condition plots"),

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

    # Legend parameters - get all Legend attributes with their actual defaults
    # This allows users to control any Legend parameter
    [
        Symbol("legend_$(attr)") => (get(LEGEND_DEFAULTS, attr, nothing), "Legend $(attr) parameter") for
        attr in propertynames(Legend)
    ]...,

    # Override specific legend parameters with custom defaults
    :legend => (true, "Whether to show the legend"),
    :legend_label => ("", "Title for the legend"),
    :legend_framevisible => (true, "Whether to show the frame of the legend"),
    :legend_position => (:lt, "Position of the legend for axislegend() (symbol like :lt, :rt, :lb, :rb, or tuple like (:left, :top), or (0.5, 0.5))"),
    :legend_channel => ([], "If plotting multiple plots, within channel to put the legend on."),
    :legend_labels => ([], "If plotting multiple plots, custom labels for conditions."),
    :legend_nbanks => (nothing, "Number of columns for the legend. If nothing, automatically determined."),
)

"""
    plot_epochs(filepath::String; 
               channel_selection::Function = channels(),
               sample_selection::Function = samples(), 
               epoch_selection::Function = epochs(),
               include_extra::Bool = false,
               layout = :single,
               kwargs...)

Load epoch data from a JLD2 file and create plots.

# Arguments
- `filepath::String`: Path to JLD2 file containing EpochData
- `channel_selection::Function`: Function that returns boolean vector for channel filtering
- `sample_selection::Function`: Function that returns boolean vector for sample filtering
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering
- `include_extra::Bool`: Whether to include extra channels
- `layout`: Layout specification (see main plot_epochs documentation)
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Load and plot from file
plot_epochs("Flank_C_3_epochs_original.jld2")

# With channel selection
plot_epochs("Flank_C_3_epochs_original.jld2", channel_selection = channels([:PO7, :PO8]), layout = :grid)
```
"""
function plot_epochs(
    filepath::String;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)
    # Load data from file
    data = load_data(filepath)
    if isnothing(data)
        @minimal_error_throw "No data found in file: $filepath"
    end

    # Handle Vector{EpochData} case (if file contains multiple epoch datasets)
    if data isa Vector && !isempty(data) && first(data) isa EpochData
        # For now, just use the first dataset (could be extended to plot all)
        data = first(data)
    end

    # Dispatch to main plot_epochs function
    return plot_epochs(
        data;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
        layout = layout,
        kwargs...,
    )
end

"""
    plot_epochs(datasets::Vector{EpochData}; 
                condition_selection::Function = conditions(),
                channel_selection::Function = channels(),
                sample_selection::Function = samples(), 
                epoch_selection::Function = epochs(),
                include_extra::Bool = false,
                layout = :single,
                kwargs...)

Plot multiple epoch datasets (conditions) with flexible layout options.
Conditions are overlaid on the same plot with different colors.

# Arguments
- `datasets::Vector{EpochData}`: Vector of epoch data structures (one per condition)
- `condition_selection::Function`: Function that returns boolean vector for condition filtering (default: `conditions()`)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: `channels()`)
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: `samples()`)
- `epoch_selection::Function`: Function that returns boolean vector for epoch filtering (default: `epochs()`)
- `include_extra::Bool`: Whether to include extra channels (default: `false`)
- `layout`: Layout specification (see single EpochData method documentation)
- `kwargs`: Additional keyword arguments

# Examples
```julia
# Plot all conditions
fig, ax = plot_epochs([dat1, dat2])

# Plot specific conditions
fig, ax = plot_epochs([dat1, dat2], condition_selection = conditions([1, 2]))
```
"""
function plot_epochs(
    datasets::Vector{EpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,
    kwargs...,
)::Tuple{Figure,Union{Axis,Vector{Axis}}}
    # Prepare kwargs and track color setting (like plot_erp)
    color_explicitly_set = haskey(kwargs, :color)
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_KWARGS, kwargs)
    plot_kwargs[:_color_explicitly_set] = color_explicitly_set

    # Use subset to filter by condition and apply other selections
    dat_subset = subset(
        datasets;
        condition_selection = condition_selection,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )

    # Check if subsetting resulted in empty data
    if isempty(dat_subset)
        n_conditions = length(datasets)
        if n_conditions == 0
            @minimal_error_throw "No data available (empty dataset)"
        else
            @minimal_error_throw "No data matched the selection criteria. Available condition indices: 1:$n_conditions"
        end
    end

    # Get channels from first dataset (all should have same channels after subsetting)
    selected_channels = channel_labels(dat_subset[1])
    extra_channels = extra_labels(dat_subset[1])
    all_plot_channels = vcat(selected_channels, extra_channels)

    # Validate we have channels to plot
    isempty(all_plot_channels) && throw(ArgumentError("No channels selected for plotting"))

    # Compute colors for each condition using same logic as plot_erp
    n_conditions = length(dat_subset)
    # For epochs, we have one color per condition (not per condition-channel combination)
    # But we need colors for both epochs and averages, so we'll compute for conditions
    condition_colors_result = _compute_dataset_colors(
        plot_kwargs[:color],
        n_conditions,
        1,  # n_channels = 1 for condition colors
        plot_kwargs[:colormap],
        plot_kwargs[:_color_explicitly_set]
    )
    # Convert gradient to vector of colors if needed
    if condition_colors_result isa AbstractVector
        condition_colors_list = condition_colors_result
    else
        # It's a gradient - extract colors from it
        condition_colors_list = [condition_colors_result[i] for i in 1:n_conditions]
    end

    # Info about what we're plotting
    @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(n_conditions) conditions"

    fig = Figure()
    axes = Axis[]

    # Apply theme font size early
    set_theme!(fontsize = plot_kwargs[:theme_fontsize])

    # Handle layout (similar to single EpochData method)
    if plot_kwargs[:average_channels]
        # Single plot averaging across channels
        ax = Axis(fig[1, 1])
        push!(axes, ax)

        # For each condition, average channels and plot
        for (cond_idx, dat) in enumerate(dat_subset)
            dat_avg = copy(dat)
            channel_average!(
                dat_avg;
                channel_selections = [channels(all_plot_channels)],
                output_labels = [:avg],
                reduce = false,
            )
            
            # Update plot_kwargs with condition-specific color (single value, not array)
            cond_plot_kwargs = merge(plot_kwargs, Dict(
                :color => condition_colors_list[cond_idx]
            ))
            
            _plot_epochs!(ax, dat_avg, [:avg], cond_plot_kwargs; label = nothing)

            # Optional ERP overlay
            if plot_kwargs[:plot_avg_trials]
                erp_dat = average_epochs(dat_avg)
                _plot_erp_average!(ax, erp_dat, [:avg], cond_plot_kwargs; label = label)
            end
        end
        
        # Add legend
        _add_epochs_legend!(ax, [:avg], dat_subset, plot_kwargs)

        # Set axis properties
        if length(all_plot_channels) == 1
            ax.title = string(all_plot_channels[1])
        else
            ax.title = "Avg: $(print_vector(all_plot_channels, max_length = 8, n_ends = 3))"
        end
        
        _set_axis_properties!(ax; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim], 
                           xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], yreversed = plot_kwargs[:yreversed])
        _set_axis_grid!(ax; 
                         xgrid = plot_kwargs[:xgrid], 
                         ygrid = plot_kwargs[:ygrid],
                         xminorgrid = plot_kwargs[:xminorgrid], 
                         yminorgrid = plot_kwargs[:yminorgrid])
        _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])

    else
        # Multi-channel layout (grid, topo, or single)
        # This follows the same pattern as single EpochData but plots all conditions
        # For now, delegate to the existing layout logic but modify _plot_epochs! calls
        # to handle multiple conditions
        
        # Determine layout type
        if layout == :single
            # Single plot with all channels overlaid
            ax = Axis(fig[1, 1])
            push!(axes, ax)
            
            # Plot all conditions for each channel
            for ch in all_plot_channels
                for (cond_idx, dat) in enumerate(dat_subset)
                    # Get label for this condition
                    label = isempty(plot_kwargs[:legend_labels]) ? dat.condition_name : plot_kwargs[:legend_labels][cond_idx]
                    if length(all_plot_channels) > 1
                        label *= " ($ch)"
                    end
                    
                    cond_plot_kwargs = merge(plot_kwargs, Dict(
                        :color => condition_colors_list[cond_idx]
                    ))
                    _plot_epochs!(ax, dat, [ch], cond_plot_kwargs; label = nothing)
                    
                    if plot_kwargs[:plot_avg_trials]
                        erp_dat = average_epochs(dat)
                        _plot_erp_average!(ax, erp_dat, [ch], cond_plot_kwargs; label = label)
                    end
                end
            end
            
            # Add legend
            _add_epochs_legend!(ax, all_plot_channels, dat_subset, plot_kwargs)
            
            # Set axis properties
            _set_axis_properties!(ax; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim], 
                               xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], yreversed = plot_kwargs[:yreversed])
            _set_axis_grid!(ax; 
                             xgrid = plot_kwargs[:xgrid], 
                             ygrid = plot_kwargs[:ygrid],
                             xminorgrid = plot_kwargs[:xminorgrid], 
                             yminorgrid = plot_kwargs[:yminorgrid])
            _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
            
        elseif layout == :grid || layout isa Vector{Int}
            # Grid layout - each channel gets its own subplot
            grid_dims = layout isa Vector{Int} ? (layout[1], layout[2]) : best_rect(length(all_plot_channels))
            _plot_epochs_grid_multi!(fig, axes, dat_subset, all_plot_channels, grid_dims, condition_colors_list, plot_kwargs)
            
        elseif layout == :topo
            # Topographic layout
            _plot_epochs_topo_multi!(fig, axes, dat_subset, all_plot_channels, condition_colors_list, plot_kwargs)
        end
    end

    # Setup interactivity if requested
    if plot_kwargs[:interactive]
        _setup_shared_interactivity!(fig, axes, :epochs)
        
        # Setup control panel for Vector{EpochData} (baseline + condition toggling)
        if length(dat_subset) > 1
            _setup_epochs_control_panel!(fig, dat_subset, axes, layout, all_plot_channels, condition_colors_list, plot_kwargs)
        end
    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    return fig, length(axes) == 1 ? first(axes) : axes
end

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

function _plot_epochs!(ax, dat, channels, plot_kwargs; label::Union{String,Nothing} = nothing)::Nothing
    # This function expects exactly one channel; callers pass [:avg] or [channel]
    @info "plot_epochs: $(print_vector(channels))"
    @assert length(channels) == 1 "_plot_epochs! expects a single channel"
    ch = channels[1]

    # Cache time vector and styles
    time_vec = dat.data[1][!, :time]
    # Color is passed as single value or from color cycle, not array
    trial_color = plot_kwargs[:color] isa Vector ? plot_kwargs[:color][1] : plot_kwargs[:color]
    trial_linewidth = plot_kwargs[:linewidth] isa Vector ? plot_kwargs[:linewidth][1] : plot_kwargs[:linewidth]

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

    lines!(ax, time_cat, y_cat, color = trial_color, linewidth = trial_linewidth, alpha = plot_kwargs[:trial_alpha], label = label)

    return nothing
end


function _plot_erp_average!(ax, erp_dat::ErpData, channels::Vector{Symbol}, plot_kwargs; label::Union{String,Nothing} = nothing)::Nothing
    @assert length(channels) == 1 "_plot_erp_average! expects a single channel"
    ch = channels[1]

    time_vec = erp_dat.data[!, :time]
    # Color is passed as single value or from color cycle, not array
    avg_color = plot_kwargs[:color] isa Vector ? plot_kwargs[:color][1] : plot_kwargs[:color]
    base_linewidth = plot_kwargs[:linewidth] isa Vector ? plot_kwargs[:linewidth][1] : plot_kwargs[:linewidth]
    avg_linewidth = base_linewidth * plot_kwargs[:avg_linewidth_multiplier]

    lines!(ax, time_vec, erp_dat.data[!, ch], color = avg_color, linewidth = avg_linewidth, alpha = 1.0, label = label)
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

"""
    _plot_epochs_grid_multi!(fig, axes, datasets, all_plot_channels, grid_dims, condition_colors, plot_kwargs)

Create a grid layout for plotting epochs from multiple conditions.
"""
function _plot_epochs_grid_multi!(
    fig::Figure,
    axes::Vector{Axis},
    datasets::Vector{EpochData},
    all_plot_channels::Vector{Symbol},
    grid_dims::Tuple{Int, Int},
    condition_colors::Vector,
    plot_kwargs::Dict,
)
    rows, cols = grid_dims

    # Calculate y-range if not provided (use first dataset as reference)
    ylim = plot_kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(datasets[1]; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    for (idx, channel) in enumerate(all_plot_channels)
        row = fld(idx-1, cols) + 1
        col = mod(idx-1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)
        
        # Plot all conditions for this channel
        for (cond_idx, dat) in enumerate(datasets)
            # Get label for this condition
            label = isempty(plot_kwargs[:legend_labels]) ? dat.condition_name : plot_kwargs[:legend_labels][cond_idx]
            
            cond_plot_kwargs = merge(plot_kwargs, Dict(
                :color => condition_colors[cond_idx]
            ))
            _plot_epochs!(ax, dat, [channel], cond_plot_kwargs; label = nothing)
            
            if plot_kwargs[:plot_avg_trials]
                erp_dat = average_epochs(dat)
                _plot_erp_average!(ax, erp_dat, [channel], cond_plot_kwargs; label = label)
            end
        end
        
        # Add legend for this channel
        _add_epochs_legend!(ax, [channel], datasets, plot_kwargs)

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

        ax.title = string(channel)
        _set_axis_properties!(ax; xlim = axis_kwargs[:xlim], ylim = axis_kwargs[:ylim],
                           xlabel = axis_kwargs[:xlabel], ylabel = axis_kwargs[:ylabel], yreversed = axis_kwargs[:yreversed])
        _set_axis_grid!(ax; 
                         xgrid = axis_kwargs[:xgrid], 
                         ygrid = axis_kwargs[:ygrid],
                         xminorgrid = axis_kwargs[:xminorgrid], 
                         yminorgrid = axis_kwargs[:yminorgrid])
        _set_origin_lines!(ax; add_xy_origin = axis_kwargs[:add_xy_origin])
    end

    # Link axes for synchronized zooming
    if length(axes) > 1
        Makie.linkaxes!(axes...)
    end
end

"""
    _plot_epochs_topo_multi!(fig, axes, datasets, all_plot_channels, condition_colors, plot_kwargs)

Create a topographic layout for plotting epochs from multiple conditions.
"""
function _plot_epochs_topo_multi!(
    fig::Figure,
    axes::Vector{Axis},
    datasets::Vector{EpochData},
    all_plot_channels::Vector{Symbol},
    condition_colors::Vector,
    plot_kwargs::Dict,
)
    # Use first dataset for layout (all should have same layout after subsetting)
    dat = datasets[1]
    
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
        
        # Plot all conditions for this channel
        for (cond_idx, dat_cond) in enumerate(datasets)
            # Get label for this condition
            label = isempty(plot_kwargs[:legend_labels]) ? dat_cond.condition_name : plot_kwargs[:legend_labels][cond_idx]
            
            cond_plot_kwargs = merge(plot_kwargs, Dict(
                :color => condition_colors[cond_idx]
            ))
            _plot_epochs!(ax, dat_cond, [ch], cond_plot_kwargs; label = nothing)
            
            if plot_kwargs[:plot_avg_trials]
                erp_dat = average_epochs(dat_cond)
                _plot_erp_average!(ax, erp_dat, [ch], cond_plot_kwargs; label = label)
            end
        end
        
        # Add legend for this channel
        _add_epochs_legend!(ax, [ch], datasets, plot_kwargs)

        # Suppress axis labels on all but the final axis
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
        hidespines!(ax)
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
    end
end

"""
    _add_epochs_legend!(ax::Axis, channels::Vector{Symbol}, datasets::Vector{EpochData}, kwargs::Dict)

Add legend to axis for epochs plot if conditions are met.
Similar to _add_legend! for ERPs but adapted for EpochData.
"""
function _add_epochs_legend!(ax::Axis, channels::Vector{Symbol}, datasets::Vector{EpochData}, kwargs::Dict)
    # Check if legend should be shown
    # Do not show if requested false, or single channel + single dataset
    if !kwargs[:legend] || (length(channels) == 1 && length(datasets) == 1)
        return nothing
    end
    if !isempty(kwargs[:legend_channel]) && isempty(intersect(kwargs[:legend_channel], channels))
        return nothing
    end
    
    # Extract legend parameters
    legend_label = kwargs[:legend_label]
    legend_position = kwargs[:legend_position]
    if kwargs[:legend_nbanks] === nothing
        kwargs[:legend_nbanks] = length(channels) > 10 ? cld(length(channels), 10) : 1
    end
    legend_kwargs = _extract_legend_kwargs(kwargs)

    # Add legend with position and optional label
    if legend_label != ""
        leg = axislegend(ax, legend_label; position = legend_position, legend_kwargs...)
    else
        leg = axislegend(ax; position = legend_position, legend_kwargs...)
    end
    
    return leg
end

"""
    _setup_epochs_control_panel!(fig::Figure, dat_subset::Vector{EpochData}, axes::Vector{Axis}, 
                                 layout, all_plot_channels::Vector{Symbol}, condition_colors_list, plot_kwargs::Dict)

Set up a control panel that opens when 'c' key is pressed.
Allows adjusting baseline and toggling conditions.
"""
function _setup_epochs_control_panel!(
    fig::Figure,
    dat_subset::Vector{EpochData},
    axes::Vector{Axis},
    layout,
    all_plot_channels::Vector{Symbol},
    condition_colors_list,
    plot_kwargs::Dict,
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)
    
    # State: baseline values and condition selections
    baseline_start_obs = Observable("")
    baseline_stop_obs = Observable("")
    condition_checked = [Observable(true) for _ in dat_subset]
    
    # Store textbox references
    start_input_ref = Ref{Union{Textbox,Nothing}}(nothing)
    stop_input_ref = Ref{Union{Textbox,Nothing}}(nothing)
    
    # Track last applied baseline to avoid re-applying if unchanged
    last_applied_baseline_interval = Ref{Union{IntervalTime,Nothing}}(nothing)
    
    # Update plot (re-plot everything with current settings)
    function update_plot!()
        try
            # Get baseline values from textboxes
            baseline_interval_new = nothing
            if start_input_ref[] !== nothing && stop_input_ref[] !== nothing
                start_str = start_input_ref[].stored_string[]
                stop_str = stop_input_ref[].stored_string[]
                if start_str != "" && stop_str != ""
                    try
                        baseline_interval_new = (parse(Float64, start_str), parse(Float64, stop_str))
                    catch e
                        @minimal_warning "Invalid baseline values: $e"
                        return
                    end
                end
            end
            
            # Apply baseline only if it changed
            if baseline_interval_new !== nothing
                current_baseline_interval_obj = IntervalTime(start=baseline_interval_new[1], stop=baseline_interval_new[2])
                if current_baseline_interval_obj != last_applied_baseline_interval[]
                    baseline!.(dat_subset, Ref(current_baseline_interval_obj))
                    last_applied_baseline_interval[] = current_baseline_interval_obj
                end
            end
            
            # Build condition mask
            condition_mask = [checked[] for checked in condition_checked]
            isempty([i for (i, m) in enumerate(condition_mask) if m]) && return
            
            # Clear axes
            for ax in axes
                empty!(ax)
            end
            
            # Filter datasets based on condition mask
            dat_to_plot = dat_subset[condition_mask]
            colors_to_plot = condition_colors_list[condition_mask]
            
            # Re-plot based on layout
            plot_kwargs_no_legend = merge(copy(plot_kwargs), Dict(:legend => false))
            
            if plot_kwargs[:average_channels]
                # Single plot averaging across channels
                ax = axes[1]
                for (cond_idx, dat) in enumerate(dat_to_plot)
                    dat_avg = copy(dat)
                    channel_average!(
                        dat_avg;
                        channel_selections = [channels(all_plot_channels)],
                        output_labels = [:avg],
                        reduce = false,
                    )
                    
                    label = isempty(plot_kwargs[:legend_labels]) ? dat.condition_name : 
                            (length(plot_kwargs[:legend_labels]) >= cond_idx ? plot_kwargs[:legend_labels][cond_idx] : dat.condition_name)
                    
                    cond_plot_kwargs = merge(plot_kwargs_no_legend, Dict(:color => colors_to_plot[cond_idx]))
                    _plot_epochs!(ax, dat_avg, [:avg], cond_plot_kwargs; label = nothing)
                    
                    if plot_kwargs[:plot_avg_trials]
                        erp_dat = average_epochs(dat_avg)
                        _plot_erp_average!(ax, erp_dat, [:avg], cond_plot_kwargs; label = label)
                    end
                end
                _add_epochs_legend!(ax, [:avg], dat_to_plot, plot_kwargs)
                
            elseif layout == :single
                ax = axes[1]
                for ch in all_plot_channels
                    for (cond_idx, dat) in enumerate(dat_to_plot)
                        label = isempty(plot_kwargs[:legend_labels]) ? dat.condition_name : 
                                (length(plot_kwargs[:legend_labels]) >= cond_idx ? plot_kwargs[:legend_labels][cond_idx] : dat.condition_name)
                        if length(all_plot_channels) > 1
                            label *= " ($ch)"
                        end
                        
                        cond_plot_kwargs = merge(plot_kwargs_no_legend, Dict(:color => colors_to_plot[cond_idx]))
                        _plot_epochs!(ax, dat, [ch], cond_plot_kwargs; label = nothing)
                        
                        if plot_kwargs[:plot_avg_trials]
                            erp_dat = average_epochs(dat)
                            _plot_erp_average!(ax, erp_dat, [ch], cond_plot_kwargs; label = label)
                        end
                    end
                end
                _add_epochs_legend!(ax, all_plot_channels, dat_to_plot, plot_kwargs)
                
            elseif layout == :grid || layout isa Vector{Int}
                grid_dims = layout isa Vector{Int} ? (layout[1], layout[2]) : best_rect(length(all_plot_channels))
                # Re-plot in existing axes (they're already cleared)
                rows, cols = grid_dims
                ylim = plot_kwargs[:ylim]
                if isnothing(ylim)
                    yr = ylimits(dat_to_plot[1]; channel_selection = channels(all_plot_channels))
                    ylim = (yr[1], yr[2])
                end
                
                for (idx, channel) in enumerate(all_plot_channels)
                    if idx <= length(axes)
                        ax = axes[idx]
                        for (cond_idx, dat) in enumerate(dat_to_plot)
                            label = isempty(plot_kwargs[:legend_labels]) ? dat.condition_name : 
                                    (length(plot_kwargs[:legend_labels]) >= cond_idx ? plot_kwargs[:legend_labels][cond_idx] : dat.condition_name)
                            
                            cond_plot_kwargs = merge(plot_kwargs_no_legend, Dict(:color => colors_to_plot[cond_idx]))
                            _plot_epochs!(ax, dat, [channel], cond_plot_kwargs; label = nothing)
                            
                            if plot_kwargs[:plot_avg_trials]
                                erp_dat = average_epochs(dat)
                                _plot_erp_average!(ax, erp_dat, [channel], cond_plot_kwargs; label = label)
                            end
                        end
                        _add_epochs_legend!(ax, [channel], dat_to_plot, plot_kwargs)
                        
                        # Set axis properties
                        axis_kwargs = merge(plot_kwargs, Dict(:ylim => ylim))
                        ax.title = string(channel)
                        _set_axis_properties!(ax; xlim = axis_kwargs[:xlim], ylim = axis_kwargs[:ylim],
                                           xlabel = axis_kwargs[:xlabel], ylabel = axis_kwargs[:ylabel], yreversed = axis_kwargs[:yreversed])
                        _set_axis_grid!(ax; 
                                         xgrid = axis_kwargs[:xgrid], 
                                         ygrid = axis_kwargs[:ygrid],
                                         xminorgrid = axis_kwargs[:xminorgrid], 
                                         yminorgrid = axis_kwargs[:yminorgrid])
                        _set_origin_lines!(ax; add_xy_origin = axis_kwargs[:add_xy_origin])
                    end
                end
                
            elseif layout == :topo
                # Re-plot in existing axes for topo layout
                # This is more complex - we need to re-plot each channel axis
                dat = dat_to_plot[1]  # Use first dataset for layout info
                ylim = plot_kwargs[:ylim]
                if isnothing(ylim)
                    yr = ylimits(dat; channel_selection = channels(all_plot_channels))
                    ylim = (yr[1], yr[2])
                end
                
                # Get positions from layout (same as in _plot_epochs_topo_multi!)
                x2 = dat.layout.data.x
                y2 = dat.layout.data.y
                minx, maxx = extrema(x2)
                miny, maxy = extrema(y2)
                xrange = maxx - minx
                yrange = maxy - miny
                margin = 0.1
                plot_w = Relative(0.15)
                plot_h = Relative(0.15)
                
                pos_map = Dict{Symbol,Tuple{Float64,Float64}}()
                for (lab, x, y) in zip(dat.layout.data.label, x2, y2)
                    nx = (x - minx) / xrange
                    ny = (y - miny) / yrange
                    pos_map[Symbol(lab)] = (nx, ny)
                end
                
                for (ch_idx, ch) in enumerate(all_plot_channels)
                    if ch_idx <= length(axes)
                        ax = axes[ch_idx]
                        pos = get(pos_map, ch, (0.5, 0.5))
                        
                        for (cond_idx, dat_cond) in enumerate(dat_to_plot)
                            label = isempty(plot_kwargs[:legend_labels]) ? dat_cond.condition_name : 
                                    (length(plot_kwargs[:legend_labels]) >= cond_idx ? plot_kwargs[:legend_labels][cond_idx] : dat_cond.condition_name)
                            
                            cond_plot_kwargs = merge(plot_kwargs_no_legend, Dict(:color => colors_to_plot[cond_idx]))
                            _plot_epochs!(ax, dat_cond, [ch], cond_plot_kwargs; label = nothing)
                            
                            if plot_kwargs[:plot_avg_trials]
                                erp_dat = average_epochs(dat_cond)
                                _plot_erp_average!(ax, erp_dat, [ch], cond_plot_kwargs; label = label)
                            end
                        end
                        _add_epochs_legend!(ax, [ch], dat_to_plot, plot_kwargs)
                        
                        # Set axis properties
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
                    end
                end
            end
            
            # Re-apply axis properties
            _set_axis_properties!.(axes; plot_kwargs...)
            if layout == :grid || layout isa Vector{Int}
                # Re-apply grid-specific properties (tick label hiding)
                grid_dims = layout isa Vector{Int} ? (layout[1], layout[2]) : best_rect(length(all_plot_channels))
                rows, cols = grid_dims
                for (idx, ax) in enumerate(axes)
                    if idx <= length(all_plot_channels)
                        row = fld(idx-1, cols) + 1
                        col = mod(idx-1, cols) + 1
                        if col != 1
                            ax.yticklabelsvisible = false
                        end
                        if row != rows
                            ax.xticklabelsvisible = false
                        end
                    end
                end
            end
        catch e
            @error "Error updating plot: $e" exception=(e, catch_backtrace())
        end
    end
    
    # # Helper to check if figure is still open
    # function is_figure_open(fig_ref::Ref{Union{Figure,Nothing}})
    #     if fig_ref[] === nothing
    #         return false
    #     end
    #     try
    #         _ = fig_ref[].scene
    #         return true
    #     catch
    #         fig_ref[] = nothing
    #         return false
    #     end
    # end
   
    # TODO: again this is too slow to be useful plus bugs to fix!!    

    # Keyboard handler for 'c' key
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.c
            # if !is_figure_open(control_fig)
                control_fig[] = Figure(title = "Epochs Control Panel", size = (300, 400))
                layout_panel = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)
                
                # Baseline section
                Label(layout_panel[1, 1], "Baseline Correction", fontsize = 14, font = :bold)
                baseline_layout = GridLayout(layout_panel[2, 1], tellwidth = false, colgap = 10)
                
                Label(baseline_layout[1, 1], "Start (ms):", width = 60)
                start_input = Textbox(baseline_layout[1, 2], placeholder = "e.g. -0.2", width = 100)
                start_input.stored_string[] = baseline_start_obs[]
                connect!(baseline_start_obs, start_input.stored_string)
                start_input_ref[] = start_input
                
                Label(baseline_layout[2, 1], "End (ms):", width = 60)
                stop_input = Textbox(baseline_layout[2, 2], placeholder = "e.g. 0.0", width = 100)
                stop_input.stored_string[] = baseline_stop_obs[]
                connect!(baseline_stop_obs, stop_input.stored_string)
                stop_input_ref[] = stop_input
                
                # Apply button
                apply_btn = Button(layout_panel[3, 1], label = "Apply Baseline", width = 200)
                on(apply_btn.clicks) do _
                    if start_input_ref[] !== nothing && hasproperty(start_input_ref[], :displayed_string)
                        start_input = start_input_ref[]
                        if start_input.displayed_string[] != start_input.stored_string[]
                            start_input.stored_string[] = start_input.displayed_string[]
                        end
                    end
                    if stop_input_ref[] !== nothing && hasproperty(stop_input_ref[], :displayed_string)
                        stop_input = stop_input_ref[]
                        if stop_input.displayed_string[] != stop_input.stored_string[]
                            stop_input.stored_string[] = stop_input.displayed_string[]
                        end
                    end
                    update_plot!()
                end
                
                # Conditions section
                Label(layout_panel[4, 1], "Conditions", fontsize = 14, font = :bold)
                conditions_layout = GridLayout(layout_panel[5, 1], tellwidth = false, rowgap = 5)
                
                for (i, dat) in enumerate(dat_subset)
                    cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][])
                    Label(conditions_layout[i, 2], dat.condition_name)
                    connect!(condition_checked[i], cb.checked)
                end
                
                # Auto-update on condition changes (re-plot)
                for checked in condition_checked
                    on(checked) do _
                        update_plot!()
                    end
                end
                
                display(control_fig[])
            # end
        end
    end
end

