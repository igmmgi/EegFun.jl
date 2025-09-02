# TODO: this still feels a bit too slow!
# TODO: proper origin ax
#
# Layout options:
# - :single (default): Single plot with all channels overlaid
# - :grid: Auto-grid layout with best_rect() dimensions
# - :topo: Topographic layout using data's layout semantics
# - [rows, cols]: Custom grid dimensions
#
# Interactive features (when interactive = true):
# - Arrow keys: Zoom in/out on X and Y axes
# - Shift + Left Click + Drag: Select time regions (works on all layouts)
# - Ctrl + Left Click + Drag: Select channels (topo/grid layouts)
# - Left Click (without Ctrl/Shift): Clear selections
# - Right Click: Print "TODO" for future functionality

# Shared interactivity functions are now included in src/eegfun.jl

# =============================================================================
# EPOCHS-SPECIFIC SELECTION STATE
# =============================================================================

# Use the shared selection state
const EpochsSelectionState = SharedSelectionState

# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================

const DEFAULT_EPOCHS_KWARGS = Dict(
    :average_channels => false,
    :xlim => nothing,
    :ylim => nothing,
    :title => nothing,
    :xlabel => "Time (S)",
    :ylabel => "mV",
    :linewidth => [1, 2],
    :color => [:grey, :red],
    :yreversed => false,
    :layout => :single,  # :single (default), :grid, or :topo
    :layout_plot_width => 0.12,
    :layout_plot_height => 0.12,
    :layout_margin => 0.02,
    :layout_show_scale => true,
    :layout_scale_position => [0.95, 0.05],
    :layout_scale_width => 0.14,
    :layout_scale_height => 0.14,
    :legend => true,
    :legend_label => "",
    :dims => nothing,
    :hidedecorations => false,
    :theme_fontsize => 24,
    :plot_avg_trials => true,                # draw ERP average overlay
    :axes_through_origin => true,
    :xgrid => false,
    :ygrid => false,
    :xminorgrid => false,
    :yminorgrid => false,
    :interactive => true,  # Enable/disable keyboard interactivity and channel selection
)

function plot_epochs(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(), 
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    layout = :single,  # :single (default), :grid, :topo, or [rows, cols]
    kwargs = Dict())::Tuple{Figure, Union{Axis, Vector{Axis}}}
    
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
    plot_kwargs = merge(copy(DEFAULT_EPOCHS_KWARGS), kwargs)

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
            # deregister_interaction!(ax, :scrollzoom)
            # deregister_interaction!(ax, :dragpan)
        end
        
        _setup_shared_interactivity!(fig, axes)
    end
    
    # Add unified selection functionality (time + channel selection)
    if plot_kwargs[:interactive]
        # Create a selection state for both time and channel selection
        selection_state = EpochsSelectionState(axes)
        
        # Set up unified selection (time selection with Shift, channel selection with Ctrl)
        if layout === :topo
            # Topographic layout - use unified selection + topo channel selection
            _setup_epochs_selection_unified!(fig, axes, selection_state, dat_subset, create_topo_layout(dat_subset.layout, all_plot_channels))
            _setup_topo_channel_selection_events!(fig, selection_state, create_topo_layout(dat_subset.layout, all_plot_channels), dat_subset, axes)
        elseif layout === :grid || typeof(layout) <: Vector{<:Integer}
            # Grid layout - use unified selection + grid channel selection
            if layout === :grid
                rows, cols = best_rect(length(all_plot_channels))
            else
                rows, cols = layout
            end
            _setup_epochs_selection_unified!(fig, axes, selection_state, dat_subset, create_grid_layout(all_plot_channels, rows, cols))
            _setup_grid_channel_selection_events!(fig, selection_state, create_grid_layout(all_plot_channels, rows, cols), dat_subset, axes)
        elseif layout === :single
            # Single plot layout - use unified selection + grid-based channel selection
            _setup_epochs_selection_unified!(fig, axes, selection_state, dat_subset, create_grid_layout(all_plot_channels, 1, 1))
            _setup_grid_channel_selection_events!(fig, selection_state, create_grid_layout(all_plot_channels, 1, 1), dat_subset, axes)
        end
    end

    # Theme adjustments
    fontsize_theme = Theme(fontsize = plot_kwargs[:theme_fontsize])
    update_theme!(fontsize_theme)

    display(fig)
    
    # Return fig and all axes (or single axis if only one)
    return fig, length(axes) == 1 ? first(axes) : axes

end

function _plot_epochs!(ax, dat, channels, kwargs)::Nothing
    # This function expects exactly one channel; callers pass [:avg] or [channel]
    println("plot_epochs! $(channels)")
    @assert length(channels) == 1 "_plot_epochs! expects a single channel"
    ch = channels[1]

    # Cache time vector and styles
    time_vec = dat.data[1][!, :time]
    trial_color = kwargs[:color][1]
    trial_linewidth = kwargs[:linewidth][1]

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


function _plot_epochs_from_erp!(ax, erp_dat::ErpData, channels::Vector{Symbol}, kwargs)::Nothing
    @assert length(channels) == 1 "_plot_epochs_from_erp! expects a single channel"
    ch = channels[1]

    time_vec = erp_dat.data[!, :time]
    avg_color = kwargs[:color][2]
    avg_linewidth = kwargs[:linewidth][2]

    lines!(ax, time_vec, erp_dat.data[!, ch], color = avg_color, linewidth = avg_linewidth)
    return nothing
end

function _plot_epochs_layout!(fig::Figure, axes::Vector{Axis}, dat::EpochData, erp_dat, all_plot_channels::Vector{Symbol}, kwargs::Dict)
    # Ensure 2D coordinates exist
    if !all(in.([:x2, :y2], Ref(propertynames(dat.layout.data))))
        polar_to_cartesian_xy!(dat.layout)
    end

    # Determine global y-lims for consistency across small axes
    ylim = kwargs[:ylim]
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
        _plot_epochs!(ax, dat, [ch], kwargs)
        erp_dat !== nothing && _plot_epochs_from_erp!(ax, erp_dat, [ch], kwargs)

        # Draw axes through origin if requested
        if get(kwargs, :axes_through_origin, false)
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
    ylim = kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(dat; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    for (idx, channel) in enumerate(all_plot_channels)
        row = fld(idx-1, cols) + 1
        col = mod(idx-1, cols) + 1
        ax = Axis(fig[row, col])
        push!(axes, ax)
        _plot_epochs!(ax, dat, [channel], kwargs)
        erp_dat !== nothing && _plot_epochs_from_erp!(ax, erp_dat, [channel], kwargs)

        # Draw axes through origin if requested
        if get(kwargs, :axes_through_origin, false)
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

# Navigation functions are now handled by _handle_shared_navigation! in shared_interactivity.jl

# Time selection helper functions are now in shared_interactivity.jl

# =============================================================================
# UNIFIED SELECTION SETUP
# =============================================================================

"""
    _setup_epochs_selection_unified!(fig::Figure, axes::Vector{Axis}, selection_state::EpochsSelectionState, data, plot_layout)

Set up unified mouse selection for epochs plots that works across all layouts.
Uses figure-level events to avoid conflicts with multiple axis handlers.
"""
function _setup_epochs_selection_unified!(fig::Figure, axes::Vector{Axis}, selection_state::EpochsSelectionState, data, plot_layout)
    # Track if Shift and Ctrl are currently pressed
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    # Use figure-level keyboard events for Shift and Ctrl tracking
    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
    end

    # Use figure-level mouse events to avoid conflicts
    on(events(fig).mousebutton) do event
        # Find which axis the mouse is over
        mouse_pos = events(fig).mouseposition[]
        active_ax = nothing
        
        for ax in axes
            if _is_mouse_in_axis(ax, mouse_pos)
                active_ax = ax
                break
            end
        end
        
        if isnothing(active_ax)
            return
        end

        mouse_x = mouseposition(active_ax)[1]

        if event.button == Mouse.left
            if event.action == Mouse.press
                if shift_pressed[] && _is_within_selection(selection_state, mouse_x)
                    _clear_shared_selection!(selection_state)
                elseif shift_pressed[]
                    _start_shared_selection!(active_ax, selection_state, mouse_x)
                end
            elseif event.action == Mouse.release && selection_state.active[]
                _finish_shared_selection!(active_ax, selection_state, mouse_x)
            end
        elseif event.button == Mouse.right && event.action == Mouse.press
            # TODO: Implement right click functionality for epochs
            @info "TODO: implement right click functionality for epochs"
        end
    end

    # Update selection rectangle while dragging using figure-level mouse position
    on(events(fig).mouseposition) do _
        if selection_state.active[]
            # Update time selection rectangle
            mouse_pos = events(fig).mouseposition[]
            active_ax = nothing
            
            for ax in axes
                if _is_mouse_in_axis(ax, mouse_pos)
                    active_ax = ax
                    break
                end
            end
            
            if !isnothing(active_ax)
                world_pos = mouseposition(active_ax)[1]
                _update_shared_selection!(active_ax, selection_state, selection_state.bounds[][1], world_pos)
            end
        end
    end
end

# =============================================================================
# MAIN INTERACTIVITY SETUP
# =============================================================================

# Keyboard interactivity is now handled by _setup_shared_interactivity!

# =============================================================================
# LAYOUT HELPER FUNCTIONS
# =============================================================================

"""
    create_topo_layout(channels, layout_data)

Create a topographic layout for the given channels.
"""
function create_topo_layout(channels::Vector{Symbol}, layout_data)
    # This is a simplified version - you might want to expand this
    return Dict(:type => :topo, :channels => channels, :layout_data => layout_data)
end

"""
    create_grid_layout(channels, rows, cols)

Create a grid layout for the given channels.
"""
function create_grid_layout(channels::Vector{Symbol}, rows::Int, cols::Int)
    return Dict(:type => :grid, :channels => channels, :rows => rows, :cols => cols)
end

# =============================================================================
# CHANNEL SELECTION EVENT FUNCTIONS
# =============================================================================

"""
    _setup_topo_channel_selection_events!(fig, selection_state, plot_layout, data, axes)

Set up channel selection events for topographic layout.
"""
function _setup_topo_channel_selection_events!(fig::Figure, selection_state::EpochsSelectionState, plot_layout, data, axes::Vector{Axis})
    # Track Ctrl key state for channel selection
    channel_selection_active = Ref(false)  # Separate state for channel selection
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    # Use figure-level keyboard events for Shift and Ctrl tracking
    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
        if key_event.action == Keyboard.release
            channel_selection_active[] = false
        end
    end   

    # Handle mouse events for channel selection
    on(events(fig).mousebutton) do event
        # Ctrl + Left Click + Drag: Select channels (existing functionality)
        if event.button == Mouse.left && ctrl_pressed[]
            if event.action == Mouse.press
                _start_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active)
            elseif event.action == Mouse.release && channel_selection_active[]
                _finish_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active, axes)
            end
        end
        
        # Left Click (without Ctrl): Clear all channel selections
        if event.button == Mouse.left && event.action == Mouse.press && !ctrl_pressed[] && !shift_pressed[]
            _clear_all_shared_channel_selections!(fig, selection_state)
        end
        
        # Right Click: Print info TODO (for future functionality)
        if event.button == Mouse.right && event.action == Mouse.press
            @info "TODO: implement right click functionality for epochs topo"
            println("TODO: implement right click functionality for epochs topo")
        end
    end
    
    # Handle mouse movement for updating channel selection rectangle
    on(events(fig).mouseposition) do _
        if channel_selection_active[]
            _update_figure_channel_selection!(fig, selection_state, plot_layout, data)
        end
    end
end

"""
    _setup_grid_channel_selection_events!(fig, selection_state, plot_layout, data, axes)

Set up channel selection events for grid layout.
"""
function _setup_grid_channel_selection_events!(fig::Figure, selection_state::EpochsSelectionState, plot_layout, data, axes::Vector{Axis})
    # Track Ctrl key state for channel selection
    ctrl_pressed = Ref(false)
    channel_selection_active = Ref(false)  # Separate state for channel selection
    shift_pressed = Ref(false)
    
    # Use figure-level keyboard events for Shift and Ctrl tracking
    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
        if key_event.action == Keyboard.release
            channel_selection_active[] = false
        end
    end

    # Handle mouse events for channel selection
    on(events(fig).mousebutton) do event
        # Ctrl + Left Click + Drag: Select channels (existing functionality)
        if event.button == Mouse.left && ctrl_pressed[]
            if event.action == Mouse.press
                _start_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active)
            elseif event.action == Mouse.release && channel_selection_active[]
                _finish_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active, axes)
            end
        end
        
        # Left Click (without Ctrl): Clear all channel selections
        if event.button == Mouse.left && event.action == Mouse.press && !ctrl_pressed[] && !shift_pressed[]
            _clear_all_shared_channel_selections!(fig, selection_state)
        end
        
        # Right Click: Print info TODO (for future functionality)
        if event.button == Mouse.right && event.action == Mouse.press
            @info "TODO: implement right click functionality for epochs grid"
            println("TODO: implement right click functionality for epochs grid")
        end
    end
    
    # Handle mouse movement for updating channel selection rectangle
    on(events(fig).mouseposition) do _
        if channel_selection_active[]
            _update_figure_channel_selection!(fig, selection_state, plot_layout, data)
        end
    end
end

