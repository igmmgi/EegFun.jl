# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ERP_IMAGE_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    
    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("Epoch", "Label for y-axis"),
    
    # Layout configuration
    :dims => (nothing, "Grid dimensions as (rows, cols). If nothing, automatically determined"),
    :hidedecorations => (false, "Whether to hide axis decorations"),
    :theme_fontsize => (24, "Font size for theme"),
    
    # Image styling
    :yreversed => (false, "Whether to reverse the y-axis"),
    :colormap => (:jet, "Colormap for the image"),
    :colorrange => (nothing, "Color range for the image. If nothing, automatically determined"),
    
    # ERP overlay
    :plot_erp => (true, "Whether to plot ERP average overlay"),
    
    # Colorbar
    :plot_colorbar => (true, "Whether to show the colorbar"),
    :colorbar_width => (30, "Width of the colorbar in pixels"),
    :colorbar_label => ("μV", "Label for the colorbar"),
    
    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    
    # Interactive features
    :interactive => (true, "Whether to enable interactive features"),
)

"""
    plot_erp_image(dat::EpochData; 
                   layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
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
function plot_erp_image(dat::EpochData; 
                       layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
                       channel_selection::Function = channels(),
                       sample_selection::Function = samples(),
                       kwargs...)
    # Merge default kwargs with provided kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_IMAGE_KWARGS, kwargs)
    
    # Use subset to get the data we want to plot (same pattern as other functions)
    dat_subset = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection
    )
    
    # Get all available channels
    selected_channels = channel_labels(dat_subset)  # Gets EEG channels from layout
    extra_channels = extra_labels(dat_subset)       # Gets extra channels (EOG, etc.)
    all_plot_channels = vcat(selected_channels, extra_channels)
    
    # Validate we have channels to plot
    if isempty(all_plot_channels)
        @minimal_error("No channels available for plotting after filtering")
    end
    
    # Create figure and apply layout system
    fig = Figure()
    plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout)
    axes, channels = apply_layout!(fig, plot_layout; plot_kwargs...)
    
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
        
        # Set color range if not provided
        if isnothing(plot_kwargs[:colorrange])
            plot_kwargs[:colorrange] = extrema(data)
        end
        
        # Create heatmap
        hm = heatmap!(ax, dat_subset.data[1].time, 1:length(dat_subset.data), transpose(data), 
                      colormap = plot_kwargs[:colormap], 
                      colorrange = plot_kwargs[:colorrange])
        push!(heatmaps, hm)
        
        # Set axis properties (only for single layout or outer edges of grid)
        if plot_layout.type == :single
            # Don't show x-axis labels on heatmap if ERP trace will be shown below
            ax.xlabel = plot_kwargs[:plot_erp] ? "" : plot_kwargs[:xlabel]
            ax.ylabel = plot_kwargs[:ylabel]
            # Hide x-axis tick labels on heatmap if ERP trace will be shown below
            if plot_kwargs[:plot_erp]
                ax.xticklabelsvisible = false
            end
            # Set title showing channels (same as plot_erp)
            ax.title = length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"
        elseif plot_layout.type == :grid
            # Use the existing grid axis properties function
            row = fld(findfirst(==(channel), channels) - 1, plot_layout.cols) + 1
            col = mod(findfirst(==(channel), channels) - 1, plot_layout.cols) + 1
            _set_grid_axis_properties!(ax, plot_layout, channel, row, col, plot_layout.rows, plot_layout.cols; 
                                     xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel])
        elseif plot_layout.type == :topo
            # Set channel name as title for topo layout
            ax.title = string(channel)
        end
        
        ax.yreversed = plot_kwargs[:yreversed]
        
        # Set limits if provided, otherwise use data limits
        if !isnothing(plot_kwargs[:xlim])
            xlims!(ax, plot_kwargs[:xlim])
        else
            xlims!(ax, extrema(dat_subset.data[1].time))
        end
        
        if !isnothing(plot_kwargs[:ylim])
            ylims!(ax, plot_kwargs[:ylim])
        else
            ylims!(ax, (0.5, length(dat_subset.data) + 0.5))
        end
        
        # Add colorbar if requested (only for single layout)
        if plot_kwargs[:plot_colorbar] && plot_layout.type == :single
            Colorbar(fig[1, 2], hm, width = plot_kwargs[:colorbar_width], label = plot_kwargs[:colorbar_label])
        end
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
    
    if length(axes) > 1
        linkaxes!(axes...)
    end

    # Display plot if requested
    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    return fig, axes
end

"""
    _setup_erp_image_interactivity!(fig::Figure, axes::Vector{Axis}, heatmaps::Vector)

Set up custom interactivity for ERP images with color scale adjustment for up/down keys.
"""
function _setup_erp_image_interactivity!(fig::Figure, axes::Vector{Axis}, heatmaps::Vector)
    _setup_keyboard_tracking(fig)
    
    # Define keyboard actions for ERP images
    keyboard_actions = Dict(
        Keyboard.up => :color_more,
        Keyboard.down => :color_less,
        Keyboard.left => :x_less,
        Keyboard.right => :x_more
    )
    
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(keyboard_actions, event.key)
            action = keyboard_actions[event.key]
            _handle_erp_image_navigation!(axes, heatmaps, action)
        end
    end
end

"""
    _handle_erp_image_navigation!(axes::Vector{Axis}, heatmaps::Vector, action::Symbol)

Handle navigation actions for ERP images.
"""
function _handle_erp_image_navigation!(axes::Vector{Axis}, heatmaps::Vector, action::Symbol)
    if action in (:color_more, :color_less) # Adjust color scale
        factor = action == :color_more ? 0.8 : 1.25
        for hm in heatmaps
            hm.colorrange[] = hm.colorrange[] .* factor
        end
    elseif action in (:x_less, :x_more) # Adjust time axis
        func = action == :x_less ? xmore! : xless!
        func(first(axes))
    end
end
