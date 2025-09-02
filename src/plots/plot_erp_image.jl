# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const DEFAULT_ERP_IMAGE_KWARGS = Dict(
    :xlim => nothing,
    :ylim => nothing,
    :xlabel => "Time (S)",
    :ylabel => "Epoch",
    :dims => nothing,
    :hidedecorations => false,
    :theme_fontsize => 24,
    :yreversed => false,
    :colormap => :jet,
    :colorrange => nothing,
    :display_plot => true,
    :plot_erp => true,
    :plot_colorbar => true,
    :colorbar_width => 30,
    :xgrid => false,
    :ygrid => false,
    :xminorgrid => false,
    :yminorgrid => false,
    :interactive => true
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

# Keyword Arguments
- `xlim`: X-axis limits (default: auto-calculated)
- `ylim`: Y-axis limits (default: auto-calculated)
- `xlabel`: X-axis label (default: "Time (S)")
- `ylabel`: Y-axis label (default: "Epoch")
- `dims`: Grid dimensions (rows, cols) (default: auto-calculated)
- `hidedecorations`: Whether to hide axis decorations (default: false)
- `theme_fontsize`: Font size for theme (default: 24)
- `yreversed`: Whether to reverse Y-axis (default: false)
- `colormap`: Colormap for the image (default: :RdBu)
- `colorrange`: Color range for the image (default: auto-calculated)
- `display_plot`: Whether to display the plot (default: true)
- `plot_erp`: Whether to plot the ERP trace below (default: true)
- `plot_colorbar`: Whether to plot the colorbar (default: true)
- `colorbar_width`: Width of the colorbar (default: 30)
"""
function plot_erp_image(dat::EpochData; 
                       layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
                       channel_selection::Function = channels(),
                       sample_selection::Function = samples(),
                       kwargs...)
    # Merge default kwargs with provided kwargs
    plot_kwargs = merge(DEFAULT_ERP_IMAGE_KWARGS, Dict(kwargs))
    
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
        error("No channels available for plotting after filtering")
    end
    
    # Create figure and apply layout system
    fig = Figure()
    plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout)
    axes, channels = apply_layout!(fig, plot_layout; plot_kwargs...)
    
    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)
    
    # Disable default interactions that conflict with our custom selection
    if get(plot_kwargs, :interactive, true)
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end
    end
    
    # Plot ERP images for each channel
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
            # For grid layout, only set labels on outer edges
            # The _set_grid_axis_properties! function already handles this correctly
            # We just need to set the actual label text for the outer edges
            row = fld(findfirst(==(channel), channels) - 1, plot_layout.cols) + 1
            col = mod(findfirst(==(channel), channels) - 1, plot_layout.cols) + 1
            
            if col == 1  # Leftmost column
                ax.ylabel = plot_kwargs[:ylabel]
            end
            if row == plot_layout.rows  # Bottom row
                ax.xlabel = plot_kwargs[:xlabel]
            end
            
            # Ensure tick labels are hidden for inner plots (in case they got reset)
            if col != 1
                ax.yticklabelsvisible = false
            end
            if row != plot_layout.rows
                ax.xticklabelsvisible = false
            end
            
            # Set channel name as title
            ax.title = string(channel)
        elseif plot_layout.type == :topo
            # Set channel name as title for topo layout
            ax.title = string(channel)
        end
        
        ax.yreversed = plot_kwargs[:yreversed]
        
        # Set limits if provided
        !isnothing(plot_kwargs[:xlim]) && xlims!(ax, plot_kwargs[:xlim])
        !isnothing(plot_kwargs[:ylim]) && ylims!(ax, plot_kwargs[:ylim])
        
        # Add colorbar if requested (only for single layout)
        if plot_kwargs[:plot_colorbar] && plot_layout.type == :single
            Colorbar(fig[1, 2], hm, width = plot_kwargs[:colorbar_width])
        end
    end
    
    # Set up interactivity AFTER heatmaps to ensure rectangles are drawn on top
    if get(plot_kwargs, :interactive, true)
        # Add keyboard interactivity
        _setup_shared_interactivity!(fig, axes)
        
        # Set up selection system (rectangles created AFTER heatmaps)
        selection_state = SharedSelectionState(axes)
        _setup_unified_selection!(fig, axes, selection_state, dat_subset, plot_layout)
        
        # Add debug output to see if time selection is working
        println("DEBUG: Set up time selection for plot_erp_image with $(length(axes)) axes")
        println("DEBUG: plot_layout.type = $(plot_layout.type)")
        println("DEBUG: dat_subset type = $(typeof(dat_subset))")
        
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

    # Display plot if requested
    if plot_kwargs[:display_plot]
        display(fig)
    end

    return fig, axes
end
