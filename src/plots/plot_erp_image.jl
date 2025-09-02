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
    :colorbar_width => 30
)

"""
    plot_erp_image(dat::EpochData; 
                   channel_selection::Function = channels(),
                   sample_selection::Function = samples(),
                   kwargs...)

Plot ERP image for specified channels and samples.

# Arguments
- `dat::EpochData`: Epoch data structure
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
    
    # Get selected channels for averaging
    all_available_channels = channel_labels(dat_subset)  # Gets EEG channels from layout
    selected_channels_mask = channel_selection(all_available_channels)
    selected_channels = all_available_channels[selected_channels_mask]
    
    # Create the ERP image data by averaging across selected channels
    data = zeros(length(dat_subset.data), nrow(dat_subset.data[1]))
    for epoch in eachindex(dat_subset.data)
        data[epoch, :] = colmeans(dat_subset.data[epoch], selected_channels)
    end
    
    # Set color range if not provided
    if isnothing(plot_kwargs[:colorrange])
        plot_kwargs[:colorrange] = extrema(data)
    end
    
    # Create figure and main heatmap
    fig = Figure()
    ax = Axis(fig[1, 1])
    hm = heatmap!(ax, dat_subset.data[1].time, 1:length(dat_subset.data), transpose(data), 
                  colormap = plot_kwargs[:colormap], 
                  colorrange = plot_kwargs[:colorrange])
    
    # Set axis properties
    ax.xlabel = plot_kwargs[:xlabel]
    ax.ylabel = plot_kwargs[:ylabel]
    ax.yreversed = plot_kwargs[:yreversed]
    
    # Set limits if provided
    !isnothing(plot_kwargs[:xlim]) && xlims!(ax, plot_kwargs[:xlim])
    !isnothing(plot_kwargs[:ylim]) && ylims!(ax, plot_kwargs[:ylim])
    
    # Add colorbar if requested
    if plot_kwargs[:plot_colorbar]
        Colorbar(fig[1, 2], hm, width = plot_kwargs[:colorbar_width])
    end

    # Add ERP trace if requested
    if plot_kwargs[:plot_erp]
        # Convert epoch data to ERP data by averaging across epochs
        erp_data = average_epochs(dat_subset)
        
        # Create ERP axis in our main figure
        ax_erp = Axis(fig[2, 1])
        
        # Get the selected channels for the ERP
        all_available_channels = channel_labels(erp_data)
        selected_channels_mask = channel_selection(all_available_channels)
        selected_channels = all_available_channels[selected_channels_mask]
        
        # Use the same ERP plotting function as plot_erp
        _plot_erp!(ax_erp, [erp_data], selected_channels; 
            xlim = plot_kwargs[:xlim],
            xlabel = plot_kwargs[:xlabel],
            ylabel = "Amplitude (μV)",
            color = :black,
            linewidth = 1
        )
    end

    # Display plot if requested
    if plot_kwargs[:display_plot]
        display(fig)
    end

    return fig, ax
end
