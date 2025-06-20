"""
    plot_erp_image(dat::EpochData, channel_predicate::Function = channels; 
                   colorrange = nothing, erp_kwargs = Dict(), colorbar_kwargs = Dict(), 
                   display_plot::Bool = true)

Plot ERP image showing trial-by-trial activity for specified channels.

# Arguments
- `dat::EpochData`: EpochData structure containing epoched EEG data
- `channel_predicate::Function`: Function that returns boolean vector for channel filtering (default: channels - all channels)
- `colorrange`: Color range for the heatmap (default: auto-calculated)
- `erp_kwargs`: Keyword arguments for ERP plotting (default: Dict())
- `colorbar_kwargs`: Keyword arguments for colorbar (default: Dict())
- `display_plot::Bool`: Whether to display the plot (default: true)

# Examples
```julia
# Plot all channels
plot_erp_image(dat)

# Plot specific channels
plot_erp_image(dat, channels([:Fp1, :Fp2]))

# Exclude reference channels
plot_erp_image(dat, channels_not([:M1, :M2]))

# Plot frontal channels only
plot_erp_image(dat, channels(1:10))

# Custom predicate
plot_erp_image(dat, x -> startswith.(string.(x), "F"))
```

# Returns
- `fig::Figure`: The Figure object containing the plot
- `ax::Axis`: The Axis object containing the plot
"""
function plot_erp_image(
    dat::EpochData,
    channel_predicate::Function = channels;
    colorrange = nothing,
    erp_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    display_plot::Bool = true,
)
    # Get the channels using the predicate
    selected_channels = channel_predicate(dat.layout.label)

    erp_default_kwargs = Dict(:plot_erp => true)
    erp_kwargs = merge(erp_default_kwargs, erp_kwargs)
    plot_erp = pop!(erp_kwargs, :plot_erp)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    data = zeros(length(dat.data), nrow(dat.data[1]))
    for epoch in eachindex(dat.data)
        data[epoch, :] = colmeans(dat.data[epoch], selected_channels)
    end
    if isnothing(colorrange)
        colorrange = extrema(data)
    end
    fig = Figure()
    ax = Axis(fig[1, 1])
    hm = heatmap!(ax, dat.data[1].time, 1:length(dat.data), transpose(data), colorrange = colorrange)
    xlims!(ax, (-0.5, 2))
    ax.xlabel = "Time (ms)"
    ax.ylabel = "Epoch"
    if plot_colorbar
        Colorbar(fig[1, 2], hm; colorbar_kwargs...)
    end

    if plot_erp
        ax = Axis(fig[2, 1])
        lines!(ax, dat.data[1].time, colmeans(data))
        xlims!(ax, (-0.5, 2))
    end

    if display_plot
        display(fig)
    end

    return fig, ax

end

# Backward compatibility - keep the old method for existing code
function plot_erp_image(
    dat::EpochData,
    channels::Vector{Symbol};
    colorrange = nothing,
    erp_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    display_plot::Bool = true,
)
    channel_predicate = x -> x .∈ Ref(channels)
    return plot_erp_image(dat, channel_predicate; 
                        colorrange = colorrange, 
                        erp_kwargs = erp_kwargs, 
                        colorbar_kwargs = colorbar_kwargs, 
                        display_plot = display_plot)
end

function plot_erp_image(dat::EpochData, channel::Symbol; kwargs...)
    plot_erp_image(dat, [channel]; kwargs...)
end



