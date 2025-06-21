"""
    plot_erp_image(dat::EpochData, channels::Function = channels();
                   kwargs = Dict())

Plot ERP image for specified channels.

# Arguments
- `dat::EpochData`: Epoch data structure
- `channels::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)

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

# Keyword Arguments
- `xlim`: X-axis limits (default: auto-calculated)
- `ylim`: Y-axis limits (default: auto-calculated)
- `xlabel`: X-axis label (default: "Time (S)")
- `ylabel`: Y-axis label (default: "Epoch")
- `dims`: Grid dimensions [rows, cols] (default: auto-calculated)
- `hidedecorations`: Whether to hide axis decorations (default: false)
- `theme_fontsize`: Font size for theme (default: 24)
- `yreversed`: Whether to reverse Y-axis (default: false)
- `colormap`: Colormap for the image (default: :RdBu)
- `colorrange`: Color range for the image (default: auto-calculated)
"""
function plot_erp_image(dat::EpochData, 
                       channels::Function = channels();
                       kwargs = Dict())
    selected_channels = channels(dat.layout.label)

    erp_default_kwargs = Dict(:plot_erp => true)
    erp_kwargs = merge(erp_default_kwargs, kwargs)
    plot_erp = pop!(erp_kwargs, :plot_erp)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    data = zeros(length(dat.data), nrow(dat.data[1]))
    for epoch in eachindex(dat.data)
        data[epoch, :] = colmeans(dat.data[epoch], selected_channels)
    end
    if isnothing(kwargs[:colorrange])
        kwargs[:colorrange] = extrema(data)
    end
    fig = Figure()
    ax = Axis(fig[1, 1])
    hm = heatmap!(ax, dat.data[1].time, 1:length(dat.data), transpose(data), kwargs...)
    xlims!(ax, (-0.5, 2))
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    if plot_colorbar
        Colorbar(fig[1, 2], hm; colorbar_kwargs...)
    end

    if plot_erp
        ax = Axis(fig[2, 1])
        lines!(ax, dat.data[1].time, colmeans(data))
        xlims!(ax, (-0.5, 2))
    end

    if kwargs[:display_plot]
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
    plot_erp_image(dat, x -> x .== channel; kwargs...)
end



