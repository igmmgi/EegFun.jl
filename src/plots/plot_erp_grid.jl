# #################################################################
# plot_grid_rect: Plot ERP data in a grid layout
"""
    plot_grid_rect(dat::ErpData; channel_predicate::Function = channels, kwargs = Dict())

Plot ERP data in a grid layout for specified channels.

# Arguments
- `dat::ErpData`: ERP data structure
- `channel_predicate::Function`: Function that returns boolean vector for channel filtering (default: channels - all channels)
- `kwargs`: Additional keyword arguments for customization

# Examples
```julia
# Plot all channels in grid
plot_grid_rect(dat)

# Plot specific channels in grid
plot_grid_rect(dat, channels([:Fp1, :Fp2, :F3, :F4]))

# Exclude reference channels
plot_grid_rect(dat, channels_not([:M1, :M2]))

# Plot frontal channels only
plot_grid_rect(dat, channels(1:10))

# Custom predicate
plot_grid_rect(dat, x -> startswith.(string.(x), "F"))
```

# Keyword Arguments
- `xlim`: X-axis limits (default: auto-calculated)
- `ylim`: Y-axis limits (default: auto-calculated)
- `xlabel`: X-axis label (default: "Time (S)")
- `ylabel`: Y-axis label (default: "mV")
- `dims`: Grid dimensions [rows, cols] (default: auto-calculated)
- `hidedecorations`: Whether to hide axis decorations (default: false)
- `theme_fontsize`: Font size for theme (default: 24)
- `yreversed`: Whether to reverse Y-axis (default: false)
"""
function plot_grid_rect(dat::ErpData; channel_predicate::Function = channels, kwargs = Dict())
    # Get the channels using the predicate
    selected_channels = channel_predicate(dat.layout.label)

    default_kwargs = Dict{Symbol,Any}(
        :xlim => nothing,
        :ylim => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :dims => nothing,
        :hidedecorations => false,
        :theme_fontsize => 24,
        :yreversed => false,
    )
    kwargs = merge(default_kwargs, kwargs)
    if isnothing(kwargs[:dims])
        dim1 = ceil(Int, sqrt(length(selected_channels)))
        dim2 = ceil(Int, length(selected_channels) ./ dim1)
        kwargs[:dims] = [dim1, dim2]
    end

    # x/y limits
    isnothing(kwargs[:xlim]) && (kwargs[:xlim] = data_limits_x(dat.data))
    isnothing(kwargs[:ylim]) && (kwargs[:ylim] = data_limits_y(dat.data, selected_channels))

    count = 1
    fig = Figure()
    for dim1 = 1:kwargs[:dims][1]
        for dim2 = 1:kwargs[:dims][2]
            # ax = Axis(fig[dim1, dim2], width = 200, height = 150)
            ax = Axis(fig[dim1, dim2])
            lines!(ax, dat.data[!, :time], dat.data[!, selected_channels[count]])
            vlines!(ax, [0], color = :black)
            hlines!(ax, [0], color = :black)
            ax.title = "$(selected_channels[count])"
            if kwargs[:hidedecorations]
                hidedecorations!(ax)
            end
            xlims!(ax, kwargs[:xlim])
            ylims!(ax, kwargs[:ylim])
            if count == 3
                ax.xlabel = kwargs[:xlabel]
                ax.ylabel = kwargs[:ylabel]
            end
            ax.yreversed = kwargs[:yreversed]
            count += 1
            if count > length(selected_channels)
                break
            end
            # colsize!(fig.layout, dim2, Relative(0.1))
        end
        # rowsize!(fig.layout, dim1, Relative(0.1))
    end
    # colgap!(fig.layout, 150)
    # rowgap!(fig.layout, 150)
    linkaxes!(filter(x -> x isa Axis, fig.content)...)
    # plot theme adjustments
    fontsize_theme = Theme(fontsize = kwargs[:theme_fontsize])
    update_theme!(fontsize_theme)
    display(fig)
    return fig, ax
end

# Backward compatibility - keep the old method for existing code
function plot_grid_rect(dat::ErpData; channels = nothing, kwargs = Dict())
    if isnothing(channels)
        return plot_grid_rect(dat, channel_predicate = channels; kwargs = kwargs)
    else
        channel_predicate = x -> x .∈ Ref(channels)
        return plot_grid_rect(dat, channel_predicate = channel_predicate; kwargs = kwargs)
    end
end

