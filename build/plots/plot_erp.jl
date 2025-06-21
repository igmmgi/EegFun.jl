# #################################################################
# plot_erp: ERP Data (Single Condition; Single Channel or Average of multiple channels)
function plot_erp!(fig, ax, dat::ErpData, channels::Function = channels(); kwargs = Dict())
    # Get all available channels (layout + additional)
    all_available_channels = _get_available_channels(dat)
    selected_channels = channels(all_available_channels)

    # For EEG channels, validate against layout
    eeg_channels = intersect(selected_channels, dat.layout.label)
    invalid_eeg_channels = setdiff(eeg_channels, dat.layout.label)
    !isempty(invalid_eeg_channels) && throw(ArgumentError("Invalid EEG channels: $(join(invalid_eeg_channels, ", "))"))
    
    # Additional channels (not in layout) are allowed
    additional_channels = setdiff(selected_channels, dat.layout.label)
    if !isempty(additional_channels)
        @info "plot_erp!: Including additional channels not in layout: $(_print_vector(additional_channels))"
    end

    default_kwargs = Dict(
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => 2,
        :color => :black,
        :linestyle => :solid,
        :colormap => :jet,
        :yreversed => false,
        :add_topoplot => true,
        :topoplot_fig => 1,
        :topoplot_size => 0.2,
        :topoplot_position => (0, 0),
        :topoplot_color => :black,
        :average_channels => false,
        :legend => true,
        :legend_label => "",
    )
    kwargs = merge(default_kwargs, kwargs)

    # plot
    if kwargs[:average_channels]
        colors = kwargs[:color]
        lines!(
            ax,
            dat.data[!, :time],
            colmeans(dat.data, selected_channels),
            color = kwargs[:color],
            linewidth = kwargs[:linewidth],
            linestyle = kwargs[:linestyle],
            label = kwargs[:legend_label],
        )
    else
        colors = Makie.cgrad(kwargs[:colormap], length(selected_channels), categorical = true)
        for (idx, channel) in enumerate(selected_channels)
            lines!(
                ax,
                dat.data[!, :time],
                dat.data[!, channel],
                color = colors[idx],
                linewidth = kwargs[:linewidth],
                linestyle = kwargs[:linestyle],
                label = string(kwargs[:legend_label], " ", channel),
            )
        end
    end

    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    if isnothing(kwargs[:title])
        ax.title = "$(_print_vector(selected_channels))"
    else
        ax.title = kwargs[:title]
    end
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]

    # TODO: improve options for legend position
    if kwargs[:legend]
        axislegend(ax, framevisible = false, position = :lt)
    end


    if kwargs[:add_topoplot]
        topo_ax = Axis(
            fig[1, 1],
            width = Relative(kwargs[:topoplot_size]),
            height = Relative(kwargs[:topoplot_size]),
            halign = kwargs[:topoplot_position][1],
            valign = kwargs[:topoplot_position][2],
        )
        layout = filter(row -> row.label in selected_channels, dat.layout)
        if kwargs[:average_channels]
            plot_layout_2d!(
                fig,
                topo_ax,
                layout,
                point_kwargs = Dict(:color => kwargs[:topoplot_color], :markersize => 18),
            )
        else
            plot_layout_2d!(
                fig,
                topo_ax,
                layout,
                point_kwargs = Dict(:colormap => colors, :color => 1:length(selected_channels), :markersize => 18),
            )
        end
    end


    # # x/y limits
    # isnothing(kwargs[:xlim]) && (kwargs[:xlim] = data_limits_x(dat.data))
    # isnothing(kwargs[:ylim]) && (kwargs[:ylim] = data_limits_y(dat.data, dat.layout.label))
    # lines!(ax, dat.data[!, :time], zeros(nrow(dat.data)), color = :black)
    # vlines!(ax, [0], color = :black)
    # hlines!(ax, [0], color = :black)
    # xlims!(ax, kwargs[:xlim])
    # ylims!(ax, kwargs[:ylim])
    # hidedecorations!(ax)
    # hidespines!(ax, :t, :r, :l, :b)

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    display(fig)
    return fig, ax

end

"""
    plot_erp(dat::ErpData, channels::Function = channels(); kwargs = Dict())

Plot ERP data for specified channels.

# Arguments
- `dat::ErpData`: ERP data structure
- `channels::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `kwargs`: Additional plotting options

# Returns
A tuple of (figure, axis) objects.

# Examples

## Basic Usage
```julia
# Plot all channels
fig, ax = plot_erp(dat)

# Plot specific channels
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2, :F3, :F4]))

# Plot single channel
fig, ax = plot_erp(dat, channels(:Fp1))
```

## Channel Filtering
```julia
# Exclude reference channels
fig, ax = plot_erp(dat, channels_not([:M1, :M2]))

# Plot frontal channels only (by index)
fig, ax = plot_erp(dat, channels(1:10))

# Plot channels starting with "F" (frontal)
fig, ax = plot_erp(dat, x -> startswith.(string.(x), "F"))

# Plot parietal channels
fig, ax = plot_erp(dat, x -> startswith.(string.(x), "P"))

# Mix of layout and additional channels
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Plotting Options
```julia
# Average multiple channels
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2, :F3, :F4]), 
    average_channels = true,
    legend_label = "Frontal Average"
)

# Custom styling
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2]), 
    color = :red,
    linewidth = 3,
    linestyle = :dash,
    title = "Frontal ERPs"
)

# Custom axis limits
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2]), 
    xlim = (-0.2, 0.8),
    ylim = (-10, 10)
)

# Without topoplot
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2]), 
    add_topoplot = false
)

# Custom topoplot position
fig, ax = plot_erp(dat, channels([:Fp1, :Fp2]), 
    topoplot_position = (0.8, 0.8),
    topoplot_size = 0.3
)
```

## Comparison Plots
```julia
# Compare original and cleaned data
fig, ax1, ax2 = plot_erp(dat_original, dat_cleaned, channels([:Fp1, :Fp2]))

# Compare multiple conditions
fig, ax = plot_erp([dat_condition1, dat_condition2], channels([:Fp1, :Fp2]))
```

## Complex Channel Selection
```julia
# All channels except reference and EOG
fig, ax = plot_erp(dat, channels_not([:M1, :M2, :vEOG, :hEOG]))

# Frontal and parietal channels
fig, ax = plot_erp(dat, x -> 
    startswith.(string.(x), "F") .| startswith.(string.(x), "P")
)

# Channels with specific indices
fig, ax = plot_erp(dat, channels([1, 5, 10, 15, 20]))
```
"""
function plot_erp(dat::ErpData, channels::Function = channels(); kwargs = Dict())
    fig = Figure()
    ax = Axis(fig[1, 1])
    fig, ax = plot_erp!(fig, ax, dat, channels; kwargs = kwargs)
    return fig, ax
end

# Backward compatibility - keep the old method for existing code
function plot_erp(dat::ErpData, channels::Vector{Symbol}; kwargs = Dict())
    channel_predicate = x -> x .∈ Ref(channels)
    return plot_erp(dat, channel_predicate; kwargs = kwargs)
end

function plot_erp(dat::ErpData, channels::Symbol; kwargs...)
    plot_erp(dat, [channels]; kwargs...)
end

function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData, channels; kwargs = Dict())
    kwargs = merge(kwargs, kwargs)
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax1.xlabelvisible = false
    ax1.xticklabelsvisible = false
    plot_erp!(fig, ax1, dat_orig, channels; kwargs = kwargs)
    ax2 = Axis(fig[2, 1])
    kwargs[:topoplot_fig] = 2
    plot_erp!(fig, ax2, dat_cleaned, channels; kwargs = kwargs)
    ax2.title = ""
    linkaxes!(ax1, ax2)
    display(fig)
    return fig, ax1, ax2
end

function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData; kwargs...)
    plot_erp(dat_orig, dat_cleaned, dat_orig.layout.label; kwargs...)
end



function plot_erp(dat_orig::Vector{ErpData}, channels; kwargs = Dict{Symbol,Any}())
    fig = Figure()
    ax1 = Axis(fig[1, 1])

    average_channels = get(kwargs, :average_channels, false)
    colormap = get(kwargs, :colormap, :jet)

    linestyles = [:solid, :dot, :dash, :dashdot, :dashdotdot]
    colors = Makie.cgrad(colormap, length(dat_orig), categorical = true)

    # Convert kwargs to Dict{Symbol,Any} if it's not already
    kwargs = Dict{Symbol,Any}(kwargs)

    for (idx, dat) in enumerate(dat_orig)
        plot_kwargs = copy(kwargs)
        plot_kwargs[:legend_label] = string("Cond: ", idx)
        if length(dat_orig) > 1 && !average_channels
            plot_kwargs[:linestyle] = linestyles[(idx-1)%length(linestyles)+1]
        elseif length(dat_orig) > 1 && average_channels
            plot_kwargs[:color] = colors[idx]
        end
        plot_erp!(fig, ax1, dat, channels; kwargs = plot_kwargs)
    end
    return fig, ax1
end
