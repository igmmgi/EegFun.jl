# #################################################################
# plot_erp: ERP Data (Single Condition; Single Channel or Average of multiple channels)
function plot_erp!(
    fig,
    ax,
    dat::ErpData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    kwargs = Dict(),
)
    # Subset data first (consistent with plot_epochs)
    dat_subset = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = true,
    )

    # Channels present after subsetting
    selected_channels = channel_labels(dat_subset)
    extra_channels = extra_labels(dat_subset)
    all_plot_channels = vcat(selected_channels, extra_channels)

    isempty(all_plot_channels) && throw(ArgumentError("No channels selected for plotting"))

    # Defaults (align with plot_epochs where applicable)
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
        :average_channels => false,
        :legend => true,
        :legend_label => "",
    )
    kwargs = merge(default_kwargs, kwargs)

    # Plot
    if kwargs[:average_channels]
        # Compute mean across selected channels once
        lines!(
            ax,
            dat_subset.data[!, :time],
            colmeans(dat_subset.data, all_plot_channels),
            color = kwargs[:color],
            linewidth = kwargs[:linewidth],
            linestyle = kwargs[:linestyle],
            label = kwargs[:legend_label],
        )
    else
        colors = Makie.cgrad(kwargs[:colormap], length(all_plot_channels), categorical = true)
        for (idx, channel) in enumerate(all_plot_channels)
            lines!(
                ax,
                dat_subset.data[!, :time],
                dat_subset.data[!, channel],
                color = colors[idx],
                linewidth = kwargs[:linewidth],
                linestyle = kwargs[:linestyle],
                label = string(kwargs[:legend_label], " ", channel),
            )
        end
    end

    # Axis properties
    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    ax.title = isnothing(kwargs[:title]) ? "$(print_vector(all_plot_channels))" : kwargs[:title]
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]

    if kwargs[:legend]
        axislegend(ax, framevisible = false, position = :lt)
    end

    return fig, ax
end

"""
    plot_erp(dat::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), kwargs = Dict())

Create a figure/axis and plot ERP data using the same selection pattern as plot_epochs.
"""
function plot_erp(dat::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), kwargs = Dict())
    fig = Figure()
    ax = Axis(fig[1, 1])
    fig, ax = plot_erp!(fig, ax, dat; channel_selection = channel_selection, sample_selection = sample_selection, kwargs = kwargs)
    return fig, ax
end

"""
    plot_erp(dat_orig::ErpData, dat_cleaned::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), kwargs = Dict())

Plot two ERP datasets on linked axes for comparison.
"""
function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), kwargs = Dict())
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax1.xlabelvisible = false
    ax1.xticklabelsvisible = false
    plot_erp!(fig, ax1, dat_orig; channel_selection = channel_selection, sample_selection = sample_selection, kwargs = kwargs)

    ax2 = Axis(fig[2, 1])
    plot_erp!(fig, ax2, dat_cleaned; channel_selection = channel_selection, sample_selection = sample_selection, kwargs = kwargs)
    ax2.title = ""

    linkaxes!(ax1, ax2)
    display(fig)
    return fig, ax1, ax2
end

"""
    plot_erp(datasets::Vector{ErpData}; channel_selection::Function = channels(), sample_selection::Function = samples(), kwargs = Dict{Symbol,Any}())

Plot multiple ERP datasets on the same axis (e.g., conditions).
"""
function plot_erp(datasets::Vector{ErpData}; channel_selection::Function = channels(), sample_selection::Function = samples(), kwargs = Dict{Symbol,Any}())
    fig = Figure()
    ax1 = Axis(fig[1, 1])

    average_channels = get(kwargs, :average_channels, false)
    colormap = get(kwargs, :colormap, :jet)

    linestyles = [:solid, :dot, :dash, :dashdot, :dashdotdot]
    colors = Makie.cgrad(colormap, length(datasets), categorical = true)

    for (idx, dat) in enumerate(datasets)
        plot_kwargs = copy(kwargs)
        plot_kwargs[:legend_label] = string("Cond: ", idx)
        if length(datasets) > 1 && !average_channels
            plot_kwargs[:linestyle] = linestyles[(idx-1)%length(linestyles)+1]
        elseif length(datasets) > 1 && average_channels
            plot_kwargs[:color] = colors[idx]
        end
        plot_erp!(fig, ax1, dat; channel_selection = channel_selection, sample_selection = sample_selection, kwargs = plot_kwargs)
    end
    return fig, ax1
end
