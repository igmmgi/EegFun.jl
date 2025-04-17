function plot_erp_image(
    dat::EpochData,
    channels::Vector{Symbol};
    colorrange = nothing,
    erp_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)

    erp_default_kwargs = Dict(:plot_erp => true)
    erp_kwargs = merge(erp_default_kwargs, erp_kwargs)
    plot_erp = pop!(erp_kwargs, :plot_erp)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    data = zeros(length(dat.data), nrow(dat.data[1]))
    for epoch in eachindex(dat.data)
        data[epoch, :] = colmeans(dat.data[epoch], channels)
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
    display(fig)
    return fig, ax
end


function plot_erp_image(dat::EpochData, channel::Symbol)
    plot_erp_image(dat, [channel])
end



