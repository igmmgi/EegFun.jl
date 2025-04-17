# #################################################################
# plot_erp: ERP Data (Single Condition; Single Channel or Average of multiple channels)
function plot_erp(fig, ax, dat::ErpData, channels::Vector{Symbol}; kwargs = Dict())

    default_kwargs = Dict(
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => 2,
        :color => :black,
        :colormap => :jet,
        :yreversed => false,
        :add_topoplot => true,
        :topoplot_fig => 1,
        :average_channels => false,
    )
    kwargs = merge(default_kwargs, kwargs)

    # plot
    if kwargs[:average_channels]
        println("averaging")
        colors = kwargs[:color]
        lines!(
            ax,
            dat.data[!, :time],
            colmeans(dat.data, channels),
            color = kwargs[:color],
            linewidth = kwargs[:linewidth],
        )
    else
        colors = Makie.cgrad(kwargs[:colormap], length(channels), categorical = true)
        for (idx, channel) in enumerate(channels)
            lines!(ax, dat.data[!, :time], dat.data[!, channel], color = colors[idx], linewidth = kwargs[:linewidth])
        end
    end

    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    if isnothing(kwargs[:title])
        ax.title = "$(print_vector_(channels))"
    else
        ax.title = kwargs[:title]
    end
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]

    if kwargs[:add_topoplot]
        # just put in top left
        topo_ax =
            Axis(fig[kwargs[:topoplot_fig], 1], width = Relative(0.2), height = Relative(0.2), halign = 0, valign = 0)
        layout = filter(row -> row.label in channels, dat.layout)
        if kwargs[:average_channels]
            plot_layout_2d!(fig, topo_ax, layout, point_kwargs = Dict(:color => kwargs[:color], :markersize => 18))
        else
            plot_layout_2d!(
                fig,
                topo_ax,
                layout,
                point_kwargs = Dict(:colormap => colors, :color => 1:length(channels), :markersize => 18),
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

function plot_erp(dat::ErpData, channels::Vector{Symbol}; kwargs = Dict())
    fig = Figure()
    ax = Axis(fig[1, 1])
    fig, ax = plot_erp(fig, ax, dat, channels; kwargs = kwargs)
    return fig, ax
end

function plot_erp(dat::ErpData, channels::Symbol; kwargs...)
    plot_erp(dat, [channels]; kwargs...)
end

function plot_erp(dat::ErpData; kwargs...)
    plot_erp(dat, dat.layout.label; kwargs...)
end



function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData, channels; kwargs = Dict())
    kwargs = merge(kwargs, kwargs)
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    plot_erp(fig, ax1, dat_orig, channels; kwargs = kwargs)
    ax2 = Axis(fig[2, 1])
    kwargs[:topoplot_fig] = 2
    plot_erp(fig, ax2, dat_cleaned, channels; kwargs = kwargs)
    linkaxes!(ax1, ax2)
    display(fig)
    return fig, ax1, ax2
end

function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData; kwargs...)
    plot_erp(dat_orig, dat_cleaned, dat_orig.layout.label; kwargs...)
end



# # Basic Tests
# # TODO: Implement proper tests
# layout = read_layout("./layouts/biosemi72.csv");
# dat = read_bdf("../Flank_C_3.bdf");
# dat = create_eeg_dataframe(dat, layout);
# filter_data!(dat, "hp", "iir", 1, order = 1)

# # Epoch Data
# epoch = extract_epochs(dat, 1, 1, -2, 4)

# # ERP Data
# erp = average_epochs(epoch)
# fig, ax = plot_erp(erp, [:Fp1, :Fp2])
# plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:average_channels => true))


# fig, ax = plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:add_topoplot => false))
# topo_ax = Axis(fig[1,1], 
#                width=Relative(0.2),
#                height=Relative(0.2),
#                halign=0.5,
#                valign=0.5)
# layout = filter(row -> row.label in [:Fp1, :Fp2], erp.layout)
# plot_layout_2d!(fig, topo_ax, layout, 
#                 point_kwargs=Dict(:colormap => :jet, 
#                                  :color => 1:2, 
#                                  :markersize => 18))


plot_erp(erp, erp)
