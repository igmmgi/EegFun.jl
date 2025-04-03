# #################################################################
# plot_grid_rect: 
function plot_grid_rect(dat::ErpData; channels = nothing, kwargs = Dict())

    isnothing(channels) && (channels = dat.layout.label)

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
        dim1 = ceil(Int, sqrt(length(channels)))
        dim2 = ceil(Int, length(channels) ./ dim1)
        kwargs[:dims] = [dim1, dim2]
    end

    # x/y limits
    isnothing(kwargs[:xlim]) && (kwargs[:xlim] = data_limits_x(dat.data))
    isnothing(kwargs[:ylim]) && (kwargs[:ylim] = data_limits_y(dat.data, dat.layout.label))

    count = 1
    fig = Figure()
    for dim1 = 1:kwargs[:dims][1]
        for dim2 = 1:kwargs[:dims][2]
            # ax = Axis(fig[dim1, dim2], width = 200, height = 150)
            ax = Axis(fig[dim1, dim2])
            lines!(ax, dat.data[!, :time], dat.data[!, channels[count]])
            vlines!(ax, [0], color = :black)
            hlines!(ax, [0], color = :black)
            ax.title = "$(channels[count])"
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
            if count > length(channels)
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

plot_grid_rect(erp, channels = [:F8, :Fp1, :Fp2, :Fpz, :Fz, :F3, :F4, :F7, :F8])