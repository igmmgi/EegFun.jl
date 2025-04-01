# TODO: butterfly plot/global field power
# TODO: spline interpolation for topoplots?

# #################################################################
# plot_epochs: Epoched Data (Single Condition; Single Channel or Average of multiple channels)
average_channels = true, function plot_epochs(dat::EpochData, channels::Vector{Symbol}; average_channels = true, kwargs = Dict())

    default_kwargs = Dict(
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => [1, 3],
        :color => [:grey, :black],
        :yreversed => false,
    )

    kwargs = merge(default_kwargs, kwargs)

    fig = Figure()
    ax = Axis(fig[1, 1])

    # plot each trial (average acrossed electordes if >1) and overall average
    avg_data = zeros(nrow(dat.data[1]))
    for trial in eachindex(dat.data)
        trial_data = colmeans(dat.data[trial], channels)
        avg_data .+= trial_data
        lines!(dat.data[trial][!, :time], trial_data, color = kwargs[:color][1], linewidth = kwargs[:linewidth][1])
    end
    avg_data ./= length(dat.data)
    lines!(dat.data[1][!, :time], avg_data, color = kwargs[:color][2], linewidth = kwargs[:linewidth][2])

    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    if isnothing(kwargs[:title])
        ax.title = length(channels) == 1 ? "Electrode: $(channels[1])" : "Electrodes Avg: $(""*join(channels,",")*"")"
    else
        ax.title = kwargs[:title]
    end
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    display(fig)
    return fig, ax

end

plot_epochs(dat::EpochData, channels::Symbol; kwargs...) = plot_epochs(dat::EpochData, [channels]; kwargs...)


# Basic Tests
# TODO: Implement proper tests
layout = read_layout("./layouts/biosemi72.csv");
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");
dat = create_eeg_dataframe(dat, layout);
filter_data!(dat, "hp", "iir", 1, order=1)

# Epoch Data
epoch = extract_epochs(dat, 1, 1, -2, 4)
plot_epochs(epoch, :Fp1) 
plot_epochs(epoch, [:Fp1, :Fp2]) 






# #################################################################
# plot_epochs: Epoched Data (Single Condition; Single Channel or Average of multiple channels)
function plot_epochs(dat::EpochData, channels::Vector{Symbol}; kwargs = Dict())

    default_kwargs = Dict(
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => [1, 3],
        :color => [:grey, :black],
        :yreversed => false,
    )

    kwargs = merge(default_kwargs, kwargs)

    fig = Figure()
    ax = Axis(fig[1, 1])

    # plot each trial (average acrossed electordes if >1) and overall average
    avg_data = zeros(nrow(dat.data[1]))
    for trial in eachindex(dat.data)
        trial_data = colmeans(dat.data[trial], channels)
        avg_data .+= trial_data
        lines!(dat.data[trial][!, :time], trial_data, color = kwargs[:color][1], linewidth = kwargs[:linewidth][1])
    end
    avg_data ./= length(dat.data)
    lines!(dat.data[1][!, :time], avg_data, color = kwargs[:color][2], linewidth = kwargs[:linewidth][2])

    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    if isnothing(kwargs[:title])
        ax.title = length(channels) == 1 ? "Electrode: $(channels[1])" : "Electrodes Avg: $(""*join(channels,",")*"")"
    else
        ax.title = kwargs[:title]
    end
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    display(fig)
    return fig, ax

end

plot_epochs(dat::EpochData, channels::Symbol; kwargs...) = plot_epochs(dat::EpochData, [channels]; kwargs...)


# #################################################################
# plot_erp: ERP Data (Single Condition; Single Channel or Average of multiple channels)
function plot_erp(fig, ax, dat::ErpData, channels::Vector{Symbol}; average_channels = false, kwargs = Dict())

    default_kwargs = Dict(
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => 4,
        :color => :black,
        :colormap => :viridis,
        :yreversed => false,
        :add_topoplot => true,
        :topoplot_fig => 1,
    )
    kwargs = merge(default_kwargs, kwargs)

    # plot
    if average_channels
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

    lines!(
        ax,
        dat.data[!, :time],
        std(Matrix(dat.data[!, dat.layout.label]), dims = 2)[:],
        color = :black,
        linewidth = kwargs[:linewidth] * 2,
    )

    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    if isnothing(kwargs[:title])
        # ax.title = length(channels) == 1 ? "Electrode: $(channels[1])" : "Electrodes Avg: $(""*join(channels,",")*"")"
        ax.title = length(channels) == 1 ? "Electrode: $(channels[1])" : "Electrodes Avg: $(print_vector_(channels))"
    else
        ax.title = kwargs[:title]
    end
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]

    if kwargs[:add_topoplot]
        # just put in top left
        topo_ax =
            Axis(fig[kwargs[:topoplot_fig], 1], width = Relative(0.2), height = Relative(0.2), halign = 0, valign = 1)
        layout = filter(row -> row.label in String.(channels), dat.layout)
        if average_channels
            plot_layout_2d(fig, topo_ax, layout, point_kwargs = Dict(:color => kwargs[:color], :markersize => 18))
        else
            plot_layout_2d(
                fig,
                topo_ax,
                layout,
                point_kwargs = Dict(:colormap => colors, :color => 1:length(channels), :markersize => 18),
            )
        end
    end

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    display(fig)
    return fig, ax

end

function plot_erp(dat::ErpData, channels::Vector{Symbol}; average_channels = false, kwargs = Dict())
    fig = Figure()
    ax = Axis(fig[1, 1])
    fig, ax = plot_erp(fig, ax, dat, channels; average_channels = average_channels, kwargs = kwargs)
    return fig, ax
end

function plot_erp(dat::ErpData, channels::Symbol; average_channels = false, kwargs...)
    plot_erp(dat, [channels], average_channels = average_channels; kwargs...)
end

function plot_erp(dat::ErpData; average_channels = false, kwargs...)
    plot_erp(dat, dat.layout.label; average_channels = average_channels, kwargs...)
end



function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData, channels; average_channels = false, kwargs = Dict())
    println("plot_erp(dat_orig, dat_cleaned, channels)")
    kwargs = merge(kwargs, kwargs)
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    plot_erp(fig, ax1, dat_orig, channels; average_channels = average_channels, kwargs = kwargs)
    ax2 = Axis(fig[2, 1])
    kwargs[:topoplot_fig] = 2
    plot_erp(fig, ax2, dat_cleaned, channels; average_channels = average_channels, kwargs = kwargs)
    linkaxes!(ax1, ax2)
    display(fig)
    return fig, ax1, ax2
end


function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData; average_channels = false, kwargs...)
    plot_erp(dat_orig, dat_cleaned, dat_orig.layout.label; average_channels = average_channels, kwargs...)
end





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


function plot_grid_topo(
    dat::ErpData;
    plot_label_position = nothing,
    plot_label = true,
    plot_label_fontsize = 16,
    xlim = nothing,
    ylim = nothing,
    plot_height = 0.05,
    plot_width = 0.05,
    hide_decorations = false,
    show_scale = false,
    scale_position = [0.95, 0.05],
)

    # x/y limits
    isnothing(xlim) && (xlim = data_limits_x(dat.data, :time))
    isnothing(ylim) && (ylim = data_limits_y(dat.data, dat.layout.label))

    if plot_label && isnothing(plot_label_position)
        plot_label_positon = [xlim[1] ylim[2]]
    end

    xminmaxrange, yminmaxrange = datarange(dat.layout.x2), datarange(dat.layout.y2)
    xpositions = (layout.x2 ./ xminmaxrange) .+ 0.5
    ypositions = (layout.y2 ./ yminmaxrange) .+ 0.5
    fig = Figure()
    for (x, y, label) in zip(xpositions, ypositions, dat.layout.label)
        ax = Axis(fig[1, 1], width = Relative(plot_width), height = Relative(plot_height), halign = x, valign = y)
        lines!(ax, dat.data[!, :time], dat.data[!, label])
        vlines!(ax, [0], color = :black)
        hlines!(ax, [0], color = :black)
        if plot_label
            text!(
                ax,
                plot_label_positon[1],
                plot_label_positon[end],
                fontsize = plot_label_fontsize,
                text = label,
                align = (:left, :top),
            )
        end
        # hide some plot stuff (but this seems v. slow!)
        if hide_decorations
            hide_decorations!(ax)
        end
        xlims!(ax, xlim)
        ylims!(ax, ylim)
    end
    if show_scale
        ax = Axis(
            fig[1, 1],
            width = Relative(plot_width),
            height = Relative(plot_height),
            halign = scale_position[1],
            valign = scale_position[2],
        )
        lines!(ax, dat.data[!, :time], zeros(nrow(dat.data)))
        vlines!(ax, [0], color = :black)
        hlines!(ax, [0], color = :black)
        xlims!(ax, xlim)
        ylims!(ax, ylim)
        hidespines!(ax, :t, :r, :l, :b)
    end
    linkaxes!(filter(x -> x isa Axis, fig.content)...)
    display(fig)
    return fig, ax
end


# #################################################################
# plot_erp_image: 

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

function plot_correlation_heatmap(corr_df::DataFrame, mask_range::Union{Nothing,Tuple{Float64,Float64}} = nothing)
    """
    Plot a heatmap of the correlation matrix using Makie.

    Parameters:
    - corr_df: A DataFrame containing the correlation matrix with row and column names.
    """
    # Extract the correlation matrix (excluding the row names column)
    corr_matrix = Matrix(corr_df[:, 2:end])

    # Mask values within the specified range
    if !isnothing(mask_range)
        min_val, max_val = mask_range
        corr_matrix[(corr_matrix.>=min_val).&(corr_matrix.<=max_val)] .= NaN
    end

    # Extract row and column names
    row_names = String.(corr_df[!, :row])
    col_names = String.(propertynames(corr_df)[2:end])

    # Create the heatmap
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "Columns",
        ylabel = "Rows",
        xticks = (1:length(col_names), col_names),
        yticks = (1:length(row_names), row_names),
    )
    heatmap!(ax, corr_matrix, colormap = :viridis, colorrange = (-1, 1))

    # Add a colorbar
    Colorbar(fig[1, 2], limits = (-1, 1), label = "Correlation")

    # Display the figure
    display(fig)
    return fig, ax
end




function plot_ica_topoplot(
    ica,
    layout;
    comps = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)
    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end
    if isnothing(comps)
        comps = 1:size(ica.mixing)[2]
    end
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)
    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    label_default_kwargs =
        Dict(:plot_labels => false, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)
    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)
    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    fig = Figure()
    dims = best_rect(length(comps))
    count = 1
    axs = []
    for dim1 = 1:dims[1]
        for dim2 = 1:dims[2]
            ax = Axis(fig[dim1, dim2])
            push!(axs, ax)
            count += 1
            if count > length(comps)
                break
            end
        end
    end
    count = 1

    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]

    for ax in axs
        ax.title = String(ica.ica_label[comps[count]])
        data = data_interpolation_topo(
            ica.mixing[:, comps[count]],
            permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
            gridscale,
        )
        gridscale = gridscale
        radius = 88 # mm
        co = contourf!(
            ax,
            range(-radius * 2, radius * 2, length = gridscale),
            range(-radius * 2, radius * 2, length = gridscale),
            data,
            colormap = :jet,
        )
        # TODO: improve colorbar stuff
        # if plot_colorbar
        #     Colorbar(ax, co; colorbar_kwargs...)
        # end
        # head shape
        plot_layout_2d(
            fig,
            ax,
            layout,
            head_kwargs = head_kwargs,
            point_kwargs = point_kwargs,
            label_kwargs = label_kwargs,
        )
        count += 1
        if count > length(comps)
            break
        end
    end
    return fig
end

function plot_ica_topoplot(
    fig,
    ax,
    ica,
    comp,
    layout;
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)

    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs =
        Dict(:plot_labels => false, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]

    if ax.title.val == ""
        ax.title = ica.ica_label[comp]
    end
    data = data_interpolation_topo(ica.mixing[:, comp], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)
    gridscale = gridscale
    radius = 88 # mm
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data,
        colormap = :jet,
    )
    # TODO: improve colorbar stuff
    if plot_colorbar
        Colorbar(fig[1, 2], co; colorbar_kwargs...)
    end
    # head shape
    plot_layout_2d(fig, ax, layout, head_kwargs = head_kwargs, point_kwargs = point_kwargs, label_kwargs = label_kwargs)
    # end
    return fig
end

function plot_ica_component_activation(ica, epochs::EpochData, component::Int)
    fig = Figure()

    # Create grid layout and make it fill the figure
    gl = fig[1, 1] = GridLayout()
    colsize!(gl, 1, Relative(1.0))

    # Get data for plotting
    activation = ica.activation[component, :]
    projection = ica.mixing[:, component] .* activation'

    # Calculate y-limits for consistent scaling
    y_min, y_max = extrema(activation)
    y_range = y_max - y_min
    limits = (y_min - 0.1 * y_range, y_max + 0.1 * y_range)

    # Plot component activation
    ax1 = Axis(gl[1, 1], title = "Component $component Activation", limits = (nothing, limits))
    lines!(ax1, epochs.time, activation)

    # Plot channel projections
    ax2 = Axis(gl[2, 1], title = "Channel Projections")
    for (i, chan) in enumerate(epochs.layout.label)
        lines!(ax2, epochs.time, projection[i, :], color = :lightgrey, linewidth = 1, alpha = 0.5)
    end

    # Set row sizes to give equal space to plots
    rowsize!(gl, 1, Relative(0.5))
    rowsize!(gl, 2, Relative(0.5))

    rowgap!(gl, 10)  # Add gap between plots

    # Link x-axes
    linkxaxes!(ax1, ax2)

    return fig
end


function plot_ica_component_activation(
    dat::ContinuousData,
    ica_result::InfoIca,
    n_visible_components::Int = 10,  # Number of components visible at once
)
    # convert ContinuousData to appropriate matrix
    dat_matrix = permutedims(Matrix(dat.data[!, ica_result.data_label]))

    # Scale dat matrix the same way as in ICA
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica_result.scale

    # Transform data to component space
    components = ica_result.unmixing * dat_matrix
    total_components = size(components, 1)

    # Create figure
    fig = Figure()

    # Set up observables for interactive plotting
    window_size = 2000
    xrange = Observable(1:window_size)  # Initial window
    xlims = @lift((dat.data.time[first($xrange)], dat.data.time[last($xrange)]))

    # Initialize y-axis limits based on initial data window
    initial_range = maximum(abs.(extrema(components[1:n_visible_components, 1:window_size])))
    ylims = Observable((-initial_range, initial_range))

    # Observable for component range
    comp_start = Observable(1)

    # Create all subplots at once
    axs = []  # Store time series axes
    lines_obs = []  # Store line observables
    topo_axs = []  # Store topography axes

    # Create menu for selecting additional channels with label
    available_channels = names(dat.data)
    selected_channel = Observable("")  # Start with no channel selected

    # Create observable for selected channel data with proper initialization
    channel_data = Observable(zeros(size(dat.data, 1)))  # Initialize with zeros matching data length
    show_channel = Observable(false)  # Track whether to show channel data

    # Create observable for channel y-scale
    channel_yscale = Observable(1.0)  # Scaling factor for channel data

    for i = 1:n_visible_components
        # Topoplot
        ax_topo = Axis(
            fig[i, 1],
            width = 75,
            height = 75,
            title = @sprintf("IC %d (%.1f%%)", i, ica_result.variance[i] * 100)
        )
        push!(topo_axs, ax_topo)
        plot_ica_topoplot(fig, ax_topo, ica_result, i, layout, colorbar_kwargs = Dict(:plot_colorbar => false))
        hidexdecorations!(ax_topo)

        # Create subplot for time series with two y-axes
        ax_time = Axis(fig[i, 2], ylabel = "ICA Amplitude")
        ax_channel = Axis(
            fig[i, 2],
            ylabel = "Channel Amplitude",
            yaxisposition = :right,
            yticklabelcolor = :grey,
            ytickcolor = :grey,
        )

        # Link x-axes
        linkxaxes!(ax_time, ax_channel)

        # Hide the right axis' spine to avoid overlap
        hidespines!(ax_channel, :l, :t, :b)

        push!(axs, ax_time)


        # Create observable for component data
        line_obs = Observable(components[i, :])
        lines!(ax_time, @lift(dat.data.time[$xrange]), @lift($line_obs[$xrange]), color = :black)

        # Plot selected channel with scaling
        lines!(
            ax_channel,
            @lift(dat.data.time[$xrange]),
            @lift($show_channel ? $channel_data[$xrange] * $channel_yscale : fill(NaN, length($xrange))),
            color = :grey,
        )

        push!(lines_obs, line_obs)


        # Set initial limits
        xlims!(ax_time, xlims[])
        ylims!(ax_time, ylims[])
        ylims!(ax_channel, ylims[])  # Will be updated when channel is selected

        # Hide x-axis decorations for all but the last plot
        if i != n_visible_components
            hidexdecorations!(ax_time, grid = false)
        end

    end

    # Link all time series axes
    linkaxes!(axs...)

    # Adjust layout
    colsize!(fig.layout, 1, Auto(150))  # Fixed width for topo plots

    # Function to update component data
    function update_components(start_idx)
        for i = 1:n_visible_components
            comp_idx = start_idx + i - 1
            if comp_idx <= total_components
                lines_obs[i][] = components[comp_idx, :]

                # Clear and redraw topography
                empty!(topo_axs[i])
                plot_ica_topoplot(
                    fig,
                    topo_axs[i],
                    ica_result,
                    comp_idx,
                    layout,
                    colorbar_kwargs = Dict(:plot_colorbar => false),
                )

                # Update title
                topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, ica_result.variance[comp_idx] * 100)
            end
        end
    end

    # # Add navigation buttons below topo plots
    topo_nav = GridLayout(fig[end+1, 1])
    prev_topo = Button(topo_nav[1, 1], label = "◄ Previous")
    next_topo = Button(topo_nav[1, 2], label = "Next ►")

    # Create a new row in the figure layout for the menu
    menu_row = fig[end+1, 1]  # Add new row at the bottom
    menu_layout = GridLayout(menu_row)  # Create layout for the menu
    Label(menu_layout[1, 1], "Additional Channel", tellwidth = true)  # Centered label
    channel_menu = Menu(menu_layout[2, 1], options = ["None"; available_channels], default = "None")

    # Connect topo navigation buttons
    on(prev_topo.clicks) do _
        new_start = max(1, comp_start[] - n_visible_components)
        comp_start[] = new_start
        update_components(new_start)
    end

    on(next_topo.clicks) do _
        new_start = min(total_components - n_visible_components + 1, comp_start[] + n_visible_components)
        comp_start[] = new_start
        update_components(new_start)
    end

    # Add keyboard controls
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            if event.key == Keyboard.left || event.key == Keyboard.right
                # Handle x-axis scrolling
                current_range = xrange[]
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - window_size)
                    xrange[] = new_start:(new_start+window_size-1)
                else  # right
                    new_start = min(size(components, 2) - window_size + 1, first(current_range) + window_size)
                    xrange[] = new_start:(new_start+window_size-1)
                end

                # Update x-axis limits for all axes
                new_xlims = (dat.data.time[first(xrange[])], dat.data.time[last(xrange[])])
                for ax in axs
                    xlims!(ax, new_xlims)
                end

            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed =
                    (Keyboard.left_shift in events(fig).keyboardstate) ||
                    (Keyboard.right_shift in events(fig).keyboardstate)
                if !shift_pressed
                    # Handle y-axis scaling
                    current_range = ylims[][2]  # Just take the positive limit since it's symmetric
                    if event.key == Keyboard.up
                        # Zoom in - decrease range by 20%
                        new_range = current_range * 0.8
                    else  # down
                        # Zoom out - increase range by 20%
                        new_range = current_range * 1.2
                    end

                    # Keep centered on zero
                    new_ylims = (-new_range, new_range)
                    ylims[] = new_ylims

                    # Update y-axis limits for all axes
                    for ax in axs
                        ylims!(ax, new_ylims)
                    end
                else



                    if event.key == Keyboard.up && shift_pressed
                        channel_yscale[] = channel_yscale[] * 1.1
                    elseif event.key == Keyboard.down && shift_pressed
                        channel_yscale[] = channel_yscale[] / 1.1
                    end
                end


            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                # Handle component scrolling
                current_start = comp_start[]
                if event.key == Keyboard.page_up
                    new_start = max(1, current_start - n_visible_components)
                else  # page_down
                    new_start = min(total_components - n_visible_components + 1, current_start + n_visible_components)
                end

                if new_start != current_start
                    comp_start[] = new_start
                    update_components(new_start)
                end
            end
        end
    end

    # Update function for channel selection
    on(channel_menu.selection) do selected
        if selected == "None"
            show_channel[] = false
            channel_data[] = zeros(size(dat.data, 1))  # Reset to zeros when no channel selected
        else
            show_channel[] = true
            channel_data[] = dat.data[!, selected]  # Update with selected channel data
        end
    end


    return fig
end
