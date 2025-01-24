# TODO: butterfly plot/global field power

#########################################
# 2D head shape
function head_shape_2d(fig, ax, layout; head_kwargs = Dict(), point_kwargs = Dict(), label_kwargs = Dict())

    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = pop!(point_kwargs, :plot_points)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = pop!(label_kwargs, :plot_labels)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    # head shape
    radius = 88 # mm
    arc!(ax, Point2f(0), radius * 2, -π, π; head_kwargs...) # head
    arc!(Point2f(radius * 2, 0), radius * 2 / 7, -π / 2, π / 2; head_kwargs...) # ear right
    arc!(Point2f(-radius * 2, 0), -radius * 2 / 7, π / 2, -π / 2; head_kwargs...) # ear left
    lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * 4; head_kwargs...) # nose

    # points
    if plot_points
        scatter!(ax, layout[!, :x2], layout[!, :y2]; point_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(ax, position = (label.x2 + xoffset, label.y2 + yoffset), label.label; label_kwargs...)
        end
    end

    # hide some plot stuff
    hidedecorations!(ax)
    hidespines!(ax)

    display(fig)
    return fig, ax

end

function head_shape_2d(layout; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    head_shape_2d(fig, ax, layout; kwargs...)
    display(fig)
    return fig, ax
end

# using LibGEOS package (currently convexhull)
function point_border(xpos, ypos, border_size)
    circle_points = 0:2*pi/361:2*pi
    xs = (border_size.*sin.(circle_points).+transpose(xpos))[:]
    ys = (border_size.*cos.(circle_points).+transpose(ypos))[:]
    xys = [[xs[i], ys[i]] for i in eachindex(xs)]
    push!(xys, xys[1])
    poly = LibGEOS.Polygon([xys])
    hull = LibGEOS.convexhull(poly)
    return hull
end

function add_topo_rois!(fig, layout, rois; border_size = 10, roi_kwargs = Dict())
    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end
    roi_default_kwargs = Dict(:color => repeat([:black], length(rois)), :linewidth => repeat([2], length(rois)))
    roi_kwargs = merge(roi_default_kwargs, roi_kwargs)
    for (idx, roi) in enumerate(rois)
        xpos = filter(row -> row.label ∈ roi, layout).x2
        ypos = filter(row -> row.label ∈ roi, layout).y2
        border = point_border(xpos, ypos, border_size)
        lines!(border, linewidth = roi_kwargs[:linewidth][idx], color = roi_kwargs[:color][idx])
    end
end

#########################################
# 3D head shape
function head_shape_3d(fig, ax, layout; point_kwargs = Dict(), label_kwargs = Dict())

    if (:x3 ∉ names(layout) || :y3 ∉ names(layout) || :z3 ∉ names(layout))
        polar_to_cartesian_xyz!(layout)
    end

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = pop!(point_kwargs, :plot_points)

    label_default_kwargs = Dict(
        :plot_labels => true,
        :fontsize => 20,
        :color => :black,
        :color => :black,
        :xoffset => 0,
        :yoffset => 0,
        :zoffset => 0,
    )
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = pop!(label_kwargs, :plot_labels)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)
    zoffset = pop!(label_kwargs, :zoffset)

    # points
    if plot_points
        scatter!(ax, layout[!, :x3], layout[!, :y3], layout[!, :z3]; point_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(
                ax,
                position = (label.x3 + xoffset, label.y3 + yoffset, label.z3 + zoffset),
                label.label;
                label_kwargs...,
            )
        end
    end

    # hide some plot stuff
    hidedecorations!(ax)
    hidespines!(ax)

    display(fig)
    return fig, ax

end

function head_shape_3d(layout; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    head_shape_3d(fig, ax, layout; kwargs...)
    display(fig)
    return fig, ax
end


##########################################
# 2D topographic plot
function plot_topoplot(
    dat;
    xlim = nothing,
    ylim = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)

    if (:x2 ∉ names(dat.layout) || :y2 ∉ names(dat.layout))
        polar_to_cartesian_xy!(dat.layout)
    end

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs =
        Dict(:plot_labels => true, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    if isnothing(xlim)
        xlim = [dat.data.time[1], dat.data.time[end]]
        xlim_idx = 1:nrow(dat.data)
    end

    # convert xlim to index
    xlim_idx = find_idx_range(dat.data.time, xlim[1], xlim[2])

    # interpolate data
    data = data_interpolation_topo(
        mean.(eachcol(dat.data[xlim_idx, dat.layout.label])),
        Matrix(dat.layout[!, [:x2, :y2]])',
        gridscale,
    )

    if isnothing(ylim)
        ylim = minimum(data[.!isnan.(data)]), maximum(data[.!isnan.(data)])
    end

    fig = Figure()
    ax = Axis(fig[1, 1])

    radius = 88 # mm
    co = contourf!(
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data,
        levels = range(ylim[1], ylim[2], div(gridscale, 2));
        topo_kwargs...,
    )

    if plot_colorbar
        Colorbar(fig[1, 2], co; colorbar_kwargs...)
    end

    # head shape
    head_shape_2d(
        fig,
        ax,
        dat.layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    return fig

end




##################################################################
# Data Browser: Continuous Data
struct ToggleButton
    label::Any
    fun::Any
end

struct Marker
    data::Any
    line::Any
    text::Any
end

function yless!(ax, yrange)
    (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
    yrange.val = yrange[][1]+100:yrange[][end]-100
    ylims!(ax, yrange.val[1], yrange.val[end])
end

function ymore!(ax, yrange)
    yrange.val = yrange[][1]-100:yrange[][end]+100
    ylims!(ax, yrange.val[1], yrange.val[end])
end

function xback!(ax, xrange, data)
    xrange.val[1] - 200 < 1 && return
    xrange[] = xrange.val .- 200
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
end

function xforward!(ax, xrange, data)
    xrange.val[1] + 200 > nrow(data) && return
    xrange[] = xrange.val .+ 200
    xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
end

function yoffset!(data::DataFrame, channel_labels, offset)
    for (idx, col) in enumerate(channel_labels)
        data[!, col] .+= offset.val[idx]
    end
end

function ycentre!(data::DataFrame, channel_labels, offset)
    for (idx, col) in enumerate(channel_labels)
        data[!, col] .-= offset.val[idx]
    end
end

function yoffset!(data, channel_labels, offset)
    for t in eachindex(data)
        for (idx, col) in enumerate(channel_labels)
            data[t][!, col] .+= offset.val[idx]
        end
    end
end

function ycentre!(data, channel_labels, offset)
    for t in eachindex(data)
        for (idx, col) in enumerate(channel_labels)
            data[t][!, col] .-= offset.val[idx]
        end
    end
end

clear_axes(ax, datas) = [delete!(ax, value) for data in datas for (key, value) in data]


function add_marker!(markers, ax, data, col; label = nothing, trial = nothing)
    if isnothing(trial)
        marker_data = data[findall(x -> x != 0, data[!, col]), [:time, col]]
    else
        marker_data = data[trial][findall(x -> x != 0, data[trial][!, col]), [:time, col]]
    end
    if isnothing(label)
        label = string.(marker_data[!, col])
    else
        label = repeat([label], nrow(marker_data))
    end
    push!(
        markers,
        Marker(
            marker_data,
            vlines!(marker_data.time, color = :grey, linewidth = 1, visible = false),
            text!(
                label,
                position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker_data.time],
                space = :data,
                align = (:center, :center),
                fontsize = 22,
                visible = false,
            ),
        ),
    )
end




function plot_lines(ax, marker, active)
    marker.line.visible = active
    marker.text.visible = active
    marker.text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker.data.time] # incase y changed
end


function plot_databrowser(dat::ContinuousData, channel_labels::Vector{<:AbstractString})

    function butterfly_plot(active)
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        if active
            ycentre!(data, channel_labels, offset)
            draw(plot_labels = false)
        elseif !active
            yoffset!(data, channel_labels, offset)
            draw(plot_labels = true)
        end
    end

    function apply_lp_filter(active)
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        data_filtered = nothing
        if active
            data_filtered = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
        end
        draw(plot_labels = true)
    end


    function plot_extreme_lines(active)
        extreme_spans.visible = active
    end

    function toggle_button_group(fig, labels)
        toggles = []
        push!(toggles, ToggleButton("Butterfly Plot", butterfly_plot))
        "triggers" in labels && push!(toggles, ToggleButton("Trigger", plot_lines))
        "is_vEOG" in labels && push!(toggles, ToggleButton("vEOG", plot_lines))
        "is_hEOG" in labels && push!(toggles, ToggleButton("hEOG", plot_lines))
        "is_extreme" in labels && push!(toggles, ToggleButton("extreme", plot_extreme_lines))
        push!(toggles, ToggleButton("LP-Filter On/Off", apply_lp_filter))

        toggle_buttons = [Toggle(fig, active = false) for toggle in toggles]
        toggle_labels = [Label(fig, toggle.label, fontsize = 22, halign = :left) for toggle in toggles]
        toggle_functions = [toggle.fun for toggle in toggles]

        return hcat(toggle_buttons, toggle_labels, toggle_functions)

    end

    # Makie Figure
    fig = Figure()
    ax = Axis(fig[1, 1])

    # controls
    # interactions(ax)
    # deregister_interaction!(ax, :rectanglezoom)
    # deregister_interaction!(ax, :dragpan)
    # deregister_interaction!(ax, :scrollzoom)
    # deregister_interaction!(ax, :limitreset)

    # data to plot
    data = copy(dat.data)
    data_filtered = nothing

    channel_data_original = Dict()
    channel_data_filtered = Dict()
    channel_data_labels = Dict()

    # default xrange/yrange
    xlimit = 5000
    xrange = Observable(1:xlimit)
    yrange = Observable(-1500:1500)
    nchannels = length(channel_labels)
    channel_labels_original = channel_labels

    if nchannels > 1
        offset = Observable(LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, nchannels + 2)[2:end-1])
    else # just centre
        offset = Observable(0.0)
    end

    @lift xlims!(ax, data.time[$xrange[1]], data.time[$xrange[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
    ax.xlabel = "Time (S)"
    ax.ylabel = "Amplitude (mV)"

    # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
    toggles = toggle_button_group(fig, names(data))
    for t = 1:length(toggles[:, 1])
        if toggles[t, 2].text.val ∈ ["Trigger", "vEOG", "hEOG"]
            on(toggles[t, 1].active) do _
                toggles[t, 3](ax, markers[t-1], toggles[t, 1].active.val)
            end
        else
            on(toggles[t, 1].active) do _
                toggles[t, 3](toggles[t, 1].active.val)
            end
        end
    end

    # menu for electrode/channel selection
    menu = hcat(
        Menu(
            fig,
            options = vcat(["All", "Left", "Right", "Central"], channel_labels_original),
            default = "All",
            direction = :down,
            fontsize = 18,
        ),
        Label(fig, "Labels", fontsize = 22, halign = :left),
    )
    on(menu[1].selection) do s
        channel_labels = [s]
        if s == "All"
            channel_labels = channel_labels_original
        elseif s == "Left"
            channel_labels = channel_labels_original[findall(occursin.(r"\d*[13579]$", channel_labels_original))]
        elseif s == "Right"
            channel_labels = channel_labels_original[findall(occursin.(r"\d*[24680]$", channel_labels_original))]
        elseif s == "Central"
            channel_labels = channel_labels_original[findall(occursin.(r"z$", channel_labels_original))]
        end
        nchannels = length(channel_labels)

        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])

        data = copy(dat.data)
        if !isnothing(data_filtered)
            apply_lp_filter(true)
        end
        if nchannels > 1
            offset = Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels + 2)[2:end-1])
        else # just centre
            offset = Observable(0.0)
        end
        yoffset!(data, channel_labels, offset)
        if !isnothing(data_filtered)
            apply_lp_filter(true)
        else
            draw(plot_labels = true)
        end
    end

    slider_extreme = Slider(fig[1, 2], range = 0:5:100, startvalue = 100, width = 100)
    slider_lp_filter = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)

    crit_val = lift(slider_extreme.value) do x
        x
    end

    slider_range = Slider(fig[3, 1], range = 100:30000, startvalue = xlimit, snap = false)
    slider_x =
        Slider(fig[2, 1], range = slider_range.value.val:nrow(data), startvalue = slider_range.value, snap = false)

    on(slider_x.value) do x
        xrange[] = max(1, x - slider_range.value.val):x-1
        xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
    end

    on(slider_range.value) do x
        xrange.val = max(1, slider_x.value.val - x):slider_x.value.val-1
        xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
    end

    # keyboard events
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat)
            event.key == Keyboard.left && xback!(ax, xrange, data)
            event.key == Keyboard.right && xforward!(ax, xrange, data)
            event.key == Keyboard.down && yless!(ax, yrange)
            event.key == Keyboard.up && ymore!(ax, yrange)
        end
        # TODO: what is best here?
        # return Consume()
        # return Consume(false)
    end

    # position GUI controls
    fig[1, 2] = grid!(
        vcat(
            toggles[:, 1:2],
            hcat(
                slider_lp_filter,
                Label(fig, @lift("LP-Filter: $($(slider_lp_filter.value)) Hz"), fontsize = 22, halign = :left),
            ),
            hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)),
            menu,
        ),
        tellheight = false,
    )
    colsize!(fig.layout, 2, Relative(1 / 6))

    ################### vertical line markers ###############################
    # Vertical line markers
    markers = []
    add_marker!(markers, ax, data, :triggers)
    # if ("is_vEOG" in names(dat.data[1]) && "is_hEOG" in names(dat.data[1]))
    if ("is_vEOG" in names(dat.data) && "is_hEOG" in names(dat.data))
        add_marker!(markers, ax, data, :is_vEOG, label = "v")
        add_marker!(markers, ax, data, :is_hEOG, label = "h")
    end

    ################### Extreme Values ###############################
    if ("is_extreme" in names(data))
        extreme = @views splitgroups(findall(x -> x != 0, data[!, :].is_extreme))
        extreme_spans = vspan!(
            ax,
            data[extreme[1], :time],
            data[extreme[2], :time],
            color = "LightGrey",
            alpha = 0.5,
            visible = false,
        )
    end

    function draw(; plot_labels = true)
        alpha_orig = 1
        linewidth_orig = 2
        plot_filtered_data = !isnothing(data_filtered)
        if plot_filtered_data
            alpha_orig = 0.5
            linewidth_orig = 1
        end
        for col in channel_labels
            # original data
            channel_data_original[col] = lines!(
                ax,
                data[!, :time],
                data[!, col],
                color = @lift(abs.(dat.data[!, col]) .>= $crit_val),
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = linewidth_orig,
                alpha = alpha_orig,
            )
            if plot_labels
                channel_data_labels[col] = text!(
                    ax,
                    @lift(data[$xrange, :time][1]),
                    @lift(data[$xrange, col][1]),
                    text = col,
                    align = (:left, :center),
                    fontsize = 18,
                )
            end
            # also show filtered data
            if plot_filtered_data
                channel_data_filtered[col] =
                    lines!(ax, data_filtered[!, :time], data_filtered[!, col], color = :blue, linewidth = 2)
            end
        end
    end

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    hideydecorations!(ax, label = true)
    yoffset!(data, channel_labels, offset)
    draw(plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

plot_databrowser(dat::ContinuousData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::ContinuousData, channel_labels::Union{<:AbstractString,Vector}) =
    plot_databrowser(dat, [channel_labels])



###########################################################

function plot_databrowser(dat::EpochData, channel_labels::Vector{<:AbstractString})

    function butterfly_plot(active)
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        if active
            ycentre!(data, channel_labels, offset)
            draw(plot_labels = false)
        elseif !active
            yoffset!(data, channel_labels, offset)
            draw(plot_labels = true)
        end
    end

    function apply_lp_filter(active)
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        data_filtered = nothing
        if active
            data_filtered = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
        end
        draw(plot_labels = true)
    end

    function update_markers!(markers)
        for marker in markers
            delete!(ax, marker.line)
            delete!(ax, marker.text)
        end
        empty!(markers)
        add_marker!(markers, ax, data, :triggers, trial = trial.val)
        if ("is_vEOG" in names(dat.data[trial.val]) && "is_hEOG" in names(dat.data[trial.val]))
            add_marker!(markers, ax, data, :is_vEOG, trial = trial.val, label = "v")
            add_marker!(markers, ax, data, :is_hEOG, trial = trial.val, label = "h")
        end
    end

    function plot_extreme_lines(active)
        if length(extreme_spans) > 0
            extreme_spans[1].visible = active
        end
    end

    function toggle_button_group(fig, labels)
        toggles = []
        push!(toggles, ToggleButton("Butterfly Plot", butterfly_plot))
        "triggers" in labels && push!(toggles, ToggleButton("Trigger", plot_lines))
        "is_vEOG" in labels && push!(toggles, ToggleButton("vEOG", plot_lines))
        "is_hEOG" in labels && push!(toggles, ToggleButton("hEOG", plot_lines))
        "is_extreme" in labels && push!(toggles, ToggleButton("extreme", plot_extreme_lines))
        push!(toggles, ToggleButton("LP-Filter On/Off", apply_lp_filter))

        toggle_buttons = [Toggle(fig, active = false) for toggle in toggles]
        toggle_labels = [Label(fig, toggle.label, fontsize = 22, halign = :left) for toggle in toggles]
        toggle_functions = [toggle.fun for toggle in toggles]

        return hcat(toggle_buttons, toggle_labels, toggle_functions)

    end

    function step_epoch_forward()
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        trial[] = min(length(dat.data), trial.val[1] + 1)
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        update_extreme_spans!()
        update_markers!(markers)
        draw()
    end

    function step_epoch_backward()
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        trial[] = max(1, trial.val[1] - 1)
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        update_extreme_spans!()
        update_markers!(markers)
        draw()
    end

    # Makie Figure
    fig = Figure()
    ax = Axis(fig[1, 1])

    # controls
    # interactions(ax)
    # deregister_interaction!(ax, :rectanglezoom)
    # deregister_interaction!(ax, :dragpan)
    # deregister_interaction!(ax, :scrollzoom)
    # deregister_interaction!(ax, :limitreset)

    # data to plot
    data = deepcopy(dat.data)
    data_filtered = nothing

    channel_data_original = Dict()
    channel_data_filtered = Dict()
    channel_data_labels = Dict()

    # default xrange/yrange
    xlimit = nrow(dat.data[1])
    xrange = Observable(1:xlimit)
    trial = Observable(1) # first trial
    yrange = Observable(-1500:1500)
    nchannels = length(channel_labels)
    channel_labels_original = channel_labels

    if nchannels > 1
        offset = Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels + 2)[2:end-1])
    else # just centre
        offset = Observable(0.0)
    end

    xlims!(ax, data[1].time[xrange.val[1]], data[1].time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
    ax.xlabel = "Time (S)"
    ax.ylabel = "Amplitude (mV)"

    # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
    toggles = toggle_button_group(fig, names(data[trial.val]))
    for t = 1:length(toggles[:, 1])
        if toggles[t, 2].text.val ∈ ["Trigger", "vEOG", "hEOG"]
            on(toggles[t, 1].active) do _
                toggles[t, 3](ax, markers[t-1], toggles[t, 1].active.val)
            end
        else
            on(toggles[t, 1].active) do _
                toggles[t, 3](toggles[t, 1].active.val)
            end
        end
    end

    # menu for electrode/channel selection
    menu = hcat(
        Menu(
            fig,
            options = vcat(["All", "Left", "Right", "Central"], channel_labels_original),
            default = "All",
            direction = :down,
            fontsize = 18,
        ),
        Label(fig, "Labels", fontsize = 22, halign = :left),
    )
    on(menu[1].selection) do s
        channel_labels = [s]
        if s == "All"
            channel_labels = channel_labels_original
        elseif s == "Left"
            channel_labels = channel_labels_original[findall(occursin.(r"\d*[13579]$", channel_labels_original))]
        elseif s == "Right"
            channel_labels = channel_labels_original[findall(occursin.(r"\d*[24680]$", channel_labels_original))]
        elseif s == "Central"
            channel_labels = channel_labels_original[findall(occursin.(r"z$", channel_labels_original))]
        end
        nchannels = length(channel_labels)
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])

        data = deepcopy(dat.data)
        if !isnothing(data_filtered)
            apply_lp_filter(true)
        end
        if nchannels > 1
            offset = Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels + 2)[2:end-1])
        else # just centre
            offset = Observable(0.0)
        end
        yoffset!(data, channel_labels, offset)
        if !isnothing(data_filtered)
            apply_lp_filter(true)
        else
            draw(plot_labels = true)
        end
    end

    menu_trial = hcat(
        Menu(fig, options = 1:length(data), default = 1, direction = :down, fontsize = 18),
        Label(fig, "Epoch", fontsize = 22, halign = :left),
    )
    on(menu_trial[1].selection) do s
        clear_axes(ax, [channel_data_original, channel_data_filtered, channel_data_labels])
        trial[] = s
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        update_extreme_spans!()
        update_markers!(markers)
        draw()
    end

    slider_extreme   = Slider(fig[1, 2], range = 0:5:100, startvalue = 200, width = 100)
    slider_lp_filter = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)

    crit_val = lift(slider_extreme.value) do x
        x
    end

    # keyboard events
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            event.key == Keyboard.left && step_epoch_backward()
            event.key == Keyboard.right && step_epoch_forward()
            event.key == Keyboard.down && yless!(ax, yrange)
            event.key == Keyboard.up && ymore!(ax, yrange)
        end
        # TODO: what is best here?
        # return Consume()
        # return Consume(false)
    end

    # position GUI controls
    fig[1, 2] = grid!(
        vcat(
            toggles[:, 1:2],
            hcat(
                slider_lp_filter,
                Label(fig, @lift("LP-Filter: $($(slider_lp_filter.value)) Hz"), fontsize = 22, halign = :left),
            ),
            hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)),
            menu,
            menu_trial,
        ),
        tellheight = false,
    )
    colsize!(fig.layout, 2, Relative(1 / 6))


    ################### vertical line markers ###############################
    # Vertical line markers
    markers = []
    update_markers!(markers)

    extreme_spans = []
    function update_extreme_spans!()
        if length(extreme_spans) > 0
            delete!(ax, extreme_spans[1])
            extreme_spans = []
        end
        tmp = findall(x -> x != 0, data[trial.val][!, :].is_extreme)
        if length(tmp) > 0
            extreme = splitgroups(tmp)
            if length(extreme) > 0
                # TODO: hard coded Toggle index!!!
                push!(
                    extreme_spans,
                    vspan!(
                        ax,
                        data[trial.val][extreme[1], :time],
                        data[trial.val][extreme[2], :time],
                        color = "LightGrey",
                        alpha = 0.5,
                        visible = toggles[5, 1].active.val,
                    ),
                )
            end
        end
    end

    #################### Extreme Values ###############################
    if ("is_extreme" in names(data[1]))
        update_extreme_spans!()
    end

    function draw(; plot_labels = true)
        alpha_orig = 1
        linewidth_orig = 2
        plot_filtered_data = !isnothing(data_filtered)
        if plot_filtered_data
            alpha_orig = 0.5
            linewidth_orig = 1
        end
        for col in channel_labels
            # original data
            channel_data_original[col] = lines!(
                ax,
                data[trial.val][!, :time],
                data[trial.val][!, col],
                color = @lift(abs.(dat.data[trial.val][!, col]) .>= $crit_val),
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = linewidth_orig,
                alpha = alpha_orig,
            )
            if plot_labels
                channel_data_labels[col] = text!(
                    ax,
                    @lift(data[trial.val][$xrange, :time][1]),
                    @lift(data[trial.val][$xrange, col][1]),
                    text = col,
                    align = (:left, :center),
                    fontsize = 20,
                )
            end
            # also show filtered data
            if plot_filtered_data
                channel_data_filtered[col] = lines!(
                    ax,
                    data_filtered[trial.val][!, :time],
                    data_filtered[trial.val][!, col],
                    color = :blue,
                    linewidth = 2,
                )
            end
        end
    end

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    hideydecorations!(ax, label = true)
    yoffset!(data, channel_labels, offset)
    draw(plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

plot_databrowser(dat::EpochData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::EpochData, channel_labels::Union{<:AbstractString,Vector}) =
    plot_databrowser(dat, [channel_labels])



# #################################################################
# plot_epochs: Epoched Data (Single Condition; Single Channel or Average of multiple channels)
function plot_epochs(dat::EpochData, channels::Union{Vector{<:AbstractString},Vector{Symbol}}; kwargs = Dict())

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

plot_epochs(dat::EpochData, channels::Union{AbstractString,Symbol}; kwargs...) =
    plot_epochs(dat::EpochData, [channels]; kwargs...)


# #################################################################
# plot_erp: ERP Data (Single Condition; Single Channel or Average of multiple channels)
function plot_erp(dat::ErpData, channels::Union{Vector{<:AbstractString},Vector{Symbol}}; kwargs = Dict())

    default_kwargs = Dict(
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => 2,
        :color => :black,
        :yreversed => false,
    )

    kwargs = merge(default_kwargs, kwargs)

    fig = Figure()
    ax = Axis(fig[1, 1])

    # plot
    lines!(ax, dat.data[!, :time], colmeans(dat.data, channels), color = kwargs[:color], linewidth = kwargs[:linewidth])

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

plot_erp(dat::ErpData, channels::Union{AbstractString,Symbol}; kwargs...) = plot_erp(dat, [channels]; kwargs...)


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
    channels::Union{Vector{<:AbstractString},Vector{Symbol}};
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


function plot_erp_image(dat::EpochData, channel::Union{AbstractString,Symbol})
    plot_erp_image(dat, [channel])
end
