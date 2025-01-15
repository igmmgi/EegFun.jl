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

    return fig

end

function head_shape_2d(layout; kwargs...)
    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    head_shape_2d(fig, ax, layout; kwargs...)
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

    return fig

end

function head_shape_3d(layout; kwargs...)
    fig = Figure()
    ax = GLMakie.Axis3(fig[1, 1])
    head_shape_3d(fig, ax, layout; kwargs...)
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

    colorbar_default_kwargs = Dict(:plot_colorbar => true, width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)

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
    ax = GLMakie.Axis(fig[1, 1])

    radius = 88 # mm
    co = contourf!(
        range(-radius * 2, radius * 2, length = grid_scale),
        range(-radius * 2, radius * 2, length = grid_scale),
        data,
        levels = range(ylim[1], ylim[2], div(grid_scale, 2)),
        topo_kwargs...,
    )

    if plot_colour_bar
        Colorbar(fig[1, 2], co, colorbar_kwargs...)
    end

    # head shape
    head_shape_2d(fig, ax, dat.layout, head_kwargs..., point_kwargs..., label_kwargs...)

    return fig

end




##################################################################
# Data Browser: Continuous Data

struct ToggleButton
    label::String
    fun::Function
end

struct Marker
    data::Any
    line::Any
    text::Any
end

function add_marker!(markers, ax, data, col; label = nothing)
    marker_data = data[findall(x -> x != 0, data[!, col]), [:time, col]]
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




function plot_databrowser(dat::ContinuousData, channel_labels::Vector{<:AbstractString})

    function butterfly_plot(active)
        clear_axes()
        if active
            ycentre!()
            draw(plot_labels = false)
        elseif !active
            yoffset!()
            draw(plot_labels = true)
        end
    end

    function apply_lp_filter(active)
        clear_axes()
        data_filtered = nothing
        if active
            data_filtered = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
        end
        draw(plot_labels = true)
    end

    function plot_lines(marker, active)
        marker.line.visible = active
        marker.text.visible = active
        marker.text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker.data.time]
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

    function step_back()
        xrange.val[1] - 200 < 1 && return
        xrange[] = xrange.val .- 200
        xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
    end

    function step_forward()
        xrange.val[1] + 200 > nrow(data) && return
        xrange[] = xrange.val .+ 200
        xlims!(ax, data.time[xrange.val[1]], data.time[xrange.val[end]])
    end

    function chans_less()
        (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
        yrange.val = yrange[][1]+100:yrange[][end]-100
        ylims!(ax, yrange.val[1], yrange.val[end])
    end

    function chans_more()
        yrange.val = yrange[][1]-100:yrange[][end]+100
        ylims!(ax, yrange.val[1], yrange.val[end])
    end

    function yoffset!()
        for (idx, col) in enumerate(channel_labels)
            data[!, col] .+= offset.val[idx]
        end
    end

    function ycentre!()
        for (idx, col) in enumerate(channel_labels)
            data[!, col] .-= offset.val[idx]
        end
    end

    function clear_axes()
        [delete!(ax, value) for (key, value) in channel_data_original]
        [delete!(ax, value) for (key, value) in channel_data_filtered]
        [delete!(ax, value) for (key, value) in channel_data_labels]
    end

    # Makie Figure
    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])

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
    xlimit = 10000
    xrange = GLMakie.Observable(1:xlimit)
    yrange = GLMakie.Observable(-1500:1500)
    nchannels = length(channel_labels)
    channel_labels_original = channel_labels
    channels_to_plot = GLMakie.Observable(channel_labels)

    if nchannels > 1
        offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
    else # just centre
        offset = GLMakie.Observable(0.0)
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
                toggles[t, 3](markers[t-1], toggles[t, 1].active.val)
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
            options = vcat(["All", "Left", "Right", "Central"], dat.layout.label),
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
        clear_axes()

        data = copy(dat.data)
        if !isnothing(data_filtered)
            apply_lp_filter(true)
        end
        if nchannels > 1
            offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
        else # just centre
            offset = GLMakie.Observable(0.0)
        end
        yoffset!()
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
            event.key == Keyboard.left && step_back()
            event.key == Keyboard.right && step_forward()
            event.key == Keyboard.down && chans_less()
            event.key == Keyboard.up && chans_more()
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
    if ("is_vEOG" in names(dat.data) && "is_hEOG" in names(dat.data))
        add_marker!(markers, ax, data, :is_vEOG)
        add_marker!(markers, ax, data, :is_hEOG)
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
    yoffset!()
    draw(plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

function plot_databrowser(dat::ContinuousData)
    plot_databrowser(dat, dat.layout.label)
end

function plot_databrowser(dat::ContinuousData, channel_labels::Union{<:AbstractString,Vector})
    plot_databrowser(dat, [channel_labels])
end







###########################################################

function plot_databrowser(dat::EpochData, channel_labels::Vector{<:AbstractString})

    function butterfly_plot(active)
        clear_axes()
        if active
            ycentre!()
            draw(plot_labels = false)
        elseif !active
            yoffset!()
            draw(plot_labels = true)
        end
    end

    function apply_lp_filter(active)
        clear_axes()
        data_filtered = nothing
        if active
            data_filtered = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
        end
        draw(plot_labels = true)
    end

    function plot_trigger_lines(active)
        trigger_lines.visible = active
        trigger_text.visible = active
        trigger_text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.val.time]
    end

    function plot_vEOG_lines(active)
        vEOG_lines.visible = active
        vEOG_text.visible = active
        vEOG_text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.time]
    end

    function plot_hEOG_lines(active)
        hEOG_lines.visible = active
        hEOG_text.visible = active
        hEOG_text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.time]
    end

    function plot_extreme_lines(active)
        if length(extreme_spans) > 0
            extreme_spans[1].visible = active
        end
    end

    function toggle_button_group(fig, labels)
        toggles = []
        push!(toggles, ToggleButton("Butterfly Plot", butterfly_plot))
        "triggers" in labels && push!(toggles, ToggleButton("Trigger", plot_trigger_lines))
        "is_vEOG" in labels && push!(toggles, ToggleButton("vEOG", plot_vEOG_lines))
        "is_hEOG" in labels && push!(toggles, ToggleButton("hEOG", plot_hEOG_lines))
        "is_extreme" in labels && push!(toggles, ToggleButton("extreme", plot_extreme_lines))
        push!(toggles, ToggleButton("LP-Filter On/Off", apply_lp_filter))

        toggle_buttons = [Toggle(fig, active = false) for toggle in toggles]
        toggle_labels = [Label(fig, toggle.label, fontsize = 22, halign = :left) for toggle in toggles]
        toggle_functions = [toggle.fun for toggle in toggles]

        return hcat(toggle_buttons, toggle_labels, toggle_functions)

    end

    function step_epoch_forward()
        clear_axes()
        trial[] = min(length(dat.data), trial.val[1] + 1)
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        update_extreme_spans!()
        update_vEOG!()
        update_hEOG!()
        draw()
    end

    function step_epoch_backward()
        clear_axes()
        trial[] = max(1, trial.val[1] - 1)
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        update_extreme_spans!()
        update_vEOG!()
        update_hEOG!()
        draw()
    end

    function chans_less()
        (yrange.val[1] + 100 >= 0 || yrange.val[end] - 100 <= 0) && return
        yrange.val = yrange[][1]+100:yrange[][end]-100
        ylims!(ax, yrange.val[1], yrange.val[end])
    end

    function chans_more()
        yrange.val = yrange[][1]-100:yrange[][end]+100
        ylims!(ax, yrange.val[1], yrange.val[end])
    end

    function yoffset!()
        for t in eachindex(data)
            for (idx, col) in enumerate(channel_labels)
                data[t][!, col] .+= offset.val[idx]
            end
        end
    end

    function ycentre!()
        for t in eachindex(data)
            for (idx, col) in enumerate(channel_labels)
                data[t][!, col] .-= offset.val[idx]
            end
        end
    end

    function clear_axes()
        [delete!(ax, value) for (key, value) in channel_data_original]
        [delete!(ax, value) for (key, value) in channel_data_filtered]
        [delete!(ax, value) for (key, value) in channel_data_labels]
    end

    # Makie Figure
    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])

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
    xrange = GLMakie.Observable(1:xlimit)
    trial = GLMakie.Observable(1) # first trial
    yrange = GLMakie.Observable(-1500:1500)
    nchannels = length(channel_labels)
    channel_labels_original = channel_labels
    channels_to_plot = GLMakie.Observable(channel_labels)

    if nchannels > 1
        offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
    else # just centre
        offset = GLMakie.Observable(0.0)
    end

    xlims!(ax, data[1].time[xrange.val[1]], data[1].time[xrange.val[end]])
    ylims!(ax, yrange.val[1], yrange.val[end])
    ax.title = "Epoch $(trial.val)/$(length(dat.data))"
    ax.xlabel = "Time (S)"
    ax.ylabel = "Amplitude (mV)"

    # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
    toggles = toggle_button_group(fig, names(data[1]))
    for t = 1:length(toggles[:, 1])
        on(toggles[t, 1].active) do _
            toggles[t, 3](toggles[t, 1].active.val)
        end
    end

    # menu for electrode/channel selection
    menu = hcat(
        Menu(
            fig,
            options = vcat(["All", "Left", "Right", "Central"], dat.layout.label),
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
        clear_axes()

        data = deepcopy(dat.data)
        if !isnothing(data_filtered)
            apply_lp_filter(true)
        end
        if nchannels > 1
            offset = GLMakie.Observable(LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels))
        else # just centre
            offset = GLMakie.Observable(0.0)
        end
        yoffset!()
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
        clear_axes()
        trial[] = s
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        update_extreme_spans!()
        update_vEOG!()
        update_hEOG!()
        draw()
    end

    slider_extreme = Slider(fig[1, 2], range = 0:5:100, startvalue = 200, width = 100)
    slider_lp_filter = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)

    crit_val = lift(slider_extreme.value) do x
        x
    end


    # keyboard events
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            #if event.action in (Keyboard.press,)
            event.key == Keyboard.left && step_epoch_backward()
            event.key == Keyboard.right && step_epoch_forward()
            event.key == Keyboard.down && chans_less()
            event.key == Keyboard.up && chans_more()
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


    #################### Triggers/Events ###############################
    trigger_data_time = @lift data[$trial][findall(x -> x != 0, data[$trial][!, :].triggers), [:time, :triggers]]
    trigger_lines = vlines!(trigger_data_time.val.time, color = :grey, linewidth = 1, visible = false)
    trigger_text = text!(
        string.(trigger_data_time.val.triggers),
        position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in trigger_data_time.val.time],
        space = :data,
        align = (:center, :center),
        fontsize = 22,
        visible = false,
    )

    ################### vEOG/hEOG ###############################
    vEOG_data_time = []
    vEOG_lines = []
    vEOG_text = []
    function update_vEOG!()
        if length(vEOG_lines) > 0
            delete!(ax, vEOG_lines)
            delete!(ax, vEOG_text)
        end
        vEOG_data_time = data[trial.val][findall(x -> x != 0, data[trial.val][!, :].is_vEOG), [:time, :is_vEOG]]
        vEOG_lines = vlines!(vEOG_data_time.time, color = :grey, linewidth = 1, visible = false)
        vEOG_text = text!(
            repeat(["v"], nrow(vEOG_data_time)),
            position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in vEOG_data_time.time],
            space = :data,
            align = (:center, :center),
            fontsize = 22,
            visible = false,
        )
        plot_vEOG_lines(toggles[3, 1].active.val)
    end
    if ("is_vEOG" in names(dat.data[1]))
        update_vEOG!()
    end

    hEOG_data_time = []
    hEOG_lines = []
    hEOG_text = []
    function update_hEOG!()
        if length(hEOG_lines) > 0
            delete!(ax, hEOG_lines)
            delete!(ax, hEOG_text)
        end
        hEOG_data_time = data[trial.val][findall(x -> x != 0, data[trial.val][!, :].is_hEOG), [:time, :is_hEOG]]
        hEOG_lines = vlines!(hEOG_data_time.time, color = :grey, linewidth = 1, visible = false)
        hEOG_text = text!(
            repeat(["h"], nrow(hEOG_data_time)),
            position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in hEOG_data_time.time],
            space = :data,
            align = (:center, :center),
            fontsize = 22,
            visible = false,
        )
        plot_hEOG_lines(toggles[4, 1].active.val)
    end
    if ("is_hEOG" in names(dat.data[1]))
        update_hEOG!()
    end

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
    yoffset!()
    draw(plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

function plot_databrowser(dat::EpochData)
    plot_databrowser(dat, dat.layout.label)
end

function plot_databrowser(dat::EpochData, channel_labels::Union{<:AbstractString,Vector})
    plot_databrowser(dat, [channel_labels])
end



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
    ax = GLMakie.Axis(fig[1, 1])

    # plot each trial (average acrossed electordes if >1) and overall average
    avg_data = zeros(nrow(dat.data[1]))
    for trial in eachindex(dat.data)
        trial_data = colmeans(dat.data[trial], channels)
        avg_data .+= trial_data
        GLMakie.lines!(
            dat.data[trial][!, :time],
            trial_data,
            color = kwargs[:color][1],
            linewidth = kwargs[:linewidth][1],
        )
    end
    avg_data ./= length(dat.data)
    GLMakie.lines!(dat.data[1][!, :time], avg_data, color = kwargs[:color][2], linewidth = kwargs[:linewidth][2])

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

    return fig

end

function plot_epochs(dat::EpochData, channels::Union{AbstractString,Symbol}; kwargs...)
    plot_epochs(dat::EpochData, [channels]; kwargs...)
end


# #################################################################
# plot_erp: ERP Data (Single Condition; Single Channel or Average of multiple channels)
function plot_erp(dat::ErpData, channels::Union{Vector{<:AbstractString},Vector{Symbol}}; kwargs = Dict())

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
    ax = GLMakie.Axis(fig[1, 1])

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

    return fig

end

function plot_erp(dat::ErpData, channels::Union{AbstractString,Symbol}; kwargs...)
    plot_erp(dat, [channels]; kwargs...)
end


# #################################################################
# plot_grid_rect: 
function plot_grid_rect(dat::ErpData; channels = nothing, kwargs = Dict())
    if isnothing(channels)
        channels = dat.layout.label
    end
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
    if (isnothing(kwargs[:xlim]))
        kwargs[:xlim] = data_limits_x(dat.data)
    end
    if (isnothing(kwargs[:ylim]))
        kwargs[:ylim] = data_limits_y(dat.data, dat.layout.label)
    end
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
    return fig
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
    if (isnothing(xlim))
        xlim = data_limits_x(dat.data, :time)
    end
    if (isnothing(ylim))
        ylim = data_limits_y(dat.data, dat.layout.label)
    end

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
    return fig
end


# #################################################################
# plot_erp_image: 

function plot_erp_image(dat::EpochData, channels::Vector{Symbol}; colorrange = nothing)
    data = zeros(length(dat.data), nrow(dat.data[1]))
    for epoch in eachindex(dat.data)
        data[epoch, :] = colmeans(dat.data[epoch], channels)
    end
    if isnothing(colorrange)
        colorrange = extrema(data)
    end
    fig = Figure()
    ax = GLMakie.Axis(fig[1, 1])
    heatmap!(ax, dat.data[1].time, 1:length(dat.data), data, colorrange = colorrange)
    ax = GLMakie.Axis(fig[2, 1])
    lines!(ax, dat.data[1].time, colmeans(data))
    xlims!(ax, (-0.5, 2))
    return fig
end

function plot_erp_image(dat::EpochData, channel)
    plot_erp_image(dat, [channel])
end
