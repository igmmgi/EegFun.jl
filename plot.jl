# TODO: butterfly plot/global field power
# TODO: spline interpolation for topoplots?

function _trigger_time_count(time, triggers)
    trigger_values = triggers[findall(diff(triggers) .>= 1).+1]
    trigger_times = time[findall(diff(triggers) .>= 1).+1]
    trigger_count = OrderedDict(i => 0 for i in sort!(collect(Set(trigger_values))))
    for val in trigger_values
        trigger_count[val] += 1
    end
    return trigger_times, trigger_values, trigger_count
end

function plot_events(dat::BioSemiBDF.BioSemiData)
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.time, dat.triggers.raw)
    plot_events(trigger_times, trigger_values, trigger_count)
end

function plot_events(dat::ContinuousData)
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.data.time, dat.data.triggers)
    plot_events(trigger_times, trigger_values, trigger_count)
end

function plot_events(trigger_times, trigger_values, trigger_count)
    fig = Figure()
    ax = Axis(fig[1, 1], yticks = (1:length(trigger_count.keys), string.(trigger_count.keys)))
    for (unique, (key, value)) in enumerate(trigger_count)
        scatter!(
            ax,
            trigger_times[trigger_values.==key],
            repeat([unique], length(trigger_values[trigger_values.==key])),
            label = "$key: $(string(value))",
        )
    end
    fig[1, 2] = Legend(fig, ax)
    ax.ylabel = "Trigger Value"
    ax.xlabel = "Time (S)"
    Legend(fig[1, 2], ax)
    display(fig)
    return fig, ax
end




#########################################
# 2D head shape
# Base method without neighbors
function plot_layout_2d(fig, ax, layout; head_kwargs = Dict(), point_kwargs = Dict(), label_kwargs = Dict())
    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end

    # Default kwargs
    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => true, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)
    plot_points = pop!(point_kwargs, :plot_points)

    label_default_kwargs = Dict(:plot_labels => true, :fontsize => 20, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)
    plot_labels = pop!(label_kwargs, :plot_labels)
    xoffset = pop!(label_kwargs, :xoffset)
    yoffset = pop!(label_kwargs, :yoffset)

    # Head shape
    radius = 88 # mm
    arc!(ax, Point2f(0), radius * 2, -π, π; head_kwargs...) # head
    arc!(Point2f(radius * 2, 0), radius * 2 / 7, -π / 2, π / 2; head_kwargs...) # ear right
    arc!(Point2f(-radius * 2, 0), -radius * 2 / 7, π / 2, -π / 2; head_kwargs...) # ear left
    lines!(ax, Point2f[(-0.05, 0.5), (0.0, 0.6), (0.05, 0.5)] .* radius * 4; head_kwargs...) # nose

    # Regular points
    if plot_points
        scatter!(ax, layout[!, :x2], layout[!, :y2]; point_kwargs...)
    end

    if plot_labels
        for label in eachrow(layout)
            text!(ax, position = (label.x2 + xoffset, label.y2 + yoffset), label.label; label_kwargs...)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)

    return fig, ax
end


function plot_layout_2d(layout; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d(fig, ax, layout; kwargs...)
    display(fig)
    return fig, ax
end





# Method with interactive neighbor visualization
function plot_layout_2d(
    fig,
    ax,
    layout,
    neighbours::OrderedDict;
    head_kwargs = Dict(),
    label_kwargs = Dict(),
)
    # Draw head shape
    plot_layout_2d(
        fig,
        ax,
        layout;
        head_kwargs = head_kwargs,
        point_kwargs = Dict(:plot_points => false),
        label_kwargs = label_kwargs,
    )

    # Setup interactive points
    positions = Observable(Point2f.(layout.x2, layout.y2))
    base_size = 15
    hover_size = 25
    sizes = Observable(fill(base_size, length(layout.label)))

    # Add interactive scatter points
    p = scatter!(ax, positions; color = :black, markersize = sizes, inspectable = true, markerspace = :pixel)

    # Initialize line segments
    linesegments = Observable(Point2f[])
    lines!(ax, linesegments, color = :gray, linewidth = 3)

    # Add hover interaction
    on(events(fig).mouseposition) do mp
        plt, i = pick(fig)
        if plt == p
            # Reset all sizes to base size
            new_sizes = fill(base_size, length(layout.label))
            new_sizes[i] = hover_size
            sizes[] = new_sizes

            # Create lines to neighboring electrodes
            hovered_pos = positions[][i]
            new_lines = Point2f[]
            for neighbor in neighbours[Symbol(layout.label[i])].electrodes
                neighbor_idx = findfirst(==(neighbor), layout.label)
                push!(new_lines, hovered_pos, positions[][neighbor_idx])
            end
            linesegments[] = new_lines
        end
    end

    return fig, ax

end

function plot_layout_2d(layout, neighbours; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_2d(fig, ax, layout, neighbours; kwargs...)
    display(fig)
    return fig, ax
end


function graham_scan(points::Vector{Point2f})
    # Need at least 3 points for a hull
    length(points) < 3 && return points
    
    # Find leftmost point
    start = 1
    for i in 2:length(points)
        if points[i][1] < points[start][1]
            start = i
        end
    end
    p0 = points[start]
    
    # Sort by angle from p0
    other_points = vcat(points[1:start-1], points[start+1:end])
    sorted = sort(other_points, 
        by = p -> atan(p[2] - p0[2], p[1] - p0[1]))
    
    # Initialize hull with first point
    hull = [p0]
    
    # Build hull
    for p in sorted
        while length(hull) > 1
            v1 = hull[end] - hull[end-1]
            v2 = p - hull[end]
            cross = v1[1] * v2[2] - v1[2] * v2[1]
            if cross > 0  # Left turn
                break
            end
            pop!(hull)
        end
        push!(hull, p)
    end
    
    return hull
end


function point_border(xpos, ypos, border_size)
    # If only one electrode, return a complete circle
    if length(xpos) == 1
        circle_points = 0:2π/36:2π
        return [Point2f(xpos[1] + border_size * sin(θ), ypos[1] + border_size * cos(θ)) for θ in circle_points]
    end
    
    # Create points around each electrode position
    circle_points = 0:2π/36:2π
    border_points = Point2f[]
    for (x, y) in zip(xpos, ypos)
        for θ in circle_points
            push!(border_points, Point2f(
                x + border_size * sin(θ),
                y + border_size * cos(θ)
            ))
        end
    end
    
    # Get hull points
    hull_points = graham_scan(border_points)
    push!(hull_points, hull_points[1])  # Close the polygon
    
    return hull_points
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


function polyxy(p::Polygon)
    coords = GeoInterface.coordinates(GeoInterface.getexterior(p))
    return first.(coords), last.(coords)
   end




function add_topo_rois!(ax, layout, rois; border_size = 10, roi_kwargs = Dict())
    if (:x2 ∉ names(layout) || :y2 ∉ names(layout))
        polar_to_cartesian_xy!(layout)
    end
    roi_default_kwargs = Dict(:color => repeat([:black], length(rois)), :linewidth => repeat([2], length(rois)))
    roi_kwargs = merge(roi_default_kwargs, roi_kwargs)
    for (idx, roi) in enumerate(rois)
        xpos = filter(row -> row.label ∈ roi, layout).x2
        ypos = filter(row -> row.label ∈ roi, layout).y2
        border = point_border(xpos, ypos, border_size)
        lines!(ax, border, linewidth = roi_kwargs[:linewidth][idx], color = roi_kwargs[:color][idx])
    end
end

#########################################
# 3D head shape
function plot_layout_3d(fig, ax, layout; point_kwargs = Dict(), label_kwargs = Dict())

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

function plot_layout_3d(layout; kwargs...)
    fig = Figure()
    ax = Axis3(fig[1, 1])
    plot_layout_3d(fig, ax, layout; kwargs...)
    display(fig)
    return fig, ax
end


function plot_layout_3d(
    fig,
    ax,
    layout,
    neighbours::OrderedDict;
    head_kwargs = Dict(),
    label_kwargs = Dict(),
)
    # Draw head shape
    plot_layout_3d(
        fig,
        ax,
        layout;
        head_kwargs = head_kwargs,
        point_kwargs = Dict(:plot_points => false),
        label_kwargs = label_kwargs,
    )

    # Setup interactive points
    positions = Observable(Point2f.(layout.x2, layout.y2))
    base_size = 15
    hover_size = 25
    sizes = Observable(fill(base_size, length(layout.label)))

    # Add interactive scatter points
    p = scatter!(ax, positions; color = :black, markersize = sizes, inspectable = true, markerspace = :pixel)

    # Initialize line segments
    linesegments = Observable(Point2f[])
    lines!(ax, linesegments, color = :gray, linewidth = 3)

    # Add hover interaction
    on(events(fig).mouseposition) do mp
        plt, i = pick(fig)
        if plt == p
            # Reset all sizes to base size
            new_sizes = fill(base_size, length(layout.label))
            new_sizes[i] = hover_size
            sizes[] = new_sizes

            # Create lines to neighboring electrodes
            hovered_pos = positions[][i]
            new_lines = Point2f[]
            for neighbor in neighbours[Symbol(layout.label[i])].electrodes
                neighbor_idx = findfirst(==(neighbor), layout.label)
                push!(new_lines, hovered_pos, positions[][neighbor_idx])
            end
            linesegments[] = new_lines
        end
    end

    return fig, ax

end

function plot_layout_3d(layout, neighbours; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    plot_layout_3d(fig, ax, layout, neighbours; kwargs...)
    display(fig)
    return fig, ax
end





##########################################

function plot_topoplot(dat; kwargs...)
    plot_topoplot(dat, dat.layout; kwargs...)
end

function plot_topoplot(dat; kwargs...)
    plot_topoplot(dat.data, dat.layout; kwargs...)
end


# 2D topographic plot
function plot_topoplot(
    dat,
    layout;
    xlim = nothing,
    ylim = nothing,
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
        xlim = [dat.time[1], dat.time[end]]
        xlim_idx = 1:nrow(dat)
    end

    # convert xlim to index
    xlim_idx = find_idx_range(dat.time, xlim[1], xlim[2])

    # interpolate data
    data = data_interpolation_topo(
        mean.(eachcol(dat[xlim_idx, layout.label])),
        permutedims(Matrix(layout[!, [:x2, :y2]])),
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
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    display(GLMakie.Screen(), fig)
    return fig, ax

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




function plot_databrowser(
    dat::ContinuousData,
    channel_labels::Vector{<:AbstractString},
    ica::Union{InfoIca,Nothing} = nothing,
)


    function butterfly_plot(active::Bool)
        if active
            offset[] = LinRange(0, 0, length(channel_labels))
        else
            offset[] = LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, length(channel_labels))
        end
        for (_, label) in channel_data_labels
            label.visible = active
        end
    end

    function apply_lp_filter(active)
        clear_axes(ax, [channel_data_original, channel_data_labels])
        if active
            data = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
        else
            data = copy(dat.data)
        end
        draw(plot_labels = true)
    end


    function plot_extreme_lines(active)
        extreme_spans.visible = active
    end

    function toggle_button_group(fig, labels)
        # Pre-define all possible toggles with their functions
        toggle_configs = [
            ("Butterfly Plot", butterfly_plot),
            ("Trigger", plot_lines, "triggers"),
            ("vEOG", plot_lines, "is_vEOG"),
            ("hEOG", plot_lines, "is_hEOG"),
            ("extreme", plot_extreme_lines, "is_extreme"),
            ("LP-Filter On/Off", apply_lp_filter),
        ]

        # Filter toggles based on available labels
        toggles = [
            (config[1:2]..., Toggle(fig, active = false)) for
            config in toggle_configs if length(config) == 2 || config[3] in labels
        ]

        # Create grid layout
        grid = [[toggle[3], Label(fig, toggle[1], fontsize = 22, halign = :left), toggle[2]] for toggle in toggles]

        return permutedims(reduce(hcat, grid))  # Return as a proper matrix
    end

    # Makie Figure
    fig = Figure()
    ax = Axis(fig[1, 1])

    # controls
    # interactions(ax)
    deregister_interaction!(ax, :rectanglezoom)

    # Add these new variables for selection
    selection_active = Observable(false)
    selection_bounds = Observable((0.0, 0.0))  # (start, end)
    selection_visible = Observable(false)
    selection = poly!(ax, Point2f[], color = (:blue, 0.3))

    # Function to clear selection
    function clear_selection()
        selection[1] = Point2f[]
        selection_bounds[] = (0.0, 0.0)
        selection_visible[] = false
    end

    # Function to get precise mouse position
    function get_mouse_pos(scene_pos)
        # Convert from screen to data coordinates
        pos = to_world(ax.scene, scene_pos)
        return Float64(pos[1])
    end

    # Function to update selection rectangle
    function update_selection_rectangle(x1, x2)
        ylims = ax.limits[][2]
        selection[1] = Point2f[
            Point2f(Float64(x1), Float64(ylims[1])),
            Point2f(Float64(x2), Float64(ylims[1])),
            Point2f(Float64(x2), Float64(ylims[2])),
            Point2f(Float64(x1), Float64(ylims[2])),
        ]
    end

    # Modify mouse event handler
    on(events(ax).mousebutton) do event
        # Check if mouse is within axis bounds
        pos = events(ax).mouseposition[]
        bbox = ax.layoutobservables.computedbbox[]

        # Check if position is within bounding box
        if bbox.origin[1] <= pos[1] <= (bbox.origin[1] + bbox.widths[1]) &&
           bbox.origin[2] <= pos[2] <= (bbox.origin[2] + bbox.widths[2])
            world_pos = get_mouse_pos(pos)
            if event.button == Mouse.left
                if event.action == Mouse.press
                    if selection_visible[] &&
                       world_pos >= min(selection_bounds[][1], selection_bounds[][2]) &&
                       world_pos <= max(selection_bounds[][1], selection_bounds[][2])
                        clear_selection()
                    else
                        selection_active[] = true
                        selection_bounds[] = (world_pos, world_pos)
                        update_selection_rectangle(world_pos, world_pos)
                    end
                elseif event.action == Mouse.release && selection_active[]
                    selection_active[] = false
                    selection_visible[] = true
                    selection_bounds[] = (selection_bounds[][1], world_pos)
                    update_selection_rectangle(selection_bounds[][1], world_pos)
                end
            elseif event.button == Mouse.right && event.action == Mouse.press
                if selection_visible[] &&
                   world_pos >= min(selection_bounds[][1], selection_bounds[][2]) &&
                   world_pos <= max(selection_bounds[][1], selection_bounds[][2])
                    show_menu()
                end
            end
        end
    end

    # Update selection rectangle while dragging
    on(events(ax).mouseposition) do pos
        if selection_active[]
            world_pos = get_mouse_pos(pos)
            update_selection_rectangle(selection_bounds[][1], world_pos)
        end
    end

    # Function to get data in selected region
    function get_selected_data()
        x_min, x_max = minmax(selection_bounds[]...)
        time_mask = (x_min .<= dat.data.time .<= x_max)
        selected_data = dat.data[time_mask, :]
        println(
            "Selected data: $(round(x_min, digits = 2)) to $(round(x_max, digits = 2)) S, size $(size(selected_data))",
        )
        return selected_data
    end

    # Function to create menu window
    function show_menu()
        menu_fig = Figure()
        menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(["Topoplot"])]
        for btn in menu_buttons
            on(btn.clicks) do n
                selected_data = get_selected_data()
                if btn.label[] == "Topoplot"
                    plot_topoplot(selected_data, dat.layout)
                end
            end
        end
        display(GLMakie.Screen(), menu_fig)
    end

    # data to plot
    data = copy(dat.data)

    channel_data_original = Dict()
    channel_data_labels = Dict()

    # default xrange/yrange
    xlimit = 5000
    xrange = Observable(1:xlimit)
    yrange = Observable(-1500:1500)
    nchannels = length(channel_labels)
    channel_labels_original = channel_labels

    if nchannels > 1
        offset = Observable((LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, nchannels + 2)[2:end-1]))
    else # just centre
        offset = Observable(zeros(length(channel_labels)))
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

        clear_axes(ax, [channel_data_original, channel_data_labels])

        data = copy(dat.data)
        if nchannels > 1
            offset = LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels + 2)[2:end-1]
        else # just centre
            offset = zeros(nchannels)
        end
        draw(plot_labels = true)
    end

    # menu for ica selection
    menu_ica = nothing
    removed_activations = nothing
    components_to_remove = nothing
    components_removed = nothing
    if !isnothing(ica)
        menu_ica = hcat(
            Menu(fig, options = vcat(["None"], ica.ica_label), default = "None", direction = :down, fontsize = 18),
            Label(fig, "ICA Components", fontsize = 22, halign = :left),
        )
        on(menu_ica[1].selection) do s
            clear_axes(ax, [channel_data_original, channel_data_labels])
            components_to_remove = extract_int(s)
            if !isnothing(components_to_remove)
                if !isnothing(components_removed)
                    # go back to original data
                    data = restore_original_data(data, ica_result, [components_removed], removed_activations)
                end
                data, removed_activations = remove_ica_components(data, ica, [components_to_remove])
                components_removed = components_to_remove
            else
                data = restore_original_data(data, ica_result, [components_removed], removed_activations)
            end
            draw(plot_labels = true)
        end
    end


    slider_extreme = Slider(fig[1, 2], range = 0:5:100, startvalue = 100, width = 100)
    slider_lp_filter = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)

    crit_val = lift(slider_extreme.value) do x
        x
    end

    slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = xlimit, snap = true)
    slider_x = Slider(fig[2, 1], range = 1:50:nrow(data), startvalue = 1, snap = true)

    on(slider_x.value) do x
        new_range = x:min(nrow(data), (x + slider_range.value.val) - 1)
        if length(new_range) > 1
            xrange[] = new_range
        end
    end

    on(slider_range.value) do x
        new_range = slider_x.value.val:min(nrow(data), x + slider_x.value.val)
        if length(new_range) > 1
            xrange[] = new_range
        end
    end

    # keyboard events
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat)
            if selection_visible[]
                # Move selection region
                if event.key == Keyboard.left
                    width = selection_bounds[][2] - selection_bounds[][1]
                    new_start = max(data.time[1], selection_bounds[][1] - width / 5)
                    selection_bounds[] = (new_start, new_start + width)
                    update_selection_rectangle(selection_bounds[][1], selection_bounds[][2])
                elseif event.key == Keyboard.right
                    width = selection_bounds[][2] - selection_bounds[][1]
                    new_start = min(data.time[end] - width, selection_bounds[][1] + width / 5)
                    selection_bounds[] = (new_start, new_start + width)
                    update_selection_rectangle(selection_bounds[][1], selection_bounds[][2])
                end
            else
                # left/right x-range
                event.key == Keyboard.left && xback!(ax, xrange, data)
                event.key == Keyboard.right && xforward!(ax, xrange, data)
            end

            # up/down y-range
            event.key == Keyboard.down && yless!(ax, yrange)
            event.key == Keyboard.up && ymore!(ax, yrange)

        end

    end

    # position GUI controls
    if isnothing(menu_ica)
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
    else
        fig[1, 2] = grid!(
            vcat(
                toggles[:, 1:2],
                hcat(
                    slider_lp_filter,
                    Label(fig, @lift("LP-Filter: $($(slider_lp_filter.value)) Hz"), fontsize = 22, halign = :left),
                ),
                hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)),
                menu,
                menu_ica,
            ),
            tellheight = false,
        )

    end
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


    # Version 1:
    function draw(; plot_labels = true)
        for (idx, col) in enumerate(channel_labels)
            channel_data_original[col] = lines!(
                ax,
                data[!, :time],
                @lift(data[!, $col] .+ $offset[idx]),  # Make y-values reactive to offset
                color = @lift(abs.(data[!, col]) .>= $crit_val),
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = 2,
            )
            if plot_labels
                channel_data_labels[col] = text!(
                    ax,
                    @lift(data[$xrange, :time][1]),
                    @lift(data[$xrange, $col][1] .+ $offset[idx]),  # Make label position reactive to offset
                    text = col,
                    align = (:left, :center),
                    fontsize = 18,
                )
            end
        end
    end

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    hideydecorations!(ax, label = true)
    draw(plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

# Convenience methods
plot_databrowser(dat::ContinuousData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::ContinuousData, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
plot_databrowser(dat::ContinuousData, channel_label::AbstractString) = plot_databrowser(dat, [channel_label])
plot_databrowser(dat::ContinuousData, channel_label::AbstractString, ica::InfoIca) =
    plot_databrowser(dat, [channel_label], ica)


###########################################################

function plot_databrowser(
    dat::EpochData,
    channel_labels::Vector{<:AbstractString},
    ica::Union{InfoIca,Nothing} = nothing,
)

    function butterfly_plot(active::Bool)
        if active
            # Set all offsets to 0 for butterfly plot
            offset[] = LinRange(0, 0, length(channel_labels))
            # Hide labels
            for (_, label) in channel_data_labels
                label.visible = false
            end
        else
            # Reset to original spacing
            offset[] = LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, length(channel_labels))
            # Show labels
            for (_, label) in channel_data_labels
                label.visible = true
            end
        end
    end


    function apply_lp_filter(active)
        clear_axes(ax, [channel_data_original, channel_data_labels])
        if active
            data = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
        else
            data = copy(dat.data)
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
        # Pre-define all possible toggles with their functions
        toggle_configs = [
            ("Butterfly Plot", butterfly_plot),
            ("Trigger", plot_lines, "triggers"),
            ("vEOG", plot_lines, "is_vEOG"),
            ("hEOG", plot_lines, "is_hEOG"),
            ("extreme", plot_extreme_lines, "is_extreme"),
            ("LP-Filter On/Off", apply_lp_filter),
        ]

        # Filter toggles based on available labels
        toggles = [
            (config[1:2]..., Toggle(fig, active = false)) for
            config in toggle_configs if length(config) == 2 || config[3] in labels
        ]

        # Create grid layout
        grid = [[toggle[3], Label(fig, toggle[1], fontsize = 22, halign = :left), toggle[2]] for toggle in toggles]

        return permutedims(reduce(hcat, grid))  # Return as a proper matrix
    end

    function step_epoch_forward()
        clear_axes(ax, [channel_data_original, channel_data_labels])
        trial[] = min(length(dat.data), trial.val[1] + 1)
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        menu_trial.i_selected[] = trial
        update_extreme_spans!()
        update_markers!(markers)
        draw()
    end

    function step_epoch_backward()
        clear_axes(ax, [channel_data_original, channel_data_labels])
        trial[] = max(1, trial.val[1] - 1)
        ax.title = "Epoch $(trial.val)/$(length(dat.data))"
        menu_trial.i_selected[] = trial
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
    channel_data_labels = Dict()

    # default xrange/yrange
    xlimit = nrow(dat.data[1])
    xrange = Observable(1:xlimit)
    trial = Observable(1) # first trial
    yrange = Observable(-1500:1500)
    nchannels = length(channel_labels)
    channel_labels_original = channel_labels

    if nchannels > 1
        offset = LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, nchannels + 2)[2:end-1]
    else # just centre
        offset = zeros(length(channel_labels))
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

        clear_axes(ax, [channel_data_original, channel_data_labels])

        data = copy(dat.data)
        if nchannels > 1
            offset = LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels + 2)[2:end-1]
        else # just centre
            offset = zeros(nchannels)
        end
        draw(plot_labels = true)


    end

    menu_trial = hcat(
        Menu(fig, options = 1:length(data), default = 1, direction = :down, fontsize = 18),
        Label(fig, "Epoch", fontsize = 22, halign = :left),
    )
    on(menu_trial[1].selection) do s
        clear_axes(ax, [channel_data_original, channel_data_labels])
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
        # Create one plot per channel and update its data
        for (idx, col) in enumerate(channel_labels)
            # Initial empty plot with label
            line = series!(
                ax,
                Point2f[],  # Will be updated by the lift
                color = :darkgrey,
                linewidth = 2,
                label = col,  # Add label to the line
                markersize = 0,  # Hide markers
            )
            channel_data_original[col] = line

            # Update plot data when range changes
            lift(xrange, crit_val) do current_range, cv
                points = [Point2f(data[i, :time], data[i, col] + offset[idx]) for i in current_range]
                line[1][] = points
            end
        end

        # Update xlims separately
        lift(xrange) do current_range
            xlims!(ax, data[current_range[1], :time], data[current_range[end], :time])
        end
    end

    # plot theme adjustments
    fontsize_theme = Theme(fontsize = 24)
    update_theme!(fontsize_theme)

    hideydecorations!(ax, label = true)
    draw(plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

# Convenience methods
plot_databrowser(dat::EpochData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::EpochData, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
plot_databrowser(dat::EpochData, channel_label::AbstractString) = plot_databrowser(dat, [channel_label])
plot_databrowser(dat::EpochData, channel_label::AbstractString, ica::InfoIca) =
    plot_databrowser(dat, [channel_label], ica)


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
function plot_erp(
    fig,
    ax,
    dat::ErpData,
    channels::Union{Vector{<:AbstractString},Vector{Symbol}};
    average_channels = false,
    kwargs = Dict(),
)

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

    if kwargs[:add_topoplot]
        # just put in top left
        topo_ax =
            Axis(fig[1, kwargs[:topoplot_fig]], width = Relative(0.2), height = Relative(0.2), halign = 0, valign = 1)
        layout = filter(row -> row.label in String.(channels), erp.layout)
        if average_channels
            head_shape_2d(fig, topo_ax, layout, point_kwargs = Dict(:color => kwargs[:color], :markersize => 18))
        else
            head_shape_2d(
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

function plot_erp(
    dat::ErpData,
    channels::Union{Vector{<:AbstractString},Vector{Symbol}};
    average_channels = false,
    kwargs = Dict(),
)
    fig = Figure()
    ax = Axis(fig[1, 1])
    fig, ax = plot_erp(fig, ax, dat, channels; average_channels = average_channels, kwargs = kwargs)
    return fig, ax
end

function plot_erp(dat::ErpData, channels::Union{AbstractString,Symbol}; average_channels = false, kwargs...)
    plot_erp(dat, [channels], average_channels; kwargs...)
end

function plot_erp(dat::ErpData; average_channels = false, kwargs...)
    plot_erp(dat, dat.layout.label; average_channels = average_channels, kwargs...)
end



function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData, channels; average_channels = false, kwargs = Dict())
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    plot_erp(fig, ax1, dat_orig, channels; average_channels = average_channels, kwargs = kwargs)
    ax2 = Axis(fig[2, 1])
    kwargs = merge(kwargs, kwargs)
    kwargs[:topoplot_fig] = 2
    plot_erp(fig, ax2, dat_cleaned, channels; average_channels = average_channels, kwargs = kwargs)
    linkaxes!(ax1, ax2)
    display(fig)
end


function plot_erp(dat_orig::ErpData, dat_cleaned::ErpData; average_channels = false, kwargs...)
    plot_erp(dat_orig, dat_cleaned, dat_orig.layout.label; average_channels = average_channels, kwargs...)
end

function best_rect(n)
    dim1 = ceil(Int, sqrt(n))
    dim2 = ceil(Int, n ./ dim1)
    return [dim1, dim2]
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
    row_names = corr_df[!, :row]
    col_names = names(corr_df)[2:end]

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
        ax.title = ica.ica_label[comps[count]]
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
        head_shape_2d(
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
    head_shape_2d(fig, ax, layout, head_kwargs = head_kwargs, point_kwargs = point_kwargs, label_kwargs = label_kwargs)
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










