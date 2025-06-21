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
    on_click = nothing,  # Callback function for when a plot is clicked
)

    # x/y limits
    isnothing(xlim) && (xlim = data_limits_x(dat.data))
    isnothing(ylim) && (ylim = data_limits_y(dat.data, dat.layout.label))

    if (:x2 ∉ propertynames(dat.layout) || :y2 ∉ propertynames(dat.layout))
        polar_to_cartesian_xy!(dat.layout)
    end

    if plot_label && isnothing(plot_label_position)
        plot_label_positon = [xlim[1] ylim[2]]
    end

    # Calculate electrode positions in data coordinate space
    xminmaxrange, yminmaxrange = datarange(dat.layout.x2), datarange(dat.layout.y2)
    # Scale the positions to match the data coordinate space while maintaining centering
    xpositions = (dat.layout.x2 ./ xminmaxrange) .* (xlim[2] - xlim[1]) .+ (xlim[1] + (xlim[2] - xlim[1])/2)
    ypositions = (dat.layout.y2 ./ yminmaxrange) .* (ylim[2] - ylim[1]) .+ (ylim[1] + (ylim[2] - ylim[1])/2)
    fig = Figure()
    
    # Store normalized positions and their corresponding channels
    position_channel_map = Dict()
    
    for (x, y, label) in zip(xpositions, ypositions, dat.layout.label)
        # Convert position to normalized coordinates for axis placement
        norm_x = (x - xlim[1]) / (xlim[2] - xlim[1])
        norm_y = (y - ylim[1]) / (ylim[2] - ylim[1])
        ax = Axis(fig[1, 1], width = Relative(plot_width), height = Relative(plot_height), halign = norm_x, valign = norm_y)
        position_channel_map[(norm_x, norm_y)] = label  # Store normalized coordinates
        
        lines!(ax, dat.data[!, :time], dat.data[!, label])
        vlines!(ax, [0], color = :black)
        hlines!(ax, [0], color = :black)
        if plot_label
            text!(
                ax,
                plot_label_positon[1],
                plot_label_positon[end],
                fontsize = plot_label_fontsize,
                text = String(label),
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

    # Add click handler to the main axis
    if !isnothing(on_click)
        # Get the main axis (the one in the grid cell)
        main_ax = fig[1, 1]
        
        # Variables to track selection
        selection_start = nothing
        selection_rects = []  # Store all selection rectangles
        selection_coords_list = []  # Store normalized coordinates of all selection rectangles
        
        # Function to convert screen coordinates to normalized coordinates
        function screen_to_normalized(pos, size)
            return (pos[1] / size[1], pos[2] / size[2])
        end
        
        # Function to check if a point is inside any selection rectangle
        function is_point_in_any_rect(point, rects)
            for rect in rects
                if is_point_in_rect(point, rect)
                    return true
                end
            end
            return false
        end
        
        # Function to check if a point is inside the selection rectangle
        function is_point_in_rect(point, rect)
            x, y = point
            x1, y1 = rect[1]
            x2, y2 = rect[2]
            return min(x1, x2) <= x <= max(x1, x2) && min(y1, y2) <= y <= max(y1, y2)
        end
        
        # Function to draw selection rectangle
        function update_selection_rect(start_pos, current_pos)
            if !isnothing(selection_start)
                fig_size = size(fig.scene)
                start_norm = screen_to_normalized(start_pos, fig_size)
                current_norm = screen_to_normalized(current_pos, fig_size)
                
                # Draw the selection rectangle with a more visible style
                rect = poly!(
                    fig.scene,
                    [start_norm, (current_norm[1], start_norm[2]), current_norm, (start_norm[1], current_norm[2])],
                    color = (:blue, 0.3),  # More opaque blue
                    strokecolor = :red,    # Red border
                    strokewidth = 2,       # Thicker border
                    overdraw = true,       # Ensure it's drawn on top
                    space = :relative      # Use relative coordinates
                )
                
                # Store the rectangle and its coordinates
                push!(selection_rects, rect)
                push!(selection_coords_list, (start_norm, current_norm))
            end
        end
        
        # Mouse button press - start selection
        on(events(fig).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                # If there are existing selections, check if click is inside any of them
                if !isempty(selection_coords_list)
                    click_pos = events(fig).mouseposition[]
                    fig_size = size(fig.scene)
                    click_norm = screen_to_normalized(click_pos, fig_size)
                    
                    # Check if click is inside any selection rectangle
                    if is_point_in_any_rect(click_norm, selection_coords_list)
                        # Select only electrodes that are inside any of the selection rectangles
                        selected_channels = Symbol[]
                        for ((x, y), channel) in position_channel_map
                            if is_point_in_any_rect((x, y), selection_coords_list)
                                push!(selected_channels, Symbol(channel))
                            end
                        end
                        
                        # Call the callback with selected channels
                        if !isempty(selected_channels)
                            on_click(selected_channels)
                        end
                        
                        # Clean up all selection rectangles
                        for rect in selection_rects
                            delete!(fig.scene, rect)
                        end
                        empty!(selection_rects)
                        empty!(selection_coords_list)
                        return
                    end
                end
                
                # Start new selection
                selection_start = events(fig).mouseposition[]
                update_selection_rect(selection_start, selection_start)
            elseif event.button == Mouse.left && event.action == Mouse.release
                if !isnothing(selection_start)
                    current_pos = events(fig).mouseposition[]
                    update_selection_rect(selection_start, current_pos)
                    selection_start = nothing
                end
            end
        end
        
        # Mouse movement - update selection rectangle
        on(events(fig).mouseposition) do pos
            if !isnothing(selection_start)
                # Remove the last rectangle if it exists
                if !isempty(selection_rects)
                    delete!(fig.scene, pop!(selection_rects))
                    pop!(selection_coords_list)
                end
                update_selection_rect(selection_start, pos)
            end
        end
    end
    
    linkaxes!(filter(x -> x isa Axis, fig.content)...)
    display(GLMakie.Screen(), fig)
    return fig, position_channel_map
end

