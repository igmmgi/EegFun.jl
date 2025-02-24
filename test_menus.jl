using GLMakie
using GeometryBasics: Point2f

# Sample data
x = 0:0.1:10
y = sin.(x)

# Create main figure
fig = Figure()
ax = Axis(fig[1, 1], title = "Original Plot")
deactivate_interaction!(ax, :rectanglezoom)
lines!(ax, x, y)

# Variables for selection
selection_active = Observable(false)
x_start = Observable(0.0)
x_end = Observable(0.0)

# Add a rectangle for the selection
selection = poly!(ax, Point2f[], color = (:blue, 0.3))

# Function to get data in selected region
function get_selected_data()
    x_min, x_max = minmax(x_start[], x_end[])
    mask = x_min .<= x .<= x_max
    return x[mask], y[mask]
end

# Function to create menu window
function show_menu()
    menu_fig = Figure(size=(200, 150))
    menu_buttons = [
        Button(menu_fig[i, 1], label = label)
        for (i, label) in enumerate(["Scatter Plot", "Histogram", "Line Fit"])
    ]
    
    for btn in menu_buttons
        on(btn.clicks) do n
            x_sel, y_sel = get_selected_data()
            
            # Create new figure for the plot
            f = Figure()
            ax2 = Axis(f[1,1])
            
            if btn.label[] == "Scatter Plot"
                scatter!(ax2, x_sel, y_sel)
                ax2.title = "Scatter Plot of Selected Region"
            elseif btn.label[] == "Histogram"
                hist!(ax2, y_sel, bins=20)
                ax2.title = "Histogram of Selected Region"
            elseif btn.label[] == "Line Fit"
                scatter!(ax2, x_sel, y_sel)
                # Simple linear fit
                coeffs = [ones(length(x_sel)) x_sel] \ y_sel
                fit_y = coeffs[1] .+ coeffs[2] .* x_sel
                lines!(ax2, x_sel, fit_y, color=:red)
                ax2.title = "Linear Fit of Selected Region"
            end
            
            display(GLMakie.Screen(), f)
            println("Created: $(btn.label[])")
        end
    end
    
    display(GLMakie.Screen(), menu_fig)
end

# Handle mouse events
on(events(ax).mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press
            selection_active[] = true
            pos = to_world(ax.scene, events(ax).mouseposition[])
            x_start[] = pos[1]
        else
            selection_active[] = false
            pos = to_world(ax.scene, events(ax).mouseposition[])
            x_end[] = pos[1]
        end
    elseif event.button == Mouse.right && event.action == Mouse.press
        pos = to_world(ax.scene, events(ax).mouseposition[])[1]
        if pos >= min(x_start[], x_end[]) && pos <= max(x_start[], x_end[])
            show_menu()
        end
    end
end

# Update selection rectangle while dragging
on(events(ax).mouseposition) do pos
    if selection_active[]
        world_pos = to_world(ax.scene, pos)
        ylims = ax.finallimits[]
        selection[1] = [
            Point2f(x_start[], ylims.origin[2]),
            Point2f(world_pos[1], ylims.origin[2]),
            Point2f(world_pos[1], ylims.origin[2] + ylims.widths[2]),
            Point2f(x_start[], ylims.origin[2] + ylims.widths[2])
        ]
    end
end

# Display the main figure
display(GLMakie.Screen(), fig)