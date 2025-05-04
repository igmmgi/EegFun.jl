function plot_ica_topoplot_single(
    fig,
    position,
    ica,
    comp_idx,
    layout;
    colorbar_kwargs = Dict(),
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
)

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30)
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)

    # Create a consistent GridLayout at this position - always the same structure
    gl = fig[position...] = GridLayout()

    # Create main axis in the first column
    ax = Axis(gl[1, 1], title = @sprintf("IC %d (%.1f%%)", comp_idx, ica.variance[comp_idx] * 100))

    # Extract layout data
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]

    # Create the topo data
    data = data_interpolation_topo(ica.mixing[:, comp_idx], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)

    # Create the plot
    radius = 88 # mm
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data;
        colormap = topo_kwargs[:colormap],
        nan_color = :transparent,
    )

    # Add head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    # Hide decorations for cleaner look
    hidexdecorations!(ax, grid = false)
    hideydecorations!(ax, grid = false)

    # Add colorbar if requested
    if plot_colorbar
        # Default values
        width = get(colorbar_kwargs, :width, 10)
        height = get(colorbar_kwargs, :height, Relative(0.8))
        ticksize = get(colorbar_kwargs, :ticklabelsize, 10)

        # Create colorbar in second column
        Colorbar(gl[1, 2], co; width = width, height = height, ticklabelsize = ticksize)

        # Set column sizes
        colsize!(gl, 1, Relative(0.85))  # Plot column
        colsize!(gl, 2, Relative(0.15))  # Colorbar column
    else
        # If no colorbar, make the plot take the full width
        colsize!(gl, 1, Relative(1.0))
    end

    return ax, co
end

function plot_ica_topoplot(
    ica,
    layout;
    comps = nothing,
    dims = nothing,
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
    use_global_scale = false,
)

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30, :colorbar_plot_numbers => [])
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)
    colorbar_plot_numbers = pop!(colorbar_kwargs, :colorbar_plot_numbers)

    # Process inputs
    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end

    if isnothing(comps)
        comps = 1:size(ica.mixing)[2]
    end

    # Create figure with reduced margins
    fig = Figure(
        figure_padding = 0,
        backgroundcolor = :white
    )

    # Calculate layout dimensions
    if isnothing(dims)
        dims = best_rect(length(comps))
    end

    # Extract layout data once for all plots
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]
    radius = 88

    # Calculate all data first
    all_data = []
    for i in eachindex(comps)
        data = data_interpolation_topo(ica.mixing[:, comps[i]], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)
        push!(all_data, data)
    end

    # Calculate levels based on use_global_scale
    levels = nothing
    if use_global_scale
        # Find global min/max, safely handling NaN values
        all_values = [v for data in all_data for v in data if !isnan(v)]

        # If we have valid values, use them; otherwise use defaults
        if !isempty(all_values)
            data_min = minimum(all_values)
            data_max = maximum(all_values)
        else
            data_min = -1.0
            data_max = 1.0
        end

        # Ensure min != max to avoid range error
        if data_min == data_max
            data_min -= 0.1
            data_max += 0.1
        end

        # Create consistent levels for all plots
        levels = range(data_min, data_max, length = 20)
    end

    # First create all the GridLayouts for consistent sizing
    grids = []
    for i in eachindex(comps)
        # Calculate grid position
        row = ceil(Int, i / dims[2])
        col = ((i - 1) % dims[2]) + 1

        # Create a GridLayout for each plot
        gl = fig[row, col] = GridLayout()
        push!(grids, gl)
    end

    # Now add the plots and colorbars
    for i in eachindex(comps)
        gl = grids[i]

        # Create main axis in the first column
        ax = Axis(gl[1, 1], title = @sprintf("IC %d (%.1f%%)", comps[i], ica.variance[comps[i]] * 100))

        # Get pre-calculated data
        data = all_data[i]

        # Create the plot with either global or individual levels
        kwargs = Dict{Symbol,Any}(:colormap => :jet, :nan_color => :transparent)
        
        # Calculate levels for this specific plot if not using global scale
        if !use_global_scale
            # Find min/max for this specific data, safely handling NaN values
            plot_values = [v for v in data if !isnan(v)]
            if !isempty(plot_values)
                data_min = minimum(plot_values)
                data_max = maximum(plot_values)
            else
                data_min = -1.0
                data_max = 1.0
            end

            # Ensure min != max to avoid range error
            if data_min == data_max
                data_min -= 0.1
                data_max += 0.1
            end

            # Create levels for this specific plot
            kwargs[:levels] = range(data_min, data_max, length = 20)
        else
            # Use the global levels we calculated earlier
            kwargs[:levels] = levels
        end

        co = contourf!(
            ax,
            range(-radius * 2, radius * 2, length = gridscale),
            range(-radius * 2, radius * 2, length = gridscale),
            data;
            kwargs...,
        )

        # Add head shape
        plot_layout_2d!(
            fig,
            ax,
            layout,
            head_kwargs = head_kwargs,
            point_kwargs = point_kwargs,
            label_kwargs = label_kwargs,
        )

        # Hide decorations for cleaner look
        hidexdecorations!(ax, grid = false)
        hideydecorations!(ax, grid = false)

        # Set up colorbar space for ALL plots for consistent sizing
        if plot_colorbar
            # Default values
            width = get(colorbar_kwargs, :width, 10)
            height = get(colorbar_kwargs, :height, Relative(0.8))
            ticksize = get(colorbar_kwargs, :ticklabelsize, 10)

            # Create colorbar or placeholder in second column
            if i in colorbar_plot_numbers
                # Only create visible colorbar for the first plot
                Colorbar(gl[1, 2], co; width = width, height = height, ticklabelsize = ticksize)
            else
                # For other plots, create an empty colorbar placeholder
                placeholder = Axis(
                    gl[1, 2];
                    width = width,
                    height = height,
                    rightspinevisible = false,
                    leftspinevisible = false,
                    topspinevisible = false,
                    bottomspinevisible = false,
                    xticklabelsvisible = false,
                    yticklabelsvisible = false,
                    xticksvisible = false,
                    yticksvisible = false,
                    xgridvisible = false,
                    ygridvisible = false,
                )
                # Hide everything to make it truly empty
                hidexdecorations!(placeholder)
                hideydecorations!(placeholder)
            end

            # Set column sizes after creating the colorbar column
            colsize!(gl, 1, Relative(0.85))  # Plot column
            colsize!(gl, 2, Relative(0.15))  # Colorbar column
        else
            # If no colorbar, make the plot take the full width
            colsize!(gl, 1, Relative(1.0))
        end
    end

    # Apply consistent sizing to all grid rows and columns
    # This ensures alignment across the entire figure
    for i = 1:dims[1]
        rowsize!(fig.layout, i, Relative(1 / dims[1]))
    end

    for j = 1:dims[2]
        colsize!(fig.layout, j, Relative(1 / dims[2]))
    end

    return fig
end

# Version for the component viewer to use
function plot_ica_topoplot(
    fig,
    ax,
    ica,
    comp_idx,
    layout;
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict(),
    topo_kwargs = Dict(),
    colorbar_kwargs = Dict(),
)

    head_default_kwargs = Dict(:color => :black, :linewidth => 2)
    head_kwargs = merge(head_default_kwargs, head_kwargs)

    point_default_kwargs = Dict(:plot_points => false, :marker => :circle, :markersize => 12, :color => :black)
    point_kwargs = merge(point_default_kwargs, point_kwargs)

    label_default_kwargs = Dict(:plot_labels => false, :fontsize => 20, :color => :black, :color => :black, :xoffset => 0, :yoffset => 0)
    label_kwargs = merge(label_default_kwargs, label_kwargs)

    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    gridscale = pop!(topo_kwargs, :gridscale)

    colorbar_default_kwargs = Dict(:plot_colorbar => true, :width => 30, :colorbar_plot_numbers => [])
    colorbar_kwargs = merge(colorbar_default_kwargs, colorbar_kwargs)
    plot_colorbar = pop!(colorbar_kwargs, :plot_colorbar)
    colorbar_plot_numbers = pop!(colorbar_kwargs, :colorbar_plot_numbers)

    # Clear the axis before drawing
    empty!(ax)

    # Extract layout data
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]

    # Create the topo data
    data = data_interpolation_topo(ica.mixing[:, comp_idx], permutedims(Matrix(tmp_layout[!, [:x2, :y2]])), gridscale)

    # Create the plot directly on the existing axis
    radius = 88 # mm
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data;
        colormap = topo_kwargs[:colormap],
        nan_color = :transparent,
    )

    # Add head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )

    # Set title - for component viewer
    ax.title = string(ica.ica_label[comp_idx])

    # Hide decorations for cleaner look
    hidexdecorations!(ax, grid = false)
    hideydecorations!(ax, grid = false)

    return co
end


# Create a state structure to hold the visualization state
mutable struct IcaComponentState
    # Data
    dat::ContinuousData
    ica_result::InfoIca
    components::Matrix{Float64}
    total_components::Int

    # View settings
    n_visible_components::Int
    window_size::Int

    # Observables
    comp_start::Observable{Int}
    xrange::Observable{UnitRange{Int}}
    ylims::Observable{Tuple{Float64,Float64}}
    channel_data::Observable{Vector{Float64}}
    show_channel::Observable{Bool}
    channel_yscale::Observable{Float64}
    use_global_scale::Observable{Bool}
    invert_scale::Observable{Bool}  # New observable for scale inversion

    # New field for specific components
    specific_components::Union{Nothing,Vector{Int}}

    # Plot elements
    axs::Vector{Axis}
    channel_axs::Vector{Union{Axis,Nothing}}  # Allow for nothing values
    topo_axs::Vector{Axis}
    lines_obs::Vector{Observable{Vector{Float64}}}

    # New field for boolean indicators
    channel_bool_indicators::Dict{Int,Any}

    function IcaComponentState(dat, ica_result, n_visible_components, window_size, specific_components = nothing)
        # Prepare data matrix
        dat_matrix = prepare_ica_data_matrix(dat, ica_result)
        components = ica_result.unmixing * dat_matrix
        total_components = size(components, 1)

        # Create observables
        comp_start = Observable(1)
        use_global_scale = Observable(false)
        invert_scale = Observable(false)  # Initialize to false

        # Find index closest to time 0 to center the initial view
        time_zero_idx = findmin(abs.(dat.data.time))[2]
        half_window = div(window_size, 2)
        start_idx = max(1, time_zero_idx - half_window)
        end_idx = min(size(components, 2), start_idx + window_size - 1)

        # Adjust start_idx if end_idx reached the boundary
        if end_idx == size(components, 2)
            start_idx = max(1, end_idx - window_size + 1)
        end

        xrange = Observable(start_idx:end_idx)

        # Set initial range based on specific components if provided
        if !isnothing(specific_components)
            comps_to_use = specific_components
        else
            comps_to_use = 1:n_visible_components
        end

        # Calculate initial y-range based on components we'll show
        initial_range = maximum(abs.(extrema(components[comps_to_use, 1:window_size])))
        ylims = Observable((-initial_range, initial_range))
        channel_data = Observable(zeros(size(dat.data, 1)))
        show_channel = Observable(false)
        channel_yscale = Observable(1.0)

        # Initialize empty plot element arrays
        axs = Vector{Axis}()
        channel_axs = Vector{Union{Axis,Nothing}}()  # Allow for nothing values
        topo_axs = Vector{Axis}()
        lines_obs = Vector{Observable{Vector{Float64}}}()
        channel_bool_indicators = Dict{Int,Any}()

        new(
            dat,
            ica_result,
            components,
            total_components,
            n_visible_components,
            window_size,
            comp_start,
            xrange,
            ylims,
            channel_data,
            show_channel,
            channel_yscale,
            use_global_scale,
            invert_scale,
            specific_components,
            axs,
            channel_axs,
            topo_axs,
            lines_obs,
            channel_bool_indicators,
        )
    end
end

"""
    plot_ica_component_activation(dat::ContinuousData, ica_result::InfoIca; kwargs...)

Create an interactive visualization of ICA components with topographic maps and time series plots.

# Arguments
- `dat::ContinuousData`: Continuous EEG data
- `ica_result::InfoIca`: ICA result object

# Keyword Arguments
- `n_visible_components::Int=10`: Number of components visible at once
- `window_size::Int=2000`: Initial window size in samples
- `topo_kwargs::Dict=Dict()`: Additional keyword arguments for topoplots

# Returns
- `fig::Figure`: The Figure object containing the plot
"""
function plot_ica_component_activation(
    dat::ContinuousData,
    ica_result::InfoIca;
    n_visible_components::Int = 10,
    window_size::Int = 2000,
    topo_kwargs = Dict(),
    specific_components = nothing,
)
    # Create state, using specific components if provided
    if !isnothing(specific_components)
        # Make sure n_visible_components matches the length of specific_components
        n_visible_components = length(specific_components)
    end

    state = IcaComponentState(dat, ica_result, n_visible_components, window_size, specific_components)

    # Create figure with reduced margins
    fig = Figure(
        figure_padding = 0,
        backgroundcolor = :white
    )

    # Setup plots
    create_component_plots!(fig, state, topo_kwargs)

    # Add controls
    add_navigation_controls!(fig, state)

    # Add navigation sliders
    add_navigation_sliders!(fig, state)

    # Add channel menu
    add_channel_menu!(fig, state)

    # Add keyboard interactions
    setup_keyboard_interactions!(fig, state)

    # Set column sizes to give more space to the time series plots
    colsize!(fig.layout, 1, Relative(0.15))  # Topoplots - now narrower
    colsize!(fig.layout, 2, Relative(0.85))  # Time series - now wider

    # Reduce spacing between plot rows AND add space before controls
    rowgap!(fig.layout, 0) # Keep 0 gap between plot rows
    rowgap!(fig.layout, state.n_visible_components, 20) # Add 20px gap after last plot row

    display(fig)
    return fig
end

# Helper function to prepare data matrix for ICA
function prepare_ica_data_matrix(dat::ContinuousData, ica_result::InfoIca)
    dat_matrix = permutedims(Matrix(dat.data[!, ica_result.data_label]))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica_result.scale
    return dat_matrix
end

# Create a simpler override specifically for the component viewer
function plot_topoplot_in_viewer!(fig, topo_ax, ica_result, comp_idx, layout; use_global_scale = false, global_min = nothing, global_max = nothing)
    # Clear the axis
    empty!(topo_ax)

    # Extract layout data and prepare the plot
    tmp_layout = layout[(layout.label.∈Ref(ica_result.data_label)), :]

    # Create the topo data
    gridscale = 300
    data = data_interpolation_topo(
        ica_result.mixing[:, comp_idx],
        permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
        gridscale,
    )

    # Calculate levels based on use_global_scale
    levels = nothing
    if use_global_scale && !isnothing(global_min) && !isnothing(global_max)
        # Use the global min/max
        data_min = global_min
        data_max = global_max
    else
        # Find min/max for this specific data, safely handling NaN values
        plot_values = [v for v in data if !isnan(v)]
        if !isempty(plot_values)
            data_min = minimum(plot_values)
            data_max = maximum(plot_values)
        else
            data_min = -1.0
            data_max = 1.0
        end
    end

    # Ensure min != max to avoid range error
    if data_min == data_max
        data_min -= 0.1
        data_max += 0.1
    end


    # Plot directly on the provided axis
    radius = 88 # mm
    co = contourf!(
        topo_ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data;
        colormap = :jet,
        nan_color = :transparent, 
        levels = range(data_min, data_max, length=20)
    )

    # Add head shape with default settings
    plot_layout_2d!(fig, topo_ax, layout, label_kwargs = Dict(:plot_labels => false), point_kwargs = Dict(:plot_points => false))

    # Hide decorations
    hidexdecorations!(topo_ax, grid = false)
    hideydecorations!(topo_ax, grid = false)

    return co
end

# Update create_component_plots! to use our new function
function create_component_plots!(fig, state, topo_kwargs = Dict())
    # Create axes for the time series plots
    for i = 1:state.n_visible_components
        # Create topo plot axis first (now on the left)
        topo_ax = Axis(fig[i, 1])
        push!(state.topo_axs, topo_ax)

        # Get the actual component number
        comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
            state.specific_components[i]
        else
            state.comp_start[] + i - 1
        end

        # Time series axis creation (now on the right)
        ax = Axis(
            fig[i, 2],
            ylabel = @sprintf("IC %d", comp_idx),
            yaxisposition = :left,
            yticklabelsvisible = false,  # Hide y-axis tick labels
            yticksvisible = true,  # Keep the tick marks themselves
            xticklabelsvisible = (i == state.n_visible_components),  # Only show x-tick labels on last plot
            xticksvisible = (i == state.n_visible_components),  # Only show x-ticks on last plot
            xgridvisible = false,  # Hide x grid
            ygridvisible = false,  # Hide y grid
            xminorgridvisible = false,  # Hide x minor grid
            yminorgridvisible = false,  # Hide y minor grid
            bottomspinevisible = true,  # Show bottom spine
            topspinevisible = true,  # Show top spine
            rightspinevisible = false,  # Hide right spine
            leftspinevisible = true,  # Show left spine
            ylabelpadding = 0.0,  # Reduce y-label padding
            yticklabelpad = 0.0,  # Reduce y-tick label padding
            yticklabelspace = 0.0,  # Reduce y-tick label space
        )
        push!(state.axs, ax)

        # Always create channel overlay axis
        ax_channel = Axis(
            fig[i, 2],
            yticklabelsvisible = false,
            yticksvisible = false,
            yaxisposition = :right,
            xaxisposition = :top,
            xticklabelsvisible = false,  # Hide x-tick labels
            xticksvisible = false,  # Hide x-ticks
            xgridvisible = false,  # Hide x grid
            ygridvisible = false,  # Hide y grid
            xminorgridvisible = false,  # Hide x minor grid
            yminorgridvisible = false,  # Hide y minor grid
            bottomspinevisible = false,  # Hide bottom spine
            topspinevisible = false,  # Hide top spine
            rightspinevisible = false,  # Hide right spine
            leftspinevisible = false,  # Hide left spine
        )
        push!(state.channel_axs, ax_channel)

        # Link axes
        linkyaxes!(ax, ax_channel)
        linkxaxes!(ax, ax_channel)

        # Channel overlay plot
        lines!(
            ax_channel,
            @lift(state.dat.data.time[$(state.xrange)]),
            @lift(
                $(state.show_channel) ? $(state.channel_data)[$(state.xrange)] .* $(state.channel_yscale) :
                zeros(Float64, length($(state.xrange)))
            ),
            color = :grey,
        )

        # Set initial x-axis limits
        xlims!(ax_channel, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))

        # Observable creation
        lines_obs = Observable(state.components[comp_idx, :])
        push!(state.lines_obs, lines_obs)

        # Component line plot
        lines!(ax, @lift(state.dat.data.time[$(state.xrange)]), @lift($(lines_obs)[$(state.xrange)]), color = :black)

        # Set initial x-axis limits
        xlims!(ax, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))

        # Create the topo plot
        if comp_idx <= state.total_components
            plot_topoplot_in_viewer!(fig, topo_ax, state.ica_result, comp_idx, state.dat.layout, use_global_scale = state.use_global_scale[])
            topo_ax.title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
        end
    end
end

# Update add_navigation_controls! to include the global scale checkbox
function add_navigation_controls!(fig, state)
    # Add navigation buttons below topo plots in column 1
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1], tellheight = false)

    # Navigation buttons in first row
    prev_topo = Button(topo_nav[1, 1], label = "◄ Previous", tellheight = false)
    next_topo = Button(topo_nav[1, 2], label = "Next ►", tellheight = false)

    # Component selection in second row
    text_label = Label(topo_nav[2, 1], "Components:", tellheight = false, halign = :right) # Align label right
    text_input = Textbox(topo_nav[2, 2], placeholder = "e.g. 1,3-5,8", tellheight = false)
    apply_button = Button(topo_nav[2, 3], label = "Apply", tellheight = false)

    # Global scale checkbox in third row
    global_scale_check = Checkbox(topo_nav[3, 1], checked = state.use_global_scale[], tellheight = false)
    Label(topo_nav[3, 2], "Use Global Scale", tellwidth = false, tellheight = false)

    # Invert scale checkbox in fourth row
    invert_scale_check = Checkbox(topo_nav[4, 1], checked = state.invert_scale[], tellheight = false)
    Label(topo_nav[4, 2], "Invert Scale", tellwidth = false, tellheight = false)

    # --- Gap settings ---
    # Add column gaps for better spacing (horizontal - unchanged)
    colgap!(topo_nav, 1, 10)
    colgap!(topo_nav, 2, 5)

    # Add row gaps for vertical spacing (slightly more increased)
    rowgap!(topo_nav, 1, 35) # Increased gap after navigation buttons
    rowgap!(topo_nav, 2, 45) # Increased gap after component selection
    rowgap!(topo_nav, 3, 35) # Increased gap after Global Scale checkbox
    # --- End Gap settings ---


    # Connect checkboxes to state
    on(global_scale_check.checked) do checked
        state.use_global_scale[] = checked
        update_components!(state)
    end

    on(invert_scale_check.checked) do checked
        state.invert_scale[] = checked
        update_components!(state)
    end

    # Connect navigation buttons
    on(prev_topo.clicks) do _
        new_start = max(1, state.comp_start[] - state.n_visible_components)
        state.comp_start[] = new_start
        update_components!(state)
    end

    on(next_topo.clicks) do _
        new_start = min(
            state.total_components - state.n_visible_components + 1,
            state.comp_start[] + state.n_visible_components,
        )
        state.comp_start[] = new_start
        update_components!(state)
    end

    # Connect apply button
    on(apply_button.clicks) do _
        text_value = text_input.displayed_string[]
        if !isempty(text_value)
            comps = parse_component_input(text_value, state.total_components)
            if !isempty(comps)
                println("Creating new plot with components: $comps")
                current_channel = state.channel_data[]
                show_channel = state.show_channel[]
                use_global = state.use_global_scale[]
                invert = state.invert_scale[]
                new_fig = plot_ica_component_activation(
                    state.dat,
                    state.ica_result,
                    specific_components = comps,
                    n_visible_components = length(comps),
                    window_size = state.window_size,
                )
            end
        end
    end
end

# Function to parse component input text into a list of component indices
function parse_component_input(text::String, total_components::Int)
    components = Int[]
    if isempty(text)
        return components
    end

    try
        # Split by comma
        parts = strip.(split(text, ','))
        for part in parts
            if occursin('-', part)
                # Handle ranges like "1-5"
                range_parts = strip.(split(part, '-'))
                if length(range_parts) == 2
                    start_num = parse(Int, range_parts[1])
                    end_num = parse(Int, range_parts[2])
                    if 1 <= start_num <= end_num <= total_components
                        append!(components, start_num:end_num)
                    end
                end
            else
                # Handle single numbers
                num = parse(Int, part)
                if 1 <= num <= total_components
                    push!(components, num)
                end
            end
        end
    catch e
        # Silently handle parsing errors
    end

    # Remove duplicates and sort
    unique!(sort!(components))
    return components
end

# Update specific components based on user input
function update_specific_components!(state, comp_indices)
    for i = 1:state.n_visible_components
        if i <= length(comp_indices)
            comp_idx = comp_indices[i]
            if comp_idx <= state.total_components
                # Update component data
                state.lines_obs[i][] = state.components[comp_idx, :]

                # Update y-axis label to show actual component number
                state.axs[i].ylabel = @sprintf("IC %d", comp_idx)

                # Clear and redraw topography
                empty!(state.topo_axs[i])
                plot_ica_topoplot(
                    state.topo_axs[i].parent,
                    state.topo_axs[i],
                    state.ica_result,
                    comp_idx,
                    state.dat.layout,
                    colorbar_kwargs = Dict(:plot_colorbar => false),
                )

                # Update title
                state.topo_axs[i].title =
                    @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
            end
        end
    end
end

# Update add_navigation_sliders! to match new layout
function add_navigation_sliders!(fig, state)
    # Create new row for position slider below the navigation buttons
    slider_row = state.n_visible_components + 3  # Move slider to a new row

    # Calculate the step size for the slider (1% of the data length)
    step_size = max(1, div(length(state.dat.data.time), 100))

    # Use a more compact style
    x_slider = Slider(
        fig[slider_row, 2],
        range = 1:step_size:length(state.dat.data.time),
        startvalue = first(state.xrange[]),
        tellwidth = false,
        tellheight = false,
        width = Auto(),
    )

    on(x_slider.value) do x
        update_view_range!(state, Int(round(x)))
    end
end

# Helper function to update view range (similar to plot_databrowser style)
function update_view_range!(state, start_pos)
    # Ensure we stay within data bounds
    if start_pos + state.window_size > length(state.dat.data.time)
        start_pos = length(state.dat.data.time) - state.window_size + 1
    end
    end_pos = start_pos + state.window_size - 1

    # Update range
    state.xrange[] = start_pos:end_pos

    # Update axis limits based on the current view range
    first_idx = clamp(first(state.xrange[]), 1, length(state.dat.data.time))
    last_idx = clamp(last(state.xrange[]), 1, length(state.dat.data.time))
    new_xlims = (state.dat.data.time[first_idx], state.dat.data.time[last_idx])

    # Update all axes including channel axes
    for ax in state.axs
        xlims!(ax, new_xlims)
    end
    for ax in state.channel_axs
        if !isnothing(ax)  # Only update if the axis exists
            xlims!(ax, new_xlims)
        end
    end

    # If we have a boolean channel selected, update its indicators
    if state.show_channel[] && !isempty(state.channel_bool_indicators)
        # Clear existing indicators
        for (i, indicator) in state.channel_bool_indicators
            if !isnothing(indicator)
                delete!(state.channel_axs[i], indicator)
            end
        end
        empty!(state.channel_bool_indicators)

        # Find the current channel
        current_channel = nothing
        for col in names(state.dat.data)
            if state.dat.data[!, col] == state.channel_data[]
                current_channel = Symbol(col)
                break
            end
        end

        # If we found a boolean channel, redraw its indicators
        if !isnothing(current_channel) && eltype(state.dat.data[!, current_channel]) == Bool
            add_boolean_indicators!(state, current_channel)
        end
    end
end

# Update add_channel_menu! to match new layout
function add_channel_menu!(fig, state)
    # Create a menu layout in column 2
    menu_row = state.n_visible_components + 2  # Keep menu in its own row

    # Create a simple grid layout
    menu_layout = GridLayout(fig[menu_row, 2], tellheight = false)

    # Create a simple label and menu
    Label(menu_layout[1, 1], "Additional Channel:", fontsize = 18, tellheight = false, width = 150)  # Fixed width for label
    channel_menu = Menu(menu_layout[1, 2], options = ["None"; names(state.dat.data)], default = "None", tellheight = false)

    on(channel_menu.selection) do selected
        update_channel_selection!(state, selected)
    end

    return menu_layout
end

# Helper function to update channel selection (matches databrowser pattern)
function update_channel_selection!(state, selected)
    # Clear previous channel visualizations from all axes
    for i = 1:state.n_visible_components
        if i <= length(state.channel_axs) &&
           haskey(state.channel_bool_indicators, i) &&
           !isnothing(state.channel_bool_indicators[i])
            delete!(state.channel_axs[i], state.channel_bool_indicators[i])
            state.channel_bool_indicators[i] = nothing
        end
    end
    empty!(state.channel_bool_indicators)

    if selected == "None"
        state.show_channel[] = false
        state.channel_data[] = zeros(Float64, size(state.dat.data, 1))
    else
        selected_sym = Symbol(selected)
        state.channel_data[] = state.dat.data[!, selected_sym]

        # Only show overlay plot for non-Boolean channels
        if eltype(state.dat.data[!, selected_sym]) == Bool
            state.show_channel[] = false
            add_boolean_indicators!(state, selected_sym)
        else
            state.show_channel[] = true
        end
    end
end

# Helper function to add boolean indicators (matching databrowser pattern)
function add_boolean_indicators!(state, channel_sym)
    # For each component axis, create a vertical line at each true position
    for i = 1:state.n_visible_components
        if i <= length(state.channel_axs)
            ax_channel = state.channel_axs[i]

            # Find all time points where the boolean is true
            true_indices = findall(state.dat.data[!, channel_sym])

            if !isempty(true_indices)
                # Get the time values for the true positions
                true_times = state.dat.data.time[true_indices]

                # Create vertical lines at each true position
                # Only create lines within the current view range
                current_range = state.xrange[]
                visible_times = true_times[true_times.>=state.dat.data.time[first(
                    current_range,
                )].&&true_times.<=state.dat.data.time[last(current_range)]]

                if !isempty(visible_times)
                    lines = vlines!(ax_channel, visible_times, color = :red, linewidth = 1)

                    # Store the reference to the lines
                    state.channel_bool_indicators[i] = lines
                end
            end
        end
    end
end

# Setup keyboard interactions
function setup_keyboard_interactions!(fig, state)
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            if event.key == Keyboard.left || event.key == Keyboard.right
                # Handle x-axis scrolling
                current_range = state.xrange[]
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - state.window_size)
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                else  # right
                    new_start =
                        min(size(state.components, 2) - state.window_size + 1, first(current_range) + state.window_size)
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                end

                # Update x-axis limits for all axes
                # Ensure the indices are within bounds
                first_idx = clamp(first(state.xrange[]), 1, length(state.dat.data.time))
                last_idx = clamp(last(state.xrange[]), 1, length(state.dat.data.time))
                new_xlims = (state.dat.data.time[first_idx], state.dat.data.time[last_idx])
                for ax in state.axs
                    xlims!(ax, new_xlims)
                end
                for ax in state.channel_axs
                    if !isnothing(ax)
                        xlims!(ax, new_xlims)
                    end
                end

            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed =
                    (Keyboard.left_shift in events(fig).keyboardstate) ||
                    (Keyboard.right_shift in events(fig).keyboardstate)

                if !shift_pressed
                    # Handle y-axis scaling
                    current_range = state.ylims[][2]  # Just take the positive limit
                    if event.key == Keyboard.up
                        # Zoom in - decrease range by 20%
                        new_range = current_range * 0.8
                    else  # down
                        # Zoom out - increase range by 20%
                        new_range = current_range * 1.2
                    end

                    # Keep centered on zero
                    state.ylims[] = (-new_range, new_range)

                    # Update y-axis limits for all axes
                    for ax in state.axs
                        ylims!(ax, state.ylims[])
                    end
                    # Also update channel axes to maintain alignment
                    for ax in state.channel_axs
                        if !isnothing(ax)
                            ylims!(ax, state.ylims[])
                        end
                    end

                    # Force a redraw of the plots with scale inversion
                    for i = 1:state.n_visible_components
                        # Get the correct component index based on whether we're using specific components
                        comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
                            state.specific_components[i]
                        else
                            state.comp_start[] + i - 1
                        end

                        if comp_idx <= state.total_components
                            component_data = state.components[comp_idx, :]
                            if state.invert_scale[]
                                component_data = -component_data
                            end
                            state.lines_obs[i][] = component_data
                        end
                    end
                else
                    # With shift - adjust ONLY channel scale without changing axis limits
                    if event.key == Keyboard.up && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] * 1.1
                    elseif event.key == Keyboard.down && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] / 1.1
                    end
                    # We don't modify any axis limits here, only the scaling factor
                    # This will affect how the channel data is plotted through the Observable
                end

            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                # Only handle page up/down if we're not using specific components
                if isnothing(state.specific_components)
                    # Handle component scrolling
                    current_start = state.comp_start[]
                    if event.key == Keyboard.page_up
                        new_start = max(1, current_start - state.n_visible_components)
                    else  # page_down
                        new_start = min(
                            state.total_components - state.n_visible_components + 1,
                            current_start + state.n_visible_components,
                        )
                    end

                    if new_start != current_start
                        state.comp_start[] = new_start
                        update_components!(state)
                    end
                end
            end
        end
    end
end

# Update component data when navigating
function update_components!(state)
    # Calculate global min/max if using global scale
    global_min = nothing
    global_max = nothing
    if state.use_global_scale[]
        all_values = Float64[]
        for i = 1:state.n_visible_components
            # Get the correct component index based on whether we're using specific components
            comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
                state.specific_components[i]
            else
                state.comp_start[] + i - 1
            end

            if comp_idx <= state.total_components
                # Extract layout data
                tmp_layout = state.dat.layout[(state.dat.layout.label.∈Ref(state.ica_result.data_label)), :]
                
                # Create the topo data
                gridscale = 300
                data = data_interpolation_topo(
                    state.ica_result.mixing[:, comp_idx],
                    permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
                    gridscale,
                )
                
                # Add non-NaN values to our collection
                append!(all_values, [v for v in data if !isnan(v)])
            end
        end
        
        if !isempty(all_values)
            global_min = minimum(all_values)
            global_max = maximum(all_values)
        end
    end

    for i = 1:state.n_visible_components
        # Get the correct component index based on whether we're using specific components
        comp_idx = if !isnothing(state.specific_components) && i <= length(state.specific_components)
            state.specific_components[i]
        else
            state.comp_start[] + i - 1
        end

        if comp_idx <= state.total_components
            # Update component data with possible inversion
            component_data = state.components[comp_idx, :]
            if state.invert_scale[]
                component_data = -component_data
            end
            state.lines_obs[i][] = component_data

            # Clear and redraw topography
            empty!(state.topo_axs[i])
            plot_topoplot_in_viewer!(
                state.topo_axs[i].parent,
                state.topo_axs[i],
                state.ica_result,
                comp_idx,
                state.dat.layout,
                use_global_scale = state.use_global_scale[],
                global_min = global_min,
                global_max = global_max
            )

            # Update title
            state.topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
        end
    end
end
