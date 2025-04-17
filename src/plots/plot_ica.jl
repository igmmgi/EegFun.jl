function plot_ica_topoplot_single(
    fig,
    position,
    ica,
    comp_idx,
    layout;
    show_colorbar = false,
    colorbar_kwargs = Dict(),
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict()
)
    # Create a consistent GridLayout at this position - always the same structure
    gl = fig[position...] = GridLayout()
    
    # Create main axis in the first column
    ax = Axis(gl[1, 1], title = String(ica.ica_label[comp_idx]))
    
    # Extract layout data
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]
    
    # Create the topo data
    gridscale = 300  # Default gridscale
    data = data_interpolation_topo(
        ica.mixing[:, comp_idx],
        permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
        gridscale
    )
    
    # Create the plot
    radius = 88 # mm
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data,
        colormap = :jet
    )
    
    # Add head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs
    )
    
    # Hide decorations for cleaner look
    hidexdecorations!(ax, grid=false)
    hideydecorations!(ax, grid=false)
    
    # Add colorbar if requested
    if show_colorbar
        # Default values
        width = get(colorbar_kwargs, :width, 10)
        height = get(colorbar_kwargs, :height, Relative(0.8))
        ticksize = get(colorbar_kwargs, :ticklabelsize, 10)
        
        # Create colorbar in second column
        Colorbar(
            gl[1, 2],
            co;
            width = width,
            height = height,
            ticklabelsize = ticksize
        )
        
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
)
    # Process inputs
    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
        polar_to_cartesian_xy!(layout)
    end
    
    if isnothing(comps)
        comps = 1:size(ica.mixing)[2]
    end
    
    # Get colorbar setting
    plot_colorbar = get(colorbar_kwargs, :plot_colorbar, true)
    
    # Create figure
    fig = Figure()
    
    # Calculate layout dimensions
    if isnothing(dims)
        dims = best_rect(length(comps))
    end
    
    # Extract layout data once for all plots
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]
    gridscale = 300
    radius = 88
    
    # Calculate global min/max for consistent colors across plots
    all_data = []
    for i in 1:length(comps)
        data = data_interpolation_topo(
            ica.mixing[:, comps[i]],
            permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
            gridscale
        )
        push!(all_data, data)
    end
    
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
    levels = range(data_min, data_max, length=20)
    
    # First create all the GridLayouts for consistent sizing
    grids = []
    for i = 1:length(comps)
        # Calculate grid position
        row = ceil(Int, i/dims[2])
        col = ((i-1) % dims[2]) + 1
        
        # Create a GridLayout for each plot
        gl = fig[row, col] = GridLayout()
        push!(grids, gl)
    end
    
    # Now add the plots and colorbars
    for i = 1:length(comps)
        gl = grids[i]
        
        # Create main axis in the first column
        ax = Axis(gl[1, 1], title = String(ica.ica_label[comps[i]]))
        
        # Get pre-calculated data
        data = all_data[i]
        
        # Create the plot with consistent levels
        co = contourf!(
            ax,
            range(-radius * 2, radius * 2, length = gridscale),
            range(-radius * 2, radius * 2, length = gridscale),
            data,
            colormap = :jet,
            levels = levels,
            nan_color = :transparent
        )
        
        # Add head shape
        plot_layout_2d!(
            fig,
            ax,
            layout,
            head_kwargs = head_kwargs,
            point_kwargs = point_kwargs,
            label_kwargs = label_kwargs
        )
        
        # Hide decorations for cleaner look
        hidexdecorations!(ax, grid=false)
        hideydecorations!(ax, grid=false)
        
        # Set up colorbar space for ALL plots for consistent sizing
        # Default values
        width = get(colorbar_kwargs, :width, 10)
        height = get(colorbar_kwargs, :height, Relative(0.8))
        ticksize = get(colorbar_kwargs, :ticklabelsize, 10)
        
        # Set column sizes for all - ensures consistent layout
        colsize!(gl, 1, Relative(0.85))  # Plot column
        
        if plot_colorbar
            if i == 1 || i == 3
                # Only create visible colorbar for the first plot
                Colorbar(
                    gl[1, 2],
                    co;
                    width = width,
                    height = height,
                    ticklabelsize = ticksize
                )
            else
                # For other plots, create an empty colorbar placeholder
                # Create an empty axis with the same dimensions but no content
                placeholder = Axis(gl[1, 2];
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
                    ygridvisible = false
                )
                # Hide everything to make it truly empty
                hidexdecorations!(placeholder)
                hideydecorations!(placeholder)
            end
            colsize!(gl, 2, Relative(0.15))  # Colorbar column
        else
            # If no colorbar wanted at all, just use the full width
            colsize!(gl, 1, Relative(1.0))
        end
    end
    
    # Apply consistent sizing to all grid rows and columns
    # This ensures alignment across the entire figure
    for i = 1:dims[1]
        rowsize!(fig.layout, i, Relative(1/dims[1]))
    end
    
    for j = 1:dims[2]
        colsize!(fig.layout, j, Relative(1/dims[2]))
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
    colorbar_kwargs = Dict(),
    head_kwargs = Dict(),
    point_kwargs = Dict(),
    label_kwargs = Dict()
)
    # Clear the axis before drawing
    empty!(ax)
    
    # Extract layout data
    tmp_layout = layout[(layout.label.∈Ref(ica.data_label)), :]
    
    # Create the topo data
    gridscale = 300  # Default gridscale
    data = data_interpolation_topo(
        ica.mixing[:, comp_idx],
        permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
        gridscale
    )
    
    # Create the plot directly on the existing axis
    radius = 88 # mm
    co = contourf!(
        ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data,
        colormap = :jet
    )
    
    # Add head shape
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs
    )
    
    # Set title - for component viewer
    ax.title = string(ica.ica_label[comp_idx])
    
    # Hide decorations for cleaner look
    hidexdecorations!(ax, grid=false)
    hideydecorations!(ax, grid=false)
    
    return co
end



# function plot_ica_component_activation(ica, epochs::EpochData, component::Int)
#     fig = Figure()
# 
#     # Create grid layout and make it fill the figure
#     gl = fig[1, 1] = GridLayout()
#     colsize!(gl, 1, Relative(1.0))
# 
#     # Get data for plotting
#     activation = ica.activation[component, :]
#     projection = ica.mixing[:, component] .* activation'
# 
#     # Calculate y-limits for consistent scaling
#     y_min, y_max = extrema(activation)
#     y_range = y_max - y_min
#     limits = (y_min - 0.1 * y_range, y_max + 0.1 * y_range)
# 
#     # Plot component activation
#     ax1 = Axis(gl[1, 1], title = "Component $component Activation", limits = (nothing, limits))
#     lines!(ax1, epochs.time, activation)
# 
#     # Plot channel projections
#     ax2 = Axis(gl[2, 1], title = "Channel Projections")
#     for (i, chan) in enumerate(epochs.layout.label)
#         lines!(ax2, epochs.time, projection[i, :], color = :lightgrey, linewidth = 1, alpha = 0.5)
#     end
# 
#     # Set row sizes to give equal space to plots
#     rowsize!(gl, 1, Relative(0.5))
#     rowsize!(gl, 2, Relative(0.5))
# 
#     rowgap!(gl, 10)  # Add gap between plots
# 
#     # Link x-axes
#     linkxaxes!(ax1, ax2)
# 
#     return fig
# end

# epoch = extract_epochs(dat, 1, 1, -2, 4)
# plot_ica_component_activation(ica_result, epoch, 1)






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
    
    # New field for specific components
    specific_components::Union{Nothing, Vector{Int}}
    
    # Plot elements
    axs::Vector{Axis}
    channel_axs::Vector{Union{Axis, Nothing}}  # Allow for nothing values
    topo_axs::Vector{Axis}
    lines_obs::Vector{Observable{Vector{Float64}}}
    
    # New field for boolean indicators
    channel_bool_indicators::Dict{Int, Any}
    
    function IcaComponentState(dat, ica_result, n_visible_components, window_size, specific_components=nothing)
        # Prepare data matrix
        dat_matrix = prepare_ica_data_matrix(dat, ica_result)
        components = ica_result.unmixing * dat_matrix
        total_components = size(components, 1)
        
        # Create observables
        comp_start = Observable(1)
        
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
        channel_axs = Vector{Union{Axis, Nothing}}()  # Allow for nothing values
        topo_axs = Vector{Axis}()
        lines_obs = Vector{Observable{Vector{Float64}}}()
        channel_bool_indicators = Dict{Int, Any}()
        
        new(
            dat, ica_result, components, total_components,
            n_visible_components, window_size,
            comp_start, xrange, ylims, channel_data, show_channel, channel_yscale,
            specific_components,
            axs, channel_axs, topo_axs, lines_obs, channel_bool_indicators
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
    specific_components = nothing
)
    # Create state, using specific components if provided
    if !isnothing(specific_components)
        # Make sure n_visible_components matches the length of specific_components
        n_visible_components = length(specific_components)
    end
    
    state = IcaComponentState(dat, ica_result, n_visible_components, window_size, specific_components)
    
    # Create figure
    fig = Figure()
    
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
function plot_topoplot_in_viewer!(
    fig,
    topo_ax, 
    ica_result, 
    comp_idx, 
    layout
)
    # Clear the axis
    empty!(topo_ax)
    
    # Extract layout data and prepare the plot
    tmp_layout = layout[(layout.label.∈Ref(ica_result.data_label)), :]
    
    # Create the topo data
    gridscale = 300
    data = data_interpolation_topo(
        ica_result.mixing[:, comp_idx],
        permutedims(Matrix(tmp_layout[!, [:x2, :y2]])),
        gridscale
    )
    
    # Plot directly on the provided axis
    radius = 88 # mm
    co = contourf!(
        topo_ax,
        range(-radius * 2, radius * 2, length = gridscale),
        range(-radius * 2, radius * 2, length = gridscale),
        data,
        colormap = :jet
    )
    
    # Add head shape with default settings
    plot_layout_2d!(fig, topo_ax, layout)
    
    # Hide decorations
    hidexdecorations!(topo_ax, grid=false)
    hideydecorations!(topo_ax, grid=false)
    
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
            yticksvisible = true  # Keep the tick marks themselves
        )
        push!(state.axs, ax)
        
        # Channel overlay axis (only created when needed)
        ax_channel = nothing
        if state.show_channel[] && !isnothing(state.channel_data[]) && 
           !isempty(state.channel_data[]) && eltype(state.channel_data[]) != Bool
            ax_channel = Axis(
                fig[i, 2],
                yticklabelsvisible = false,
                yticksvisible = false,
                yaxisposition = :right,
                xaxisposition = :top  # Add this to ensure proper time axis alignment
            )
            push!(state.channel_axs, ax_channel)
            
            # Link axes
            linkyaxes!(ax, ax_channel)
            linkxaxes!(ax, ax_channel)  # Add this to ensure time axis stays synchronized
            
            # Channel overlay plot
            lines!(
                ax_channel,
                @lift(state.dat.data.time[$(state.xrange)]),
                @lift($(state.show_channel) ? $(state.channel_data)[$(state.xrange)] .* $(state.channel_yscale) : zeros(length($(state.xrange)))),
                color = :red
            )
            
            # Set initial x-axis limits
            xlims!(ax_channel, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))
        else
            # Push a placeholder nothing for consistency
            push!(state.channel_axs, nothing)
        end
        
        # Observable creation
        lines_obs = Observable(state.components[comp_idx, :])
        push!(state.lines_obs, lines_obs)
        
        # Component line plot
        lines!(
            ax,
            @lift(state.dat.data.time[$(state.xrange)]),
            @lift($(lines_obs)[$(state.xrange)])
        )
        
        # Set initial x-axis limits
        xlims!(ax, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))
        
        # Hide x-axis decorations for all plots except bottom axis of last plot
        if i != state.n_visible_components
            hidexdecorations!(ax, grid=false)
        end
        
        # Create the topo plot using our new simpler function
        if comp_idx <= state.total_components
            # Use the simpler function that avoids gridposition issues
            plot_topoplot_in_viewer!(fig, topo_ax, state.ica_result, comp_idx, state.dat.layout)
            
            # Update title
            topo_ax.title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
        end
    end
    
    # Set column sizes to give more space to the time series plots
    colsize!(fig.layout, 1, Relative(0.15))  # Topoplots - now narrower
    colsize!(fig.layout, 2, Relative(0.85))  # Time series - now wider
end

# Update add_navigation_controls! to match new layout
function add_navigation_controls!(fig, state)
    # Add navigation buttons below topo plots in column 1
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1])  # Only in column 1
    prev_topo = Button(topo_nav[1, 1], label = "◄ Previous")
    next_topo = Button(topo_nav[1, 2], label = "Next ►")
    
    # Add component selection text box
    text_label = Label(topo_nav[2, 1], "Components:")
    text_input = Textbox(topo_nav[2, 2], placeholder = "e.g. 1,3-5,8")
    
    # Connect navigation buttons
    on(prev_topo.clicks) do _
        new_start = max(1, state.comp_start[] - state.n_visible_components)
        state.comp_start[] = new_start
        update_components!(state)
    end
    
    on(next_topo.clicks) do _
        new_start = min(state.total_components - state.n_visible_components + 1, 
                        state.comp_start[] + state.n_visible_components)
        state.comp_start[] = new_start
        update_components!(state)
    end
    
    # Add a submit button next to the textbox to explicitly apply changes
    apply_button = Button(topo_nav[2, 3], label = "Apply")
    on(apply_button.clicks) do _
        # Get the raw string value directly from the textbox
        text_value = text_input.displayed_string[]
        
        # Only process if we have a valid text input
        if !isempty(text_value)
            # Convert to component indices
            comps = parse_component_input(text_value, state.total_components)
            
            if !isempty(comps)
                println("Creating new plot with components: $comps")
                
                # Close current figure and create a new one with just these components
                # Get current settings to preserve them
                current_channel = state.channel_data[]
                show_channel = state.show_channel[]
                
                # Create a new figure with exactly these components
                new_fig = plot_ica_component_activation(
                    state.dat, 
                    state.ica_result,
                    specific_components=comps, 
                    n_visible_components=length(comps),
                    window_size=state.window_size
                )
                
                # The current figure will be replaced by the new one in the display
            else
                println("No valid components found in input: $text_value")
            end
        else
            println("Empty text input")
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
                state.topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
            end
        end
    end
end

# Update add_navigation_sliders! to match new layout
function add_navigation_sliders!(fig, state)
    # Create new row for position slider below the navigation buttons
    slider_row = state.n_visible_components + 3
    
    # Calculate the step size for the slider (1% of the data length)
    step_size = max(1, div(length(state.dat.data.time), 100))
    
    # Use a more consistent style matching plot_databrowser
    x_slider = Slider(
        fig[slider_row, 2],  # Now in column 2 where the time series are
        range = 1:step_size:length(state.dat.data.time),
        startvalue = first(state.xrange[]),
        tellwidth = false,
        width = Auto()
    )
    
    # Connect slider to state using the same pattern as in plot_databrowser
    on(x_slider.value) do x
        # Update view range
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
    menu_row = state.n_visible_components + 2
    
    # Create a simple grid layout
    menu_layout = GridLayout(fig[menu_row, 2])  # Now in column 2 where the time series are
    
    # Create a simple label and menu
    Label(menu_layout[1, 1], "Additional Channel:", fontsize = 18)
    
    # Use a standard menu without fixed width
    channel_menu = Menu(
        menu_layout[1, 2],
        options = ["None"; names(state.dat.data)],
        default = "None",
    )
    
    # Connect menu selection callback
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
        state.channel_data[] = zeros(size(state.dat.data, 1))
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
                visible_times = true_times[true_times .>= state.dat.data.time[first(current_range)] .&& 
                                         true_times .<= state.dat.data.time[last(current_range)]]
                
                if !isempty(visible_times)
                    lines = vlines!(
                        ax_channel,
                        visible_times,
                        color = :red,
                        linewidth = 1
                    )
                    
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
                    new_start = min(size(state.components, 2) - state.window_size + 1, 
                                   first(current_range) + state.window_size)
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
                
            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed = (Keyboard.left_shift in events(fig).keyboardstate) || 
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
                        ylims!(ax, state.ylims[])
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
                # Handle component scrolling
                current_start = state.comp_start[]
                if event.key == Keyboard.page_up
                    new_start = max(1, current_start - state.n_visible_components)
                else  # page_down
                    new_start = min(state.total_components - state.n_visible_components + 1, 
                                   current_start + state.n_visible_components)
                end
                
                if new_start != current_start
                    state.comp_start[] = new_start
                    update_components!(state)
                end
            end
        end
    end
end

# Update component data when navigating
function update_components!(state)
    for i = 1:state.n_visible_components
        comp_idx = state.comp_start[] + i - 1
        if comp_idx <= state.total_components
            # Update component data
            state.lines_obs[i][] = state.components[comp_idx, :]
            
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
            state.topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
        end
    end
end
