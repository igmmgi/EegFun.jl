

"""
    plot_ica_component_activation(dat::ContinuousData, ica_result::InfoIca; kwargs...)

Creates an interactive visualization of ICA components with topographic maps and time series.

# Arguments
- `dat::ContinuousData`: Continuous EEG data
- `ica_result::InfoIca`: ICA result object

# Keyword Arguments
- `n_visible_components::Int=10`: Number of components visible at once
- `window_size::Int=2000`: Initial window size in samples
- `topo_kwargs::Dict=Dict()`: Additional keyword arguments for topoplots

# Returns
- `fig::Figure`: The Figure object containing the interactive plot

# Example
```julia
fig = plot_ica_component_activation(dat, ica_result)
```
"""
function plot_ica_component_activation(
    dat::ContinuousData,
    ica_result::InfoIca;
    n_visible_components::Int = 10,
    window_size::Int = 2000,
    topo_kwargs = Dict()
)
    # Create state object to hold all visualization state
    state = create_ica_browser_state(dat, ica_result, n_visible_components, window_size)
    
    # Create figure and main layout
    fig = Figure(size = (1200, 800), fontsize = 18)
    
    # Setup subplots and axes
    axs, topo_axs = setup_ica_layout(fig, state)
    
    # Setup UI controls
    setup_ica_controls(fig, axs, topo_axs, state)
    
    # Initial draw
    update_component_plots(topo_axs, axs, state)
    
    display(fig)
    return fig
end

# Helper state structure to hold all visualization state
mutable struct IcaBrowserState
    # Data
    dat::ContinuousData
    ica_result::InfoIca
    components::Matrix{Float64}
    
    # View state
    n_visible_components::Int
    total_components::Int
    window_size::Int
    comp_start::Observable{Int}
    xrange::Observable{UnitRange{Int}}
    ylims::Observable{Tuple{Float64,Float64}}
    
    # Channel visualization
    channel_data::Observable{Vector{Float64}}
    show_channel::Observable{Bool}
    channel_yscale::Observable{Float64}
    
    # Constructor
    function IcaBrowserState(dat, ica_result, n_visible_components, window_size)
        # Prepare data
        dat_matrix = prepare_ica_data(dat, ica_result)
        components = ica_result.unmixing * dat_matrix
        
        # Initialize observables
        comp_start = Observable(1)
        xrange = Observable(1:window_size)
        initial_range = maximum(abs.(extrema(components[1:n_visible_components, 1:window_size])))
        ylims = Observable((-initial_range, initial_range))
        
        # Channel state
        channel_data = Observable(zeros(size(dat.data, 1)))
        show_channel = Observable(false)
        channel_yscale = Observable(1.0)
        
        new(
            dat, 
            ica_result, 
            components,
            n_visible_components,
            size(components, 1),
            window_size,
            comp_start,
            xrange,
            ylims,
            channel_data,
            show_channel,
            channel_yscale
        )
    end
end

function create_ica_browser_state(dat, ica_result, n_visible_components, window_size)
    return IcaBrowserState(dat, ica_result, n_visible_components, window_size)
end

function prepare_ica_data(dat::ContinuousData, ica_result::InfoIca)
    dat_matrix = permutedims(Matrix(dat.data[!, ica_result.data_label]))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica_result.scale
    return dat_matrix
end

function setup_ica_layout(fig, state)
    # Create arrays to store axes
    axs = Vector{Axis}(undef, state.n_visible_components)
    topo_axs = Vector{Axis}(undef, state.n_visible_components)
    
    # Create subplots
    for i = 1:state.n_visible_components
        # Topo plot
        topo_axs[i] = Axis(
            fig[i, 1], 
            width = Relative(1),
            height = Relative(1),
            title = @sprintf("IC %d (%.1f%%)", i, state.ica_result.variance[i] * 100)
        )
        hidexdecorations!(topo_axs[i])
        hideydecorations!(topo_axs[i])
        
        # Time series plot
        axs[i] = Axis(fig[i, 2], ylabel = "ICA Amplitude")
        
        # Only add x-axis label and ticks to the bottom plot
        if i == state.n_visible_components
            axs[i].xlabel = "Time"
            axs[i].xticklabelsvisible = true
            axs[i].yticklabelsvisible = true
        else
            axs[i].xticklabelsvisible = false
            hidexdecorations!(axs[i], grid = false)
        end
    end
    
    # Link all time series axes for consistent zooming/panning
    linkaxes!(axs...)
    
    # Set column sizes
    colsize!(fig.layout, 1, Auto(150))  # Fixed width for topo plots
    
    return axs, topo_axs
end

function setup_ica_controls(fig, axs, topo_axs, state)
    # Add navigation controls
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1])
    prev_button = Button(topo_nav[1, 1], label = "◄ Previous")
    next_button = Button(topo_nav[1, 2], label = "Next ►")
    
    # Channel selection controls
    menu_row = GridLayout(fig[state.n_visible_components+1, 2])
    Label(menu_row[1, 1], "Additional Channel", tellwidth = true)
    available_channels = names(state.dat.data)
    channel_menu = Menu(menu_row[1, 2], options = ["None"; available_channels], default = "None")
    
    # Setup button callbacks
    on(prev_button.clicks) do _
        new_start = max(1, state.comp_start[] - state.n_visible_components)
        state.comp_start[] = new_start
        update_component_plots(topo_axs, axs, state)
    end
    
    on(next_button.clicks) do _
        new_start = min(
            state.total_components - state.n_visible_components + 1, 
            state.comp_start[] + state.n_visible_components
        )
        state.comp_start[] = new_start
        update_component_plots(topo_axs, axs, state)
    end
    
    # Channel menu callback
    on(channel_menu.selection) do selected
        if selected == "None"
            state.show_channel[] = false
            state.channel_data[] = zeros(size(state.dat.data, 1))
        else
            state.show_channel[] = true
            state.channel_data[] = state.dat.data[!, selected]
        end
    end
    
    # Keyboard interactions
    setup_keyboard_interactions(fig, axs, topo_axs, state)
end

function setup_keyboard_interactions(fig, axs, topo_axs, state)
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press,)
            if event.key == Keyboard.left || event.key == Keyboard.right
                # Handle x-axis scrolling
                current_range = state.xrange[]
                if event.key == Keyboard.left
                    new_start = max(1, first(current_range) - state.window_size)
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                else  # right
                    new_start = min(
                        size(state.components, 2) - state.window_size + 1, 
                        first(current_range) + state.window_size
                    )
                    state.xrange[] = new_start:(new_start+state.window_size-1)
                end
                
                # Update x-axis limits for all axes
                new_xlims = (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])])
                for ax in axs
                    xlims!(ax, new_xlims)
                end
                
            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed = (Keyboard.left_shift in events(fig).keyboardstate) || 
                                (Keyboard.right_shift in events(fig).keyboardstate)
                
                if !shift_pressed
                    # Handle y-axis scaling
                    current_range = state.ylims[][2]  # Take positive limit (symmetric)
                    if event.key == Keyboard.up
                        # Zoom in - decrease range by 20%
                        new_range = current_range * 0.8
                    else  # down
                        # Zoom out - increase range by 20%
                        new_range = current_range * 1.2
                    end
                    
                    # Keep centered on zero
                    state.ylims[] = (-new_range, new_range)
                    
                    # Update y-axis limits
                    for ax in axs
                        ylims!(ax, state.ylims[])
                    end
                else
                    # Handle channel scaling with shift key
                    if event.key == Keyboard.up
                        state.channel_yscale[] = state.channel_yscale[] * 1.1
                    elseif event.key == Keyboard.down
                        state.channel_yscale[] = state.channel_yscale[] / 1.1
                    end
                end
            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                # Handle component scrolling
                current_start = state.comp_start[]
                if event.key == Keyboard.page_up
                    new_start = max(1, current_start - state.n_visible_components)
                else  # page_down
                    new_start = min(
                        state.total_components - state.n_visible_components + 1, 
                        current_start + state.n_visible_components
                    )
                end
                
                if new_start != current_start
                    state.comp_start[] = new_start
                    update_component_plots(topo_axs, axs, state)
                end
            end
        end
    end
end

function update_component_plots(topo_axs, axs, state)
    for i = 1:state.n_visible_components
        comp_idx = state.comp_start[] + i - 1
        if comp_idx <= state.total_components
            # Get component data for this position
            comp_data = state.components[comp_idx, :]
            
            # Clear existing content
            empty!(topo_axs[i])
            empty!(axs[i])
            
            # Draw topography
            plot_ica_topoplot(
                topo_axs[i].parent, 
                topo_axs[i], 
                state.ica_result, 
                comp_idx, 
                state.dat.layout, 
                colorbar_kwargs = Dict(:plot_colorbar => false)
            )
            
            # Draw time series
            lines!(
                axs[i], 
                state.dat.data.time, 
                comp_data, 
                color = :black
            )
            
            # Draw channel data if visible
            if state.show_channel[]
                lines!(
                    axs[i],
                    state.dat.data.time,
                    state.channel_data[] * state.channel_yscale[],
                    color = :grey
                )
            end
            
            # Set limits
            xlims!(axs[i], state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])])
            ylims!(axs[i], state.ylims[])
            
            # Update title
            topo_axs[i].title = @sprintf("IC %d (%.1f%%)", comp_idx, state.ica_result.variance[comp_idx] * 100)
        end
    end
end

p