# TODO: butterfly plot/global field power
# TODO: spline interpolation for topoplots?

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
    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
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
    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300, :size => 1)
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
            ax = Axis(
                fig[dim1, dim2],
                width = Relative(topo_kwargs[:size]),
                height = Relative(topo_kwargs[:size]),
                halign = 0.5,
                valign = 0.5,
            )
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
        #  end
        # head shape
        plot_layout_2d!(
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

    if (:x2 ∉ propertynames(layout) || :y2 ∉ propertynames(layout))
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
    plot_layout_2d!(
        fig,
        ax,
        layout,
        head_kwargs = head_kwargs,
        point_kwargs = point_kwargs,
        label_kwargs = label_kwargs,
    )
    # end
    return fig
end



# layout = read_layout("./layouts/biosemi72.csv");
# dat = read_bdf("../Flank_C_3.bdf");
# dat = create_eeg_dataframe(dat, layout);
# filter_data!(dat, "hp", "iir", 1, order = 1)
# rereference!(dat, :avg)
# diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG);
# diff_channel!(dat, :F9, :F10, :hEOG);
# # autodetect EOG signals
# detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
# detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
# is_extreme_value!(dat, dat.layout.label, 50);
# dat_ica = filter_data(dat, "hp", "iir", 1, order = 1)
# good_samples = findall(dat_ica.data[!, :is_extreme_value] .== false)
# good_channels = setdiff(dat_ica.layout.label, [:PO9])
# dat_for_ica = create_ica_data_matrix(dat_ica.data, good_channels, samples_to_include = good_samples)
# ica_result = infomax_ica(dat_for_ica, good_channels, n_components = length(good_channels) - 1, params = IcaPrms())

# plot_ica_topoplot(ica_result, dat.layout)
# plot_ica_topoplot(ica_result, dat.layout, comps = 1:10)
# plot_ica_topoplot(ica_result, dat.layout, comps = 1:3)
# plot_ica_topoplot(ica_result, dat.layout, comps = 1)
# plot_ica_topoplot(ica_result, dat.layout, comps = 1:15)
# plot_ica_topoplot(ica_result, dat.layout, comps = [1,3, 10])

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
    
    # Plot elements
    axs::Vector{Axis}
    channel_axs::Vector{Axis}  # Store channel axes separately
    topo_axs::Vector{Axis}
    lines_obs::Vector{Observable{Vector{Float64}}}
    
    function IcaComponentState(dat, ica_result, n_visible_components, window_size)
        # Prepare data matrix
        dat_matrix = prepare_ica_data_matrix(dat, ica_result)
        components = ica_result.unmixing * dat_matrix
        total_components = size(components, 1)
        
        # Create observables
        comp_start = Observable(1)
        # Find index closest to time 0 to center the initial view
        time_zero_idx = findmin(abs.(dat.data.time))[2]
        half_window = div(window_size, 2)
        start_idx = 1 #max(1, time_zero_idx - half_window)
        end_idx = min(size(components, 2), start_idx + window_size - 1)
        # Adjust start_idx if end_idx reached the boundary
        if end_idx == size(components, 2)
            start_idx = max(1, end_idx - window_size + 1)
        end
        println(start_idx, end_idx)
        xrange = Observable(start_idx:end_idx)
        
        initial_range = maximum(abs.(extrema(components[1:n_visible_components, 1:window_size])))
        ylims = Observable((-initial_range, initial_range))
        channel_data = Observable(zeros(size(dat.data, 1)))
        show_channel = Observable(false)
        channel_yscale = Observable(1.0)
        
        # Initialize empty plot element arrays
        axs = Vector{Axis}()
        channel_axs = Vector{Axis}()
        topo_axs = Vector{Axis}()
        lines_obs = Vector{Observable{Vector{Float64}}}()
        
        new(
            dat, ica_result, components, total_components,
            n_visible_components, window_size,
            comp_start, xrange, ylims, channel_data, show_channel, channel_yscale,
            axs, channel_axs, topo_axs, lines_obs
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
    topo_kwargs = Dict()
)
    # Create state
    state = IcaComponentState(dat, ica_result, n_visible_components, window_size)
    
    # Create figure
    fig = Figure()
    
    # Setup plots
    create_component_plots!(fig, state, topo_kwargs)
    
    # Add controls
    add_navigation_controls!(fig, state)
    
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

# Create the component plots (topoplots and time series)
function create_component_plots!(fig, state, topo_kwargs)
    # Set defaults for topo plots
    topo_default_kwargs = Dict(:colormap => :jet, :gridscale => 300, :size => 1)
    topo_kwargs = merge(topo_default_kwargs, topo_kwargs)
    
    # Runtime calculated values
    xlims = @lift((state.dat.data.time[first($(state.xrange))], state.dat.data.time[last($(state.xrange))]))
    
    # Force immediate calculation of xlims to ensure proper initial view
    initial_xlims = (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])])
    
    # Create plots for each component
    for i = 1:state.n_visible_components
        # Create topoplot
        ax_topo = Axis(
            fig[i, 1],
            width = Relative(1),
            height = Relative(1),
            title = @sprintf("IC %d (%.1f%%)", i, state.ica_result.variance[i] * 100)
        )
        push!(state.topo_axs, ax_topo)
        
        # Add topoplot
        plot_ica_topoplot(fig, ax_topo, state.ica_result, i, state.dat.layout, 
                           colorbar_kwargs = Dict(:plot_colorbar => false))
        hidexdecorations!(ax_topo)
        hideydecorations!(ax_topo)

        # Create time series plot with two y-axes
        ax_time = Axis(fig[i, 2], ylabel = "ICA Amplitude")
        ax_channel = Axis(
            fig[i, 2],
            ylabel = "Channel Amplitude",
            yaxisposition = :right,
            yticklabelcolor = :grey,
            ytickcolor = :grey,
            xgridvisible = false,
            ygridvisible = false,
        )
        
        # Link x-axes
        linkxaxes!(ax_time, ax_channel)
        
        # Hide the right axis' spine
        hidespines!(ax_channel, :l, :t, :b)
        
        # Only add x-axis title and tick labels to bottom plot
        if i == state.n_visible_components
            ax_time.xlabel = "Time"
            ax_time.xticklabelsvisible = true
            ax_time.yticklabelsvisible = true
        else
            ax_time.xticklabelsvisible = false
            hidexdecorations!(ax_channel, grid = false)
        end   
    # Force immediate calculation of xlims to ensure proper initial view
    initial_xlims = (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])])

        xlims!(ax_channel, initial_xlims)
        
        push!(state.axs, ax_time)
        push!(state.channel_axs, ax_channel)
        
        # Create observable for component data and add plot
        line_obs = Observable(state.components[i, :])
        lines!(ax_time, 
               @lift(state.dat.data.time[$(state.xrange)]), 
               @lift($(line_obs)[$(state.xrange)]), 
               color = :black, 
               linewidth = 2)
        
        # Add channel plot with visibility control
        lines!(
            ax_channel,
            @lift(state.dat.data.time[$(state.xrange)]),
            @lift($(state.show_channel) ? 
                 $(state.channel_data)[$(state.xrange)] * $(state.channel_yscale) : 
                 fill(NaN, length($(state.xrange)))),
            color = :grey,
        )
        
        push!(state.lines_obs, line_obs)
        
        # Set initial limits
        xlims!(ax_time, initial_xlims)
        ylims!(ax_time, state.ylims[])
        ylims!(ax_channel, state.ylims[])
    end
    
    # Link all time series axes
    linkaxes!(state.axs...)
    
    # Adjust layout
    colsize!(fig.layout, 1, Auto(150))  # Fixed width for topo plots
end

# Add navigation buttons
function add_navigation_controls!(fig, state)
    # Add navigation buttons below topo plots
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1])
    prev_topo = Button(topo_nav[1, 1], label = "◄ Previous")
    next_topo = Button(topo_nav[1, 2], label = "Next ►")
    
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
end

# Add channel selection menu
function add_channel_menu!(fig, state)
    # Create a new row in the figure layout for the menu
    menu_row = fig[state.n_visible_components+2, 1]  # Add new row at the bottom
    menu_layout = GridLayout(menu_row)  # Create layout for the menu
    Label(menu_layout[1, 1], "Additional Channel", tellwidth = true)  # Label
    available_channels = names(state.dat.data)
    channel_menu = Menu(menu_layout[2, 1], options = ["None"; available_channels], default = "None")
    
    # Menu selection callback
    on(channel_menu.selection) do selected
        if selected == "None"
            state.show_channel[] = false
            state.channel_data[] = zeros(size(state.dat.data, 1))
        else
            state.show_channel[] = true
            state.channel_data[] = state.dat.data[!, Symbol(selected)]
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
                else
                    # With shift - adjust channel scale
                    if event.key == Keyboard.up && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] * 1.1
                    elseif event.key == Keyboard.down && shift_pressed
                        state.channel_yscale[] = state.channel_yscale[] / 1.1
                    end

                    # Update channel axis limits to match scaling change
                    for i = 1:state.n_visible_components
                        # Get channel axis from our stored array
                        if i <= length(state.channel_axs)
                            # Use a different range for the channel axis
                            channel_range = state.ylims[][2] * state.channel_yscale[]
                            ylims!(state.channel_axs[i], (-channel_range, channel_range))
                        end
                    end
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

plot_ica_component_activation(dat, ica_result, n_visible_components = 10)