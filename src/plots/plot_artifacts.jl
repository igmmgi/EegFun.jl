# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

"""
    _create_artifact_controls(fig, n_epochs, epochs_with_artifacts)

Create the control buttons and labels for artifact navigation.
"""
function _create_artifact_controls(fig, n_epochs, epochs_with_artifacts)
    controls_layout = GridLayout(fig, tellheight = false, tellwidth = false)
    
    # Epoch navigation
    back_button = Button(controls_layout[1, 1], label = "◀ Previous")
    epoch_label = Label(controls_layout[1, 2], "Epoch 1 / $n_epochs", justification = :center)
    forward_button = Button(controls_layout[1, 3], label = "Next ▶")
    
    # Artifact navigation
    back_artifact_button = Button(controls_layout[1, 4], label = "◀ Previous Artifact")
    artifact_label = Label(
        controls_layout[1, 5],
        "Artifact 1 / $(length(epochs_with_artifacts))",
        justification = :center,
    )
    forward_artifact_button = Button(controls_layout[1, 6], label = "Next Artifact ▶")
    
    return (; back_button, epoch_label, forward_button, back_artifact_button, artifact_label, forward_artifact_button) # named tuple
end

"""
    _create_rejected_color_map(rejected_channels, colormap_name)

Create a color mapping for rejected channels.
"""
function _create_rejected_color_map(rejected_channels, colormap_name)
    rejected_color_map = Dict{Symbol,Any}()
    colormap = cgrad(colormap_name)
    max_colors = length(colormap.colors)
    for (i, ch) in enumerate(rejected_channels)
        color_idx = ((i - 1) % max_colors) + 1
        rejected_color_map[ch] = colormap[color_idx]
    end
    return rejected_color_map
end

"""
    _get_channel_styling(ch, is_rejected, is_selected, rejected_color_map, plot_kwargs; is_repaired=false)

Get the styling parameters for a channel based on its state.
"""
function _get_channel_styling(ch, is_rejected, is_selected, rejected_color_map, plot_kwargs; is_repaired=false)
    color = is_rejected ? rejected_color_map[ch] : :black
    alpha = is_rejected ? plot_kwargs[:alpha_rejected] : plot_kwargs[:alpha_normal]
    base_linewidth = is_rejected ? plot_kwargs[:linewidth_rejected] : plot_kwargs[:linewidth_normal]
    linewidth = is_selected ? base_linewidth * 2 : base_linewidth
    
    if is_rejected
        label_text = is_repaired ? "$ch (repaired)" : "$ch (rejected)"
    else
        label_text = "$ch"
    end
    label_text = is_selected ? rich(label_text, font = :bold) : label_text
    
    return (; color, alpha, linewidth, label_text) # named tuple
end

"""
    _plot_channels!(ax, epoch, selected_channels, rejected_channels, selected_channels_set, 
                   rejected_color_map, plot_kwargs; is_repaired=false)

Plot all channels for a single epoch on the given axis.
"""
function _plot_channels!(ax, epoch, selected_channels, rejected_channels, selected_channels_set, rejected_color_map, plot_kwargs; is_repaired=false)
    for ch in selected_channels
        is_rejected = ch in rejected_channels
        is_selected = ch in selected_channels_set
        
        styling = _get_channel_styling(ch, is_rejected, is_selected, rejected_color_map, plot_kwargs; is_repaired=is_repaired)
        
        lines!(
            ax,
            epoch.time,
            epoch[!, ch];
            color = styling.color,
            alpha = styling.alpha,
            linewidth = styling.linewidth,
            label = styling.label_text,
        )
    end
end

"""
    _create_legend!(ax, selected_channels, current_legend)

Create or update the legend for the given axis.
"""
function _create_legend!(ax, selected_channels, current_legend)
    n_channels = length(selected_channels)
    n_cols = n_channels > 10 ? cld(n_channels, 20) : 1
    current_legend[] = axislegend(ax, position = (:right, :top), nbanks = n_cols)
end

"""
    _setup_axis_grid!(ax, plot_kwargs)

Apply grid settings to the axis.
"""
function _setup_axis_grid!(ax, plot_kwargs)
    ax.xgridvisible = plot_kwargs[:xgrid]
    ax.ygridvisible = plot_kwargs[:ygrid]
    ax.xminorgridvisible = plot_kwargs[:xminorgrid]
    ax.yminorgridvisible = plot_kwargs[:yminorgrid]
end

"""
    _setup_axis_limits!(ax, plot_kwargs)

Apply axis limits to the axis.
"""
function _setup_axis_limits!(ax, plot_kwargs)
    if plot_kwargs[:xlimits] !== nothing
        xmin, xmax = plot_kwargs[:xlimits]
        xlims!(ax, xmin, xmax)
    end
    
    if plot_kwargs[:ylimits] !== nothing
        ymin, ymax = plot_kwargs[:ylimits]
        ylims!(ax, ymin, ymax)
    end
end

"""
    _setup_origin_lines!(ax, plot_kwargs)

Add origin lines at x=0 and y=0 to the axis.
"""
function _setup_origin_lines!(ax, plot_kwargs)
    if plot_kwargs[:axes_through_origin]
        hlines!(ax, 0, color = :gray, linewidth = 0.5, alpha = 0.7)
        vlines!(ax, 0, color = :gray, linewidth = 0.5, alpha = 0.7)
    end
end

"""
    _create_channel_selection_handlers(fig, ax, epochs, selected_channels, selected_channels_set, 
                                     plot_kwargs, epoch_idx, update_plot!)

Create event handlers for channel selection via shift+click.
"""
function _create_channel_selection_handlers(fig, ax, epochs, selected_channels, selected_channels_set, plot_kwargs, epoch_idx, update_plot!)

    on(events(ax).mousebutton, priority = 0) do event
        if event.button == Mouse.left && event.action == Mouse.press && _is_shift_held(fig)
            mouse_pos = _get_mouse_position(ax)
            mouse_pos === nothing && return Consume(false)

            mouse_time, mouse_amp = mouse_pos
            current_epoch = epochs.data[epoch_idx[]]
            
            closest_channel, min_distance = _find_closest_channel(mouse_time, mouse_amp, current_epoch, selected_channels)
            
            if closest_channel !== nothing && min_distance < plot_kwargs[:selection_threshold]
                _toggle_channel_selection(closest_channel, selected_channels_set)
                update_plot!()
            end
        end
        return Consume(false)
    end
end

"""
    _create_navigation_handlers(controls, epoch_idx, artifact_idx, n_epochs, epochs_with_artifacts)

Create event handlers for navigation buttons.
"""
function _create_navigation_handlers(controls, epoch_idx, artifact_idx, n_epochs, epochs_with_artifacts)

    on(controls.back_button.clicks) do _
        epoch_idx[] = max(1, epoch_idx[] - 1)
    end
    
    on(controls.forward_button.clicks) do _
        epoch_idx[] = min(n_epochs, epoch_idx[] + 1)
    end
    
    on(controls.back_artifact_button.clicks) do _
        artifact_idx[] = max(1, artifact_idx[] - 1)
    end
    
    on(controls.forward_artifact_button.clicks) do _
        artifact_idx[] = min(length(epochs_with_artifacts), artifact_idx[] + 1)
    end
end

# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ARTIFACT_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Display plot."),
    
    # line settings
    :colormap_name => (:jet, "Colormap name (e.g., :Set1_9, :tab20)."),
    :linewidth_normal => (1, "Line width for normal channels."),
    :linewidth_rejected => (2, "Line width for rejected channels."),
    :alpha_normal => (0.2, "Transparency for normal channels."),
    :alpha_rejected => (1.0, "Transparency for rejected channels."),

    # selection
    :selection_threshold => (50, "Distance threshold for channel selection."),

    # Axis limits
    :xlimits => (nothing, "X-axis limits as (min, max) tuple or nothing for auto-scaling"),
    :ylimits => (nothing, "Y-axis limits as (min, max) tuple or nothing for auto-scaling"),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    
    # Origin lines
    :axes_through_origin => (true, "Whether to add origin lines at x=0 and y=0"),
)

"""
    plot_artifact_detection(epochs::EpochData, artifacts::EpochRejectionInfo; channel_selection::Function=channels(), kwargs...)

Interactive plot of artifact detection results with Previous/Next buttons for epoch navigation.

# Arguments
- `epochs::EpochData`: The epoch data
- `artifacts::EpochRejectionInfo`: Artifact detection results
- `channel_selection::Function`: Channel predicate for selecting channels to plot (default: all layout channels)

# Keyword Arguments
$(generate_kwargs_doc(PLOT_ARTIFACT_KWARGS))

# Returns
- `Figure`: Interactive Makie figure with navigation buttons

# Examples
```julia
# Basic interactive plot
fig = plot_artifact_detection(epochs, artifacts)

# With specific channels
fig = plot_artifact_detection(epochs, artifacts, 
                            channel_selection = channels([:Fp1, :Fp2, :Cz]))

# With custom axis limits
fig = plot_artifact_detection(epochs, artifacts, xlimits = (-0.2, 0.8), ylimits = (-100, 100))
```
"""
function plot_artifact_detection(
    epochs::EpochData,
    artifacts::EpochRejectionInfo;
    channel_selection::Function = channels(),
    kwargs...,
)
    # Merge user kwargs with defaults and validate
    plot_kwargs = _merge_plot_kwargs(PLOT_ARTIFACT_KWARGS, kwargs)
    
    # Get channels to plot
    selected_channels = get_selected_channels(epochs, channel_selection, include_meta = false, include_extra = false)

    # Create figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)")

    # Setup axis styling
    _setup_axis_grid!(ax, plot_kwargs)
    _setup_axis_limits!(ax, plot_kwargs)
    _setup_origin_lines!(ax, plot_kwargs)

    # Get epochs with artifacts
    epochs_with_artifacts = unique([r.epoch for r in artifacts.rejected_epochs])
    sort!(epochs_with_artifacts)

    # Create controls
    n_epochs = length(epochs.data)
    controls = _create_artifact_controls(fig[2, 1], n_epochs, epochs_with_artifacts)

    # Set row sizes
    rowsize!(fig.layout, 1, Relative(0.9))  # Plot area
    rowsize!(fig.layout, 2, Relative(0.1))  # Controls area

    # Create observables
    epoch_idx = Observable(1)
    artifact_idx = Observable(1)
    current_legend = Ref{Union{Nothing,Legend}}(nothing)
    selected_channels_set = Set{Symbol}()

    # Create color map for rejected channels
    rejected_channels = [r.label for r in artifacts.rejected_epochs]
    rejected_color_map = _create_rejected_color_map(rejected_channels, plot_kwargs[:colormap_name])

    # Function to update plot based on epoch
    function update_plot!(epoch_idx_val)

        empty!(ax)

        # Delete existing legend if present
        if current_legend[] !== nothing
            delete!(current_legend[])
            current_legend[] = nothing
        end

        # Get current epoch data and rejected channels
        epoch = epochs.data[epoch_idx_val]
        epoch_rejected_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx_val]

        # Update title
        ax.title = "Artifact Detection - Epoch $(epoch_idx_val)"

        # Plot all channels
        _plot_channels!(ax, epoch, selected_channels, epoch_rejected_channels, selected_channels_set, 
                      rejected_color_map, plot_kwargs)

        # Create legend
        _create_legend!(ax, selected_channels, current_legend)
    end

    # Create event handlers
    _create_channel_selection_handlers(fig, ax, epochs, selected_channels, selected_channels_set, 
                                     plot_kwargs, epoch_idx, () -> update_plot!(epoch_idx[]))
    
    _create_navigation_handlers(controls, epoch_idx, artifact_idx, n_epochs, epochs_with_artifacts)

    # Update plot and controls when epoch changes
    on(epoch_idx) do idx
        update_plot!(idx)
        controls.epoch_label.text = "Epoch $idx / $n_epochs"
    end

    # Update plot and controls when artifact changes
    on(artifact_idx) do idx
        if !isempty(epochs_with_artifacts)
            epoch_idx[] = epochs_with_artifacts[idx]
            controls.artifact_label.text = "Artifact $idx / $(length(epochs_with_artifacts))"
        end
    end

    # Initialize plot and controls
    update_plot!(epoch_idx[])
    controls.epoch_label.text = "Epoch 1 / $n_epochs"
    controls.artifact_label.text = !isempty(epochs_with_artifacts) ? 
        "Artifact 1 / $(length(epochs_with_artifacts))" : "No artifacts found"

    plot_kwargs[:display_plot] && display_figure(fig)
    return fig
end

"""
    plot_artifact_repair(epochs_original::EpochData, epochs_repaired::EpochData, artifacts::EpochRejectionInfo; channel_selection::Function=channels(), kwargs...)

Interactive plot comparison between original and repaired epochs with navigation buttons.

# Arguments
- `epochs_original::EpochData`: Original epoch data
- `epochs_repaired::EpochData`: Repaired epoch data  
- `artifacts::EpochRejectionInfo`: Artifact detection results
- `channel_selection::Function`: Channel predicate for selecting channels to plot (default: all layout channels)

# Keyword Arguments
$(generate_kwargs_doc(PLOT_ARTIFACT_KWARGS))

# Returns
- `Figure`: Interactive Makie figure with navigation buttons showing before/after comparison

# Examples
```julia
# Basic interactive repair comparison
fig = plot_artifact_repair(epochs_orig, epochs_repaired, artifacts)

# With specific channels
fig = plot_artifact_repair(epochs_orig, epochs_repaired, artifacts,
                          channel_selection = channels([:Fp1, :Fp2, :Cz]))

# With custom axis limits
fig = plot_artifact_repair(epochs_orig, epochs_repaired, artifacts, xlimits = (-0.2, 0.8), ylimits = (-100, 100))
```
"""
function plot_artifact_repair(
    epochs_original::EpochData,
    epochs_repaired::EpochData,
    artifacts::EpochRejectionInfo;
    channel_selection::Function = channels(),
    kwargs...,
)
    # Merge user kwargs with defaults and validate
    plot_kwargs = _merge_plot_kwargs(PLOT_ARTIFACT_KWARGS, kwargs)
    
    # Get channels to plot
    selected_channels = get_selected_channels(epochs_original, channel_selection, include_meta = false, include_extra = false)

    # Create figure and axes
    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)")
    ax2 = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)")

    # Setup axis styling
    _setup_axis_grid!(ax1, plot_kwargs)
    _setup_axis_limits!(ax1, plot_kwargs)
    _setup_origin_lines!(ax1, plot_kwargs)
    
    _setup_axis_grid!(ax2, plot_kwargs)
    _setup_axis_limits!(ax2, plot_kwargs)
    _setup_origin_lines!(ax2, plot_kwargs)

    # Get epochs with artifacts
    epochs_with_artifacts = unique([r.epoch for r in artifacts.rejected_epochs])
    sort!(epochs_with_artifacts)

    # Create controls
    n_epochs = length(epochs_original.data)
    controls = _create_artifact_controls(fig[3, 1], n_epochs, epochs_with_artifacts)

    # Set row sizes
    rowsize!(fig.layout, 1, Relative(0.45))  # Original plot
    rowsize!(fig.layout, 2, Relative(0.45))  # Repaired plot
    rowsize!(fig.layout, 3, Relative(0.1))   # Controls area

    # Create observables
    epoch_idx = Observable(1)
    artifact_idx = Observable(1)
    current_legend1 = Ref{Union{Nothing,Legend}}(nothing)
    current_legend2 = Ref{Union{Nothing,Legend}}(nothing)
    selected_channels_set = Set{Symbol}()

    # Create color map for rejected channels
    rejected_channels = [r.label for r in artifacts.rejected_epochs]
    rejected_color_map = _create_rejected_color_map(rejected_channels, plot_kwargs[:colormap_name])

    # Function to update comparison plot
    function update_comparison_plot!(epoch_idx_val)

        empty!(ax1)
        empty!(ax2)

        # Delete existing legends if present
        if current_legend1[] !== nothing
            delete!(current_legend1[])
            current_legend1[] = nothing
        end
        if current_legend2[] !== nothing
            delete!(current_legend2[])
            current_legend2[] = nothing
        end

        # Get epoch data
        epoch_orig = epochs_original.data[epoch_idx_val]
        epoch_repaired = epochs_repaired.data[epoch_idx_val]

        # Find rejected channels for this epoch
        epoch_rejected_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx_val]

        # Update titles
        ax1.title = "Original Data - Epoch $(epoch_idx_val)"
        ax2.title = "Repaired Data - Epoch $(epoch_idx_val)"

        # Plot channels for both axes
        _plot_channels!(ax1, epoch_orig, selected_channels, epoch_rejected_channels, selected_channels_set, 
                      rejected_color_map, plot_kwargs)
        
        # For repaired plot, show which channels were originally rejected (now repaired)
        _plot_channels!(ax2, epoch_repaired, selected_channels, epoch_rejected_channels, selected_channels_set, 
                      rejected_color_map, plot_kwargs; is_repaired=true)

        # Create legends
        _create_legend!(ax1, selected_channels, current_legend1)
        _create_legend!(ax2, selected_channels, current_legend2)

        # Link axes
        linkxaxes!(ax1, ax2)
    end

    # Create event handlers for both axes
    _create_channel_selection_handlers(fig, ax1, epochs_original, selected_channels, selected_channels_set, 
                                     plot_kwargs, epoch_idx, () -> update_comparison_plot!(epoch_idx[]))
    _create_channel_selection_handlers(fig, ax2, epochs_original, selected_channels, selected_channels_set, 
                                     plot_kwargs, epoch_idx, () -> update_comparison_plot!(epoch_idx[]))
    

    # Create event handlers
    _create_navigation_handlers(controls, epoch_idx, artifact_idx, n_epochs, epochs_with_artifacts)

    # Update plot and controls when epoch changes
    on(epoch_idx) do idx
        update_comparison_plot!(idx)
        controls.epoch_label.text = "Epoch $idx / $n_epochs"
    end

    # Update plot and controls when artifact changes
    on(artifact_idx) do idx
        if !isempty(epochs_with_artifacts)
            epoch_idx[] = epochs_with_artifacts[idx]
            controls.artifact_label.text = "Artifact $idx / $(length(epochs_with_artifacts))"
        end
    end

    # Initialize plot and controls
    update_comparison_plot!(epoch_idx[])
    controls.epoch_label.text = "Epoch 1 / $n_epochs"
    controls.artifact_label.text = !isempty(epochs_with_artifacts) ? 
        "Artifact 1 / $(length(epochs_with_artifacts))" : "No artifacts found"

    plot_kwargs[:display_plot] && display_figure(fig)
    return fig
end
