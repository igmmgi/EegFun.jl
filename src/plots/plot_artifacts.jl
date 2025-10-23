# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ARTIFACT_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :colormap_name => (:Set1_9, "Colormap name for artifact channel colors (e.g., :Set1_9, :tab20)."),
    :display_plot => (true, "Whether to display the plot."),
    :plot_height => (600, "Height of the plot area in pixels."),
    :controls_height => (80, "Height of the controls area in pixels."),
    :button_width => (100, "Width of navigation buttons."),
    :artifact_button_width => (130, "Width of artifact navigation buttons."),
    :label_width => (100, "Width of epoch/artifact labels."),
    :artifact_label_width => (130, "Width of artifact label."),
    :linewidth_normal => (1, "Line width for normal channels."),
    :linewidth_rejected => (2, "Line width for rejected channels."),
    :linewidth_selected => (3, "Line width for selected channels."),
    :linewidth_selected_rejected => (4, "Line width for selected rejected channels."),
    :alpha_normal => (0.2, "Transparency for normal channels."),
    :alpha_rejected => (1.0, "Transparency for rejected channels."),
    :legend_position => ((:right, :top), "Position of the legend."),
    :legend_nbanks => (5, "Number of columns in legend for many channels."),
    :selection_threshold => (50, "Distance threshold for channel selection via mouse click."),
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
fig = plot_artifact_detection_interactive(epochs, artifacts)

# With specific channels
fig = plot_artifact_detection_interactive(epochs, artifacts, 
                                        channel_selection = channels([:Fp1, :Fp2, :Cz]))
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

    # Create figure/ax with space for controls
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)")

    # Create horizontal layout for controls 
    n_epochs = length(epochs.data)
    controls_layout = GridLayout(fig[2, 1], tellheight = false, tellwidth = false)

    # Get epochs with artifacts
    epochs_with_artifacts = unique([r.epoch for r in artifacts.rejected_epochs])
    sort!(epochs_with_artifacts)

    # Add controls with artifact navigation
    back_button = Button(controls_layout[1, 1], label = "◀ Previous", width = plot_kwargs[:button_width])
    epoch_label =
        Label(controls_layout[1, 2], "Epoch 1 / $n_epochs", width = plot_kwargs[:button_width], justification = :center)
    forward_button = Button(controls_layout[1, 3], label = "Next ▶", width = plot_kwargs[:button_width])

    # Add artifact navigation buttons
    back_artifact_button =
        Button(controls_layout[1, 4], label = "◀ Previous Artifact", width = plot_kwargs[:artifact_button_width])
    artifact_label = Label(
        controls_layout[1, 5],
        "Artifact 1 / $(length(epochs_with_artifacts))",
        width = plot_kwargs[:artifact_button_width],
        justification = :center,
    )
    forward_artifact_button =
        Button(controls_layout[1, 6], label = "Next Artifact ▶", width = plot_kwargs[:artifact_button_width])

    # Center the entire control group
    for i = 1:6
        colsize!(controls_layout, i, Auto())
    end

    # Set row sizes
    rowsize!(controls_layout, 1, Auto())
    rowsize!(fig.layout, 1, Relative(0.9))  # Plot area
    rowsize!(fig.layout, 2, Relative(0.1))  # Controls area

    # Create observables for epoch and artifact indices
    epoch_idx = Observable(1)
    artifact_idx = Observable(1)

    # Keep reference to legend for cleanup
    current_legend = Ref{Union{Nothing,Legend}}(nothing)

    # Track selected channels for highlighting
    selected_channels_set = Set{Symbol}()

    # Function to update plot based on epoch
    function update_plot!(ax, epoch_idx_val)
        # Clear the axis
        empty!(ax)

        # Delete existing legend if present
        if current_legend[] !== nothing
            delete!(current_legend[])
            current_legend[] = nothing
        end

        # Get current epoch data
        epoch = epochs.data[epoch_idx_val]
        time_points = epoch.time

        # Find rejected channels for this epoch
        rejected_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx_val]

        # Create color mapping for rejected channels using Makie categorical colors
        rejected_color_map = Dict{Symbol,Any}()
        if !isempty(rejected_channels)
            # Use a colormap that cycles for any number of channels
            colormap_name = plot_kwargs[:colormap_name]
            # Get the actual number of colors in the colormap
            colormap = cgrad(colormap_name)
            max_colors = length(colormap.colors)

            for (i, ch) in enumerate(rejected_channels)
                color_idx = ((i - 1) % max_colors) + 1
                rejected_color_map[ch] = colormap.colors[color_idx]
            end
        end

        # Update title
        ax.title = "Artifact Detection - Epoch $(epoch_idx_val)"

        # Configure grid and axes
        ax.xgridvisible = plot_kwargs[:xgrid]
        ax.ygridvisible = plot_kwargs[:ygrid]
        ax.xminorgridvisible = plot_kwargs[:xminorgrid]
        ax.yminorgridvisible = plot_kwargs[:yminorgrid]

        # Add origin lines if requested
        if plot_kwargs[:axes_through_origin]
            hlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 1)
            vlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 1)
        end

        # Plot each channel
        for ch in selected_channels
            if hasproperty(epoch, ch)
                # Determine styling based on rejection and selection status
                is_rejected = ch in rejected_channels
                is_selected = ch in selected_channels_set

                color = is_rejected ? rejected_color_map[ch] : :black
                alpha = is_rejected ? plot_kwargs[:alpha_rejected] : plot_kwargs[:alpha_normal]
                linewidth =
                    is_selected ?
                    is_rejected ? plot_kwargs[:linewidth_selected_rejected] : plot_kwargs[:linewidth_selected] :
                    is_rejected ? plot_kwargs[:linewidth_rejected] : plot_kwargs[:linewidth_normal]

                # Create label with bold formatting for selected channels
                label_text = is_rejected ? "$ch (rejected)" : "$ch"
                label_text = is_selected ? rich(label_text, font = :bold) : label_text

                lines!(
                    ax,
                    time_points,
                    epoch[!, ch],
                    color = color,
                    alpha = alpha,
                    linewidth = linewidth,
                    label = label_text,
                )
            end
        end

        # Add legend and store reference with multiple columns for many channels
        n_channels = length(selected_channels)
        n_cols = n_channels > 10 ? cld(n_channels, 20) : 1
        current_legend[] =
            axislegend(ax, position = plot_kwargs[:legend_position], nbanks = plot_kwargs[:legend_nbanks])
    end

    # Add shift+click functionality to select/deselect channels
    on(events(ax).mousebutton, priority = 0) do event
        if event.button == Mouse.left && event.action == Mouse.press
            # Check if shift is held
            if Keyboard.left_shift in events(fig).keyboardstate || Keyboard.right_shift in events(fig).keyboardstate
                # Get mouse position in data coordinates
                pos = mouseposition(ax.scene)
                if pos[1] !== nothing && pos[2] !== nothing
                    mouse_time = pos[1]
                    mouse_amp = pos[2]

                    # Get current epoch data
                    current_epoch = epochs.data[epoch_idx[]]
                    time_points = current_epoch.time

                    # Find the closest time point
                    if !isempty(time_points)
                        closest_time_idx = argmin(abs.(time_points .- mouse_time))

                        # Check which channel line is closest to mouse position
                        min_distance = Inf
                        closest_channel = nothing

                        for ch in selected_channels
                            if hasproperty(current_epoch, ch)
                                channel_amp = current_epoch[closest_time_idx, ch]
                                distance = abs(channel_amp - mouse_amp)
                                if distance < min_distance
                                    min_distance = distance
                                    closest_channel = ch
                                end
                            end
                        end

                        # Toggle selection if we found a close channel
                        if closest_channel !== nothing && min_distance < 50  # Threshold for selection
                            if closest_channel in selected_channels_set
                                delete!(selected_channels_set, closest_channel)
                                @info "Deselected channel: $closest_channel"
                            else
                                push!(selected_channels_set, closest_channel)
                                @info "Selected channel: $closest_channel"
                            end

                            # Refresh the plot
                            update_plot!(ax, epoch_idx[])
                        end
                    end
                end
            end
        end
        return Consume(false)
    end

    # Button click handlers
    on(back_button.clicks) do _
        ;
        epoch_idx[] = max(1, epoch_idx[] - 1)
    end
    on(forward_button.clicks) do _
        ;
        epoch_idx[] = min(n_epochs, epoch_idx[] + 1)
    end
    on(back_artifact_button.clicks) do _
        ;
        artifact_idx[] = max(1, artifact_idx[] - 1)
    end
    on(forward_artifact_button.clicks) do _
        ;
        artifact_idx[] = min(length(epochs_with_artifacts), artifact_idx[] + 1)
    end

    # Update plot and controls when epoch changes
    on(epoch_idx) do idx
        update_plot!(ax, idx)
        epoch_label.text = "Epoch $idx / $n_epochs"
    end

    # Update plot and controls when artifact changes
    on(artifact_idx) do idx
        if !isempty(epochs_with_artifacts)
            epoch_idx[] = epochs_with_artifacts[idx]
            artifact_label.text = "Artifact $idx / $(length(epochs_with_artifacts))"
        end
    end

    # Initialize plot and controls
    update_plot!(ax, epoch_idx[])
    epoch_label.text = "Epoch 1 / $n_epochs"
    artifact_label.text =
        !isempty(epochs_with_artifacts) ? "Artifact 1 / $(length(epochs_with_artifacts))" : "No artifacts found"

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
fig = plot_repair_comparison_interactive(epochs_orig, epochs_repaired, artifacts)

# With specific channels
fig = plot_repair_comparison_interactive(epochs_orig, epochs_repaired, artifacts,
                                       channel_selection = channels([:Fp1, :Fp2, :Cz]))
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
    selected_channels =
        get_selected_channels(epochs_original, channel_selection, include_meta = false, include_extra = false)

    # Create figure with space for controls
    fig = Figure()

    # Create axes for original and repaired data (full width)
    ax1 = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)")
    ax2 = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Amplitude (μV)")

    # Create horizontal layout for controls
    n_epochs = length(epochs_original.data)
    controls_layout = GridLayout(fig[3, 1], tellheight = false, tellwidth = false)

    # Get epochs with artifacts
    epochs_with_artifacts = unique([r.epoch for r in artifacts.rejected_epochs])
    sort!(epochs_with_artifacts)

    # Add controls with artifact navigation
    back_button = Button(controls_layout[1, 1], label = "◀ Previous", width = plot_kwargs[:button_width])
    epoch_label =
        Label(controls_layout[1, 2], "Epoch 1 / $n_epochs", width = plot_kwargs[:button_width], justification = :center)
    forward_button = Button(controls_layout[1, 3], label = "Next ▶", width = plot_kwargs[:button_width])

    # Add artifact navigation buttons
    back_artifact_button =
        Button(controls_layout[1, 4], label = "◀ Previous Artifact", width = plot_kwargs[:artifact_button_width])
    artifact_label = Label(
        controls_layout[1, 5],
        "Artifact 1 / $(length(epochs_with_artifacts))",
        width = plot_kwargs[:artifact_button_width],
        justification = :center,
    )
    forward_artifact_button =
        Button(controls_layout[1, 6], label = "Next Artifact ▶", width = plot_kwargs[:artifact_button_width])

    # Center the entire control group
    for i = 1:6
        colsize!(controls_layout, i, Auto())
    end

    # Set row sizes
    rowsize!(controls_layout, 1, Auto())
    rowsize!(fig.layout, 1, Relative(0.45))  # Original plot
    rowsize!(fig.layout, 2, Relative(0.45))  # Repaired plot
    rowsize!(fig.layout, 3, Relative(0.1))   # Controls area

    # Create observables for epoch and artifact indices
    epoch_idx = Observable(1)
    artifact_idx = Observable(1)

    # Keep references to legends for cleanup
    current_legend1 = Ref{Union{Nothing,Legend}}(nothing)
    current_legend2 = Ref{Union{Nothing,Legend}}(nothing)

    # Track selected channels for highlighting
    selected_channels_set = Set{Symbol}()

    # Function to update plot based on epoch
    function update_comparison_plot!(ax1, ax2, epoch_idx_val)
        # Clear both axes
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
        time_points = epoch_orig.time

        # Find rejected channels for this epoch
        rejected_channels = [r.label for r in artifacts.rejected_epochs if r.epoch == epoch_idx_val]

        # Create color mapping for rejected channels using Makie categorical colors
        rejected_color_map = Dict{Symbol,Any}()
        if !isempty(rejected_channels)
            # Use a colormap that cycles for any number of channels
            colormap_name = plot_kwargs[:colormap_name]
            # Get the actual number of colors in the colormap
            colormap = cgrad(colormap_name)
            max_colors = length(colormap.colors)

            for (i, ch) in enumerate(rejected_channels)
                color_idx = ((i - 1) % max_colors) + 1
                rejected_color_map[ch] = colormap.colors[color_idx]
            end
        end

        # Update titles
        ax1.title = "Original Data - Epoch $(epoch_idx_val)"
        ax2.title = "Repaired Data - Epoch $(epoch_idx_val)"

        # Configure grid and axes for both plots
        for ax in [ax1, ax2]
            ax.xgridvisible = plot_kwargs[:xgrid]
            ax.ygridvisible = plot_kwargs[:ygrid]
            ax.xminorgridvisible = plot_kwargs[:xminorgrid]
            ax.yminorgridvisible = plot_kwargs[:yminorgrid]

            # Add origin lines if requested
            if plot_kwargs[:axes_through_origin]
                hlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 1)
                vlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 1)
            end
        end

        # Plot each channel
        for ch in selected_channels
            if hasproperty(epoch_orig, ch) && hasproperty(epoch_repaired, ch)
                # Determine styling based on rejection and selection status
                is_rejected = ch in rejected_channels
                is_selected = ch in selected_channels_set

                # Original plot styling
                orig_color = is_rejected ? rejected_color_map[ch] : :black
                orig_alpha = is_rejected ? 1.0 : 0.2
                orig_linewidth = is_selected ? is_rejected ? 4 : 3 : is_rejected ? 2 : 1
                rep_linewidth = is_selected ? 3 : 1

                # Create labels with bold formatting for selected channels
                orig_label_text = is_rejected ? "$ch (rejected)" : "$ch"
                rep_label_text = is_rejected ? "$ch (repaired)" : "$ch"

                if is_selected
                    orig_label_text = rich(orig_label_text, font = :bold)
                    rep_label_text = rich(rep_label_text, font = :bold)
                end

                lines!(
                    ax1,
                    time_points,
                    epoch_orig[!, ch],
                    color = orig_color,
                    alpha = orig_alpha,
                    linewidth = orig_linewidth,
                    label = orig_label_text,
                )

                # Repaired plot
                lines!(
                    ax2,
                    time_points,
                    epoch_repaired[!, ch],
                    color = orig_color,
                    alpha = orig_alpha,
                    linewidth = orig_linewidth,
                    label = rep_label_text,
                )
            end
        end

        # Add legends and store references with multiple columns for many channels
        n_channels = length(selected_channels)
        n_cols = n_channels > 10 ? cld(n_channels, 10) : 1
        current_legend1[] = axislegend(ax1, position = (:right, :top), nbanks = n_cols)
        current_legend2[] = axislegend(ax2, position = (:right, :top), nbanks = n_cols)

        # Link axes
        linkaxes!(ax1, ax2)
    end

    # Add shift+click functionality to both axes for channel selection
    for current_ax in [ax1, ax2]
        on(events(current_ax).mousebutton, priority = 0) do event
            if event.button == Mouse.left && event.action == Mouse.press
                # Check if shift is held
                if Keyboard.left_shift in events(fig).keyboardstate || Keyboard.right_shift in events(fig).keyboardstate
                    # Get mouse position in data coordinates
                    pos = mouseposition(current_ax.scene)
                    if pos[1] !== nothing && pos[2] !== nothing
                        mouse_time = pos[1]
                        mouse_amp = pos[2]

                        # Get current epoch data
                        current_epoch = epochs_original.data[epoch_idx[]]
                        time_points = current_epoch.time

                        # Find the closest time point
                        if !isempty(time_points)
                            closest_time_idx = argmin(abs.(time_points .- mouse_time))

                            # Check which channel line is closest to mouse position
                            min_distance = Inf
                            closest_channel = nothing

                            for ch in selected_channels
                                if hasproperty(current_epoch, ch)
                                    channel_amp = current_epoch[closest_time_idx, ch]
                                    distance = abs(channel_amp - mouse_amp)
                                    if distance < min_distance
                                        min_distance = distance
                                        closest_channel = ch
                                    end
                                end
                            end

                            # Toggle selection if we found a close channel
                            if closest_channel !== nothing && min_distance < plot_kwargs[:selection_threshold]  # Threshold for selection
                                if closest_channel in selected_channels_set
                                    delete!(selected_channels_set, closest_channel)
                                    @info "Deselected channel: $closest_channel"
                                else
                                    push!(selected_channels_set, closest_channel)
                                    @info "Selected channel: $closest_channel"
                                end

                                # Refresh the plot
                                update_comparison_plot!(ax1, ax2, epoch_idx[])
                            end
                        end
                    end
                end
            end
            return Consume(false)
        end
    end

    # Button click handlers
    on(back_button.clicks) do _
        ;
        epoch_idx[] = max(1, epoch_idx[] - 1)
    end
    on(forward_button.clicks) do _
        ;
        epoch_idx[] = min(n_epochs, epoch_idx[] + 1)
    end
    on(back_artifact_button.clicks) do _
        ;
        artifact_idx[] = max(1, artifact_idx[] - 1)
    end
    on(forward_artifact_button.clicks) do _
        ;
        artifact_idx[] = min(length(epochs_with_artifacts), artifact_idx[] + 1)
    end

    # Update plot and controls when epoch changes
    on(epoch_idx) do idx
        update_comparison_plot!(ax1, ax2, idx)
        epoch_label.text = "Epoch $idx / $n_epochs"
    end

    # Update plot and controls when artifact changes
    on(artifact_idx) do idx
        if !isempty(epochs_with_artifacts)
            epoch_idx[] = epochs_with_artifacts[idx]
            artifact_label.text = "Artifact $idx / $(length(epochs_with_artifacts))"
        end
    end

    # Initialize plot and controls
    update_comparison_plot!(ax1, ax2, epoch_idx[])
    epoch_label.text = "Epoch 1 / $n_epochs"

    # Initialize artifact controls
    if !isempty(epochs_with_artifacts)
        artifact_label.text = "Artifact 1 / $(length(epochs_with_artifacts))"
    else
        artifact_label.text = "No artifacts found"
    end

    plot_kwargs[:display_plot] && display_figure(fig)
    return fig
end
