# =============================================================================
# CONTROL PANEL FOR PLOT_ERP
# =============================================================================

"""
    _setup_erp_control_panel!(fig::Figure, datasets::Vector{ErpData}, axes::Vector{Axis}, 
                               plot_layout::PlotLayout, plot_kwargs::Dict, 
                               channel_selection::Function, sample_selection::Function,
                               baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real},Nothing})

Set up a control panel that opens when 'c' key is pressed.
Allows adjusting baseline and toggling conditions.
"""
function _setup_erp_control_panel!(
    fig::Figure,
    dat_subset::Vector{ErpData},
    axes::Vector{Axis},
    plot_layout::PlotLayout,
    plot_kwargs::Dict,
    channel_selection::Function,
    sample_selection::Function,
    baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real},Nothing},
    initial_legends::Dict{Int, Union{Legend, Nothing}} = Dict{Int, Union{Legend, Nothing}}(),
)
    # Track if control panel has been created
    control_fig = Ref{Union{Figure,Nothing}}(nothing)

    # Get time range from first dataset for baseline limits
    time_vec = first(dat_subset).data.time
    time_min = minimum(time_vec)
    time_max = maximum(time_vec)

    # Current baseline values (observables)
    baseline_start_val = baseline_interval !== nothing && baseline_interval isa Tuple ? baseline_interval[1] : nothing
    baseline_stop_val = baseline_interval !== nothing && baseline_interval isa Tuple ? baseline_interval[2] : nothing
    baseline_start = Observable(baseline_start_val !== nothing ? string(baseline_start_val) : "")
    baseline_stop = Observable(baseline_stop_val !== nothing ? string(baseline_stop_val) : "")
    
    # Condition checkboxes (all checked by default)
    condition_checked = [Observable(true) for _ in dat_subset]

    # Store original channel list from plot_layout
    layout_channels = plot_layout.channels
    
    # Store line references organized by dataset index and axis
    # line_refs[dataset_idx][ax_idx] = Vector of Lines objects for that dataset on that axis
    line_refs = [Dict{Int, Vector{Lines}}() for _ in dat_subset]
    
    # Store original colors for each dataset (by original dataset index)
    # original_colors[dataset_idx] = color for that dataset
    original_colors = Vector{Any}(undef, length(dat_subset))
    original_colors_set = false
    
    # Store legend references for each axis
    # legend_refs[ax_idx] = Legend object for that axis
    # Start with the initial legends passed from plot_erp
    legend_refs = copy(initial_legends)
    
    # Store text box references for baseline inputs (needed in update_plot!)
    baseline_start_input_ref = Ref{Union{Textbox,Nothing}}(nothing)
    baseline_stop_input_ref = Ref{Union{Textbox,Nothing}}(nothing)
    
    # Function to collect legend references from axes (fallback if initial_legends is empty)
    function collect_legend_refs!()
        # Only collect if we don't already have references
        if isempty(legend_refs)
            for (ax_idx, ax) in enumerate(axes)
                try
                    if hasproperty(ax.scene, :children) && ax.scene.children isa Dict
                        for child in values(ax.scene.children)
                            if child isa Legend
                                legend_refs[ax_idx] = child
                                break  # Only one legend per axis
                            end
                        end
                    end
                catch
                end
            end
        end
    end
    
    # Function to collect line references from axes
    function collect_line_refs!()
        # Clear existing references
        for refs in line_refs
            empty!(refs)
        end
        
        for (ax_idx, ax) in enumerate(axes)
            try
                # Try multiple ways to access scene children
                scene_children = nothing
                if hasproperty(ax.scene, :children)
                    scene_children = ax.scene.children
                elseif hasproperty(ax, :scene) && hasproperty(ax.scene, :plots)
                    scene_children = ax.scene.plots
                end
                
                if scene_children !== nothing
                    children_to_check = []
                    if scene_children isa Dict
                        children_to_check = collect(values(scene_children))
                    elseif scene_children isa AbstractVector
                        children_to_check = scene_children
                    end
                    
                    for child in children_to_check
                        if child isa Lines && !(child isa Legend)
                            # Try to match line to dataset by checking its label
                            # Lines are created with labels like "condition_name" or "condition_name (channel)"
                            label_str = nothing
                            if hasproperty(child, :label)
                                label_val = child.label
                                if label_val isa Observable
                                    label_str = label_val[]
                                else
                                    label_str = string(label_val)
                                end
                            end
                            
                            if label_str !== nothing
                                for (dataset_idx, dat) in enumerate(dat_subset)
                                    expected_label_start = dat.condition_name
                                    # Check if this line belongs to this dataset
                                    if startswith(label_str, expected_label_start)
                                        if !haskey(line_refs[dataset_idx], ax_idx)
                                            line_refs[dataset_idx][ax_idx] = []
                                        end
                                        push!(line_refs[dataset_idx][ax_idx], child)
                                        
                                        # Store the original color for this dataset (only from first channel found)
                                        # We want the color from the first channel to use for all channels of this dataset
                                        if !original_colors_set && !isassigned(original_colors, dataset_idx) && hasproperty(child, :color)
                                            color_val = child.color
                                            if color_val isa Observable
                                                stored_color = color_val[]
                                            else
                                                stored_color = color_val
                                            end
                                            # Convert to RGB if it's a color type, to ensure consistency
                                            if stored_color isa ColorTypes.Colorant
                                                original_colors[dataset_idx] = stored_color
                                            else
                                                original_colors[dataset_idx] = stored_color
                                            end
                                        end
                                        break  # Found match, move to next line
                                    end
                                end
                            end
                        end
                    end
                end
            catch e
                @debug "Error collecting line references from axis $ax_idx: $e"
            end
        end
    end
    
    # Collect line references and legend references after initial plot
    collect_line_refs!()
    collect_legend_refs!()
    original_colors_set = true
    
    # Debug: print how many line references we found and colors stored
        @debug "Collected line references: $(sum(length, line_refs)) total lines across $(length(dat_subset)) datasets"
    @debug "Collected legend references: $(length(legend_refs)) legends"
    @debug "Stored colors: $(sum(i -> isassigned(original_colors, i) && original_colors[i] !== nothing, 1:length(original_colors))) out of $(length(original_colors))"
    
    # Function to update line visibility based on checkboxes
    function update_visibility!()
        for (dataset_idx, checked_obs) in enumerate(condition_checked)
            visible = checked_obs[]
            # Only update if we have line references for this dataset
            if length(line_refs[dataset_idx]) > 0
                for (ax_idx, lines_vec) in line_refs[dataset_idx]
                    for line in lines_vec
                        try
                            if visible
                                line.visible = true
                            else
                                line.visible = false
                            end
                        catch e
                            # If line was deleted or visibility update fails, try hide/show
                            try
                                if visible
                                    show!(line)
                                else
                                    hide!(line)
                                end
                            catch e2
                                @debug "Could not update visibility for line: $e2"
                            end
                        end
                    end
                end
            else
                @debug "No line references found for dataset $dataset_idx (condition: $(dat_subset[dataset_idx].condition_name))"
            end
        end
    end
    
    # Function to update the plot (for baseline changes)
    function update_plot!()
        try
            @info "update_plot! called"
            # Parse baseline interval
            # Read directly from text boxes if available, otherwise use observables
            baseline_interval_new = nothing
            if baseline_start_input_ref[] !== nothing && baseline_stop_input_ref[] !== nothing
                start_input = baseline_start_input_ref[]
                stop_input = baseline_stop_input_ref[]
                
                # Try to get the current text value
                # stored_string updates on Enter/blur, but we want current displayed text
                # Try multiple properties to get the current value
                start_str = ""
                stop_str = ""
                
                # First try stored_string (what user has confirmed with Enter)
                if hasproperty(start_input, :stored_string)
                    start_str = start_input.stored_string[]
                end
                if hasproperty(stop_input, :stored_string)
                    stop_str = stop_input.stored_string[]
                end
                
                # If stored_string is empty, try displayed_string or text property
                if start_str == "" && hasproperty(start_input, :displayed_string)
                    start_str = start_input.displayed_string[]
                elseif start_str == "" && hasproperty(start_input, :text)
                    start_str = start_input.text[]
                end
                
                if stop_str == "" && hasproperty(stop_input, :displayed_string)
                    stop_str = stop_input.displayed_string[]
                elseif stop_str == "" && hasproperty(stop_input, :text)
                    stop_str = stop_input.text[]
                end
                
                @info "Reading from textboxes: start='$start_str', stop='$stop_str'"
            else
                # Fallback to observables if text boxes not yet created
                start_str = baseline_start[]
                stop_str = baseline_stop[]
                @info "Textboxes not created yet, using observables: start='$start_str', stop='$stop_str'"
            end
            
            if start_str != "" && stop_str != ""
                try
                    start_val = parse(Float64, start_str)
                    stop_val = parse(Float64, stop_str)
                    baseline_interval_new = (start_val, stop_val)
                    @info "Applying baseline: $baseline_interval_new"
                catch e
                    @warn "Invalid baseline values: start='$start_str', stop='$stop_str' - $e"
                    return  # Don't update if parsing fails
                end
            else
                @info "No baseline specified (empty text boxes): start='$start_str', stop='$stop_str'"
            end

            # Get selected conditions
            selected_conditions = [i for (i, checked) in enumerate(condition_checked) if checked[]]
            if isempty(selected_conditions)
                @warn "No conditions selected"
                return  # Don't update if no conditions selected
            end
            @info "Selected conditions: $selected_conditions"
            
            # Filter by selected conditions and make a copy
            dat_subset_selected = [copy(dat) for dat in dat_subset[selected_conditions]]
            
            # Apply new baseline if specified (apply directly to dat_subset)
            if baseline_interval_new !== nothing
                @info "Applying baseline! to $(length(dat_subset_selected)) datasets with interval $baseline_interval_new"
                try
                    baseline!.(dat_subset_selected, Ref(baseline_interval_new))
                    @info "Baseline applied successfully"
                catch e
                    @error "Error applying baseline: $e"
                    rethrow(e)
                end
            end
            
            # Extract channel labels (same as original plot)
            selected_channels = channel_labels(dat_subset_selected)
            extra_channels = extra_labels(dat_subset_selected)
            all_plot_channels = vcat(selected_channels, extra_channels)

            # Delete old legends BEFORE re-plotting (so they don't get recreated)
            for (ax_idx, ax) in enumerate(axes)
                # Search for ALL legends in the scene and collect them
                legends_to_delete = []
                try
                    if hasproperty(ax.scene, :children) && ax.scene.children isa Dict
                        for (key, child) in ax.scene.children
                            if child isa Legend
                                push!(legends_to_delete, child)
                            end
                        end
                    end
                catch
                end
                
                # Also add stored reference if we have one
                if haskey(legend_refs, ax_idx) && legend_refs[ax_idx] !== nothing
                    leg_ref = legend_refs[ax_idx]
                    if !(leg_ref in legends_to_delete)
                        push!(legends_to_delete, leg_ref)
                    end
                end
                
                # Delete all legends we found
                for leg in legends_to_delete
                    try
                        delete!(leg)
                    catch
                    end
                end
                
                # Clear stored reference
                if haskey(legend_refs, ax_idx)
                    legend_refs[ax_idx] = nothing
                end
            end

            # Clear all plot elements from axes (but keep axis structure)
            for ax in axes
                empty!(ax)
            end
            
            # Clear line references (will be recollected)
            for refs in line_refs
                empty!(refs)
            end

            # Re-plot WITHOUT creating new legend, using original colors
            # Map selected datasets to their original colors (preserving original order/colors)
            # Start with a copy of plot_kwargs and override specific values
            plot_kwargs_no_legend = copy(plot_kwargs)
            plot_kwargs_no_legend[:legend] = false
            
            # Build color array matching the structure expected by _plot_erp!
            # Colors are indexed as: (dataset_idx - 1) * n_channels + channel_idx
            if original_colors_set
                # Check if all selected conditions have valid colors
                valid_colors = true
                for i in selected_conditions
                    if !(1 <= i <= length(original_colors)) || !isassigned(original_colors, i) || original_colors[i] === nothing
                        valid_colors = false
                        break
                    end
                end
                
                if valid_colors
                    n_channels = length(all_plot_channels)
                    expanded_colors = []
                    for orig_idx in selected_conditions
                        orig_color = original_colors[orig_idx]
                        # Add this color for each channel (preserving original color per dataset)
                        # This ensures all channels of the same dataset use the same color
                        for _ in 1:n_channels
                            push!(expanded_colors, orig_color)
                        end
                    end
                    # Override color and mark as explicitly set to prevent recomputation
                    plot_kwargs_no_legend[:color] = expanded_colors
                    plot_kwargs_no_legend[:_color_explicitly_set] = true
                    @debug "Using stored colors: $(length(expanded_colors)) colors for $(length(selected_conditions)) datasets with $n_channels channels"
                    @debug "Color values: $(expanded_colors[1:min(4, length(expanded_colors))])"
                else
                    @debug "Not using stored colors - some missing or invalid"
                end
            end
            
            @info "Re-plotting with $(length(dat_subset_selected)) datasets"
            if plot_layout.type == :single
                ax_returned, leg = _plot_erp!(axes[1], dat_subset_selected, all_plot_channels; fig=fig, plot_kwargs_no_legend...)
                if leg !== nothing
                    legend_refs[1] = leg
                end
                @info "Re-plotted single layout"
            else
                for (ax_idx, (ax, channel)) in enumerate(zip(axes, layout_channels))
                    if channel in all_plot_channels
                        # For grid/topo with single channel, colors are just per dataset
                        if original_colors_set
                            # Check if all selected conditions have valid colors
                            valid_colors = true
                            for i in selected_conditions
                                if !(1 <= i <= length(original_colors)) || !isassigned(original_colors, i) || original_colors[i] === nothing
                                    valid_colors = false
                                    break
                                end
                            end
                            
                            if valid_colors
                                channel_colors = [original_colors[i] for i in selected_conditions]
                                plot_kwargs_no_legend[:color] = channel_colors
                                plot_kwargs_no_legend[:_color_explicitly_set] = true
                            end
                        end
                        ax_returned, leg = _plot_erp!(ax, dat_subset_selected, [channel]; fig=fig, plot_kwargs_no_legend...)
                        if leg !== nothing
                            legend_refs[ax_idx] = leg
                        end
                    end
                end
                @info "Re-plotted grid/topo layout"
            end
            
            # Collect line references again after re-plotting
            collect_line_refs!()
            
            # Update visibility based on checkboxes
            update_visibility!()
            
            # Update legend to match current plot (create new with current datasets)
            for (ax_idx, ax) in enumerate(axes)
                # Create new legend with current datasets and channels
                if haskey(legend_refs, ax_idx) && legend_refs[ax_idx] !== nothing
                    # Legend was already created by _plot_erp! and stored
                    # No need to recreate
                else
                    # Create legend if needed
                    channels_to_plot = plot_layout.type == :single ? all_plot_channels : [layout_channels[ax_idx]]
                    leg = _add_legend!(ax, channels_to_plot, dat_subset_selected, plot_kwargs)
                    if leg !== nothing
                        legend_refs[ax_idx] = leg
                    end
                end
            end
            
            @info "Plot update complete"
        catch e
            @error "Error in update_plot!: $e" exception=(e, catch_backtrace())
        end
    end

    # Keyboard handler for 'c' key
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.c
            # Only create if it doesn't exist yet
            if control_fig[] === nothing
                # Create control panel
                control_fig[] = Figure(title = "ERP Control Panel", size = (300, 400))
                control_layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10, colgap = 10)

                # Baseline section
                Label(control_layout[1, 1], "Baseline Correction", fontsize = 14, font = :bold, tellwidth = false)
                
                baseline_layout = GridLayout(control_layout[2, 1], tellwidth = false, colgap = 10)
                Label(baseline_layout[1, 1], "Start:", width = 60, tellwidth = false)
                baseline_start_input = Textbox(baseline_layout[1, 2], placeholder = "e.g. -0.2", width = 100)
                baseline_start_input.stored_string[] = baseline_start[]
                baseline_start_input_ref[] = baseline_start_input  # Store reference for update_plot!
                
                Label(baseline_layout[2, 1], "Stop:", width = 60, tellwidth = false)
                baseline_stop_input = Textbox(baseline_layout[2, 2], placeholder = "e.g. 0.0", width = 100)
                baseline_stop_input.stored_string[] = baseline_stop[]
                baseline_stop_input_ref[] = baseline_stop_input  # Store reference for update_plot!

                # Connect text inputs to observables
                on(baseline_start_input.stored_string) do val
                    baseline_start[] = val
                end
                on(baseline_stop_input.stored_string) do val
                    baseline_stop[] = val
                end

                # Conditions section
                Label(control_layout[3, 1], "Conditions", fontsize = 14, font = :bold, tellwidth = false)
                
                conditions_layout = GridLayout(control_layout[4, 1], tellwidth = false, rowgap = 5)
                condition_checkboxes = Checkbox[]
                for (i, dat) in enumerate(dat_subset)
                    cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][], tellwidth = false)
                    Label(conditions_layout[i, 2], dat.condition_name, tellwidth = false)
                    connect!(condition_checked[i], cb.checked)
                    push!(condition_checkboxes, cb)
                end

                # Apply button
                apply_button = Button(control_layout[5, 1], label = "Apply", width = 200)
                on(apply_button.clicks) do _
                    # Force textboxes to commit their current values to stored_string
                    if baseline_start_input_ref[] !== nothing
                        # Try to get current displayed value and commit it
                        start_input = baseline_start_input_ref[]
                        if hasproperty(start_input, :displayed_string)
                            start_input.stored_string[] = start_input.displayed_string[]
                        elseif hasproperty(start_input, :text)
                            start_input.stored_string[] = start_input.text[]
                        end
                    end
                    if baseline_stop_input_ref[] !== nothing
                        stop_input = baseline_stop_input_ref[]
                        if hasproperty(stop_input, :displayed_string)
                            stop_input.stored_string[] = stop_input.displayed_string[]
                        elseif hasproperty(stop_input, :text)
                            stop_input.stored_string[] = stop_input.text[]
                        end
                    end
                    update_plot!()
                end

                # Auto-update visibility on condition checkbox changes (no re-plot needed!)
                for checked in condition_checked
                    on(checked) do _
                        update_visibility!()
                    end
                end

                # Display the control panel
                display(control_fig[])
            end
        end
    end
end

