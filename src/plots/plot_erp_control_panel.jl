# =============================================================================
# CONTROL PANEL FOR PLOT_ERP
# =============================================================================

"""
    _setup_erp_control_panel!(fig::Figure, dat_subset::Vector{ErpData}, axes::Vector{Axis}, 
                               plot_layout::PlotLayout, plot_kwargs::Dict,
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
    baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real},Nothing},
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)
    layout_channels = plot_layout.channels
    
    # State: baseline values and condition selections
    baseline_start_obs = Observable(baseline_interval !== nothing && baseline_interval isa Tuple ? string(baseline_interval[1]) : "")
    baseline_stop_obs = Observable(baseline_interval !== nothing && baseline_interval isa Tuple ? string(baseline_interval[2]) : "")
    condition_checked = [Observable(true) for _ in dat_subset]
    
    # Store textbox references
    start_input_ref = Ref{Union{Textbox,Nothing}}(nothing)
    stop_input_ref = Ref{Union{Textbox,Nothing}}(nothing)
    
    # Update plot (re-plot everything with current settings)
    function update_plot!()
        try
            # Get baseline values from textboxes
            baseline_interval_new = nothing
            if start_input_ref[] !== nothing && stop_input_ref[] !== nothing
                start_str = start_input_ref[].stored_string[]
                stop_str = stop_input_ref[].stored_string[]
                if start_str != "" && stop_str != ""
                    try
                        baseline_interval_new = (parse(Float64, start_str), parse(Float64, stop_str))
                    catch e
                        @warn "Invalid baseline values: $e"
                        return
                    end
                end
            end
            
            # Apply baseline if specified
            if baseline_interval_new !== nothing
                baseline!.(dat_subset, Ref(baseline_interval_new))
            end
            dat_to_plot = dat_subset
            
            # Get channels
            selected_channels = channel_labels(dat_to_plot)
            extra_channels = extra_labels(dat_to_plot)
            all_plot_channels = vcat(selected_channels, extra_channels)
            
            # Clear axes
            for ax in axes
                empty!(ax)
            end
            
            # Build condition mask
            condition_mask = [checked[] for checked in condition_checked]
            
            # Re-plot with condition mask
            plot_kwargs_no_legend = merge(copy(plot_kwargs), Dict(:legend => false))
            
            if plot_layout.type == :single
                _plot_erp!(axes[1], dat_to_plot, all_plot_channels; condition_mask=condition_mask, plot_kwargs_no_legend...)
            else
                for (ax, channel) in zip(axes, layout_channels)
                    if channel in all_plot_channels
                        _plot_erp!(ax, dat_to_plot, [channel]; condition_mask=condition_mask, plot_kwargs_no_legend...)
                    end
                end
            end
            
            # Re-apply axis properties
            _apply_axis_properties!.(axes; plot_kwargs...)
            _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...)
        catch e
            @error "Error updating plot: $e" exception=(e, catch_backtrace())
        end
    end
    
    # Keyboard handler for 'c' key
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.c
            # Check if control panel exists and is still open
            if control_fig[] !== nothing
                try
                    # Try to access the scene to see if it's still valid
                    if !isopen(control_fig[].scene)
                        control_fig[] = nothing
                    end
                catch
                    # If we can't access it, it's probably closed
                    control_fig[] = nothing
                end
            end
            
            if control_fig[] === nothing
                control_fig[] = Figure(title = "ERP Control Panel", size = (300, 400))
                layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)
                
                # Baseline section
                Label(layout[1, 1], "Baseline Correction", fontsize = 14, font = :bold)
                baseline_layout = GridLayout(layout[2, 1], tellwidth = false, colgap = 10)
                
                Label(baseline_layout[1, 1], "Start:", width = 60)
                start_input = Textbox(baseline_layout[1, 2], placeholder = "e.g. -0.2", width = 100)
                start_input.stored_string[] = baseline_start_obs[]
                connect!(baseline_start_obs, start_input.stored_string)
                start_input_ref[] = start_input
                
                Label(baseline_layout[2, 1], "Stop:", width = 60)
                stop_input = Textbox(baseline_layout[2, 2], placeholder = "e.g. 0.0", width = 100)
                stop_input.stored_string[] = baseline_stop_obs[]
                connect!(baseline_stop_obs, stop_input.stored_string)
                stop_input_ref[] = stop_input
                
                # Conditions section
                Label(layout[3, 1], "Conditions", fontsize = 14, font = :bold)
                conditions_layout = GridLayout(layout[4, 1], tellwidth = false, rowgap = 5)
                
                for (i, dat) in enumerate(dat_subset)
                    cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][])
                    Label(conditions_layout[i, 2], dat.condition_name)
                    connect!(condition_checked[i], cb.checked)
                end
                
                # Apply button
                apply_btn = Button(layout[5, 1], label = "Apply Baseline", width = 200)
                on(apply_btn.clicks) do _
                    if start_input_ref[] !== nothing && hasproperty(start_input_ref[], :displayed_string)
                        start_input = start_input_ref[]
                        if start_input.displayed_string[] != start_input.stored_string[]
                            start_input.stored_string[] = start_input.displayed_string[]
                        end
                    end
                    if stop_input_ref[] !== nothing && hasproperty(stop_input_ref[], :displayed_string)
                        stop_input = stop_input_ref[]
                        if stop_input.displayed_string[] != stop_input.stored_string[]
                            stop_input.stored_string[] = stop_input.displayed_string[]
                        end
                    end
                    update_plot!()
                end
                
                # Auto-update on condition changes (re-plot)
                for checked in condition_checked
                    on(checked) do _
                        update_plot!()
                    end
                end
                
                display(control_fig[])
            end
        end
    end
end
