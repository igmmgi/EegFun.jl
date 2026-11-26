# Reusable Makie UI controls, BroadcastChannel observables, and menu helpers.

######################################
######### Internal Functions #########
######################################

function _build_checkbox_grid!(grid, items, cols, checked_fn; fontsize = 14, label_fn = string)
    checkboxes = Checkbox[]
    for (i, item) in enumerate(items)
        row = ((i - 1) ÷ cols) + 1
        col_base = ((i - 1) % cols) * 2
        cb = Checkbox(grid[row, col_base+1], checked = checked_fn(item, i), width = 16, height = 16)
        push!(checkboxes, cb)
        Label(grid[row, col_base+2], label_fn(item), fontsize = fontsize, halign = :left)
    end
    actual_cols = min(cols, length(items))
    for c = 1:actual_cols
        col_base = (c - 1) * 2
        colgap!(grid, col_base + 1, 5)
        c < actual_cols && colgap!(grid, col_base + 2, 20)
    end
    rowgap!(grid, 5)
    return checkboxes
end

function _popup_button(fig, label, popup_fn)
    btn = Button(fig, label = label, width = Auto(), fontsize = 18)
    on(btn.clicks) do _
        popup_fn()
    end
    return btn
end

function _add_group_buttons!(fig, grid_row, cbs, items, groups)
    area = fig[grid_row, 1] = GridLayout()
    for (j, (lbl, pred)) in enumerate(groups)
        btn = Button(area[1, j], label = lbl)
        on(btn.clicks) do _
            for (i, cb) in enumerate(cbs)
                cb.checked[] = pred(items[i])
            end
        end
    end
end

function _radio_buttons!(btns::Vector, selected_obs::Observable, values::Vector)
    for (btn, val) in zip(btns, values)
        on(btn.clicks) do _
            selected_obs[] = val
            for (b, v) in zip(btns, values)
                b.buttoncolor[] = v == val ? :lightblue : :white
            end
        end
    end
    for (b, v) in zip(btns, values)
        b.buttoncolor[] = v == selected_obs[] ? :lightblue : :white
    end
end

function _show_checkbox_dialog(
    title_str::String,
    items,
    cols::Int,
    checked_fn,
    groups,
    on_apply!::Function;
    fig_size::Tuple{Int,Int} = (900, 700),
    apply_width::Int = 200,
)
    menu_fig = _eegfun_figure(size = fig_size)
    Label(menu_fig[1, 1], title_str, fontsize = 18, halign = :center)
    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, items, cols, checked_fn)
    !isempty(groups) && _add_group_buttons!(menu_fig, 3, cbs, items, groups)
    action_row = isempty(groups) ? 3 : 4
    action_area = menu_fig[action_row, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = apply_width)
    on(_ -> on_apply!(cbs), btn_apply.clicks)
    _display_in_screen(menu_fig)
    return menu_fig
end

function _sync_lines!(lines::Vector)
    length(lines) > 1 || return
    syncing = Ref(false)
    for line in lines
        other_lines = [l for l in lines if l !== line]
        on(line.visible) do visible_val
            syncing[] && return
            syncing[] = true
            for other_line in other_lines
                other_line.visible[] != visible_val && (other_line.visible = visible_val)
            end
            syncing[] = false
        end
    end
end

function _toggle_legends!(legend_refs)
    isnothing(legend_refs) && return
    for leg in legend_refs
        if !isnothing(leg) && isdefined(leg, :blockscene)
            leg.blockscene.visible[] = !leg.blockscene.visible[]
        end
    end
end

function _filter_visible_conditions(data, condition_checked_ref)
    if isnothing(condition_checked_ref[]) || !(data isa AbstractVector)
        return data
    end
    condition_checked = condition_checked_ref[]
    length(condition_checked) == length(data) || return data
    visible_data = [data[i] for i in eachindex(data) if condition_checked[i][]]
    return length(visible_data) == 1 ? visible_data[1] : visible_data
end

function _get_selection_bounds(selection_state, data)
    x_min, x_max = minmax(selection_state.bounds[]...)
    return (data, x_min, x_max)
end

function _show_waveform_topo_menu!(selection_state, data, condition_checked_ref; window_title::String = "")
    !isempty(window_title) && _set_window_title(window_title)
    menu_fig = _eegfun_figure(size = (300, 300))
    data_to_plot = _filter_visible_conditions(data, condition_checked_ref)
    plot_types =
        data_to_plot isa AbstractVector && length(data_to_plot) > 1 ?
        ["Topoplot", "Topoplot (average)", "Topoplot 3D", "Topoplot 3D (average)"] : ["Topoplot", "Topoplot 3D"]

    for (idx, plot_type) in enumerate(plot_types)
        btn = Button(menu_fig[idx, 1], label = plot_type)
        on(btn.clicks) do _
            original_data, x_min, x_max = _get_selection_bounds(selection_state, data)
            time_sel = x -> (x.time .>= x_min) .& (x.time .<= x_max)
            vis_data = _filter_visible_conditions(original_data, condition_checked_ref)
            erp_data = if vis_data isa EpochData
                average_epochs(vis_data)
            elseif vis_data isa AbstractVector{<:EpochData}
                erps = [average_epochs(d) for d in vis_data]
                length(erps) == 1 ? erps[1] : erps
            else
                vis_data
            end

            if btn.label[] == "Topoplot"
                plot_topography(erp_data, sample_selection = time_sel)
            elseif btn.label[] == "Topoplot (average)"
                plot_topography(erp_data isa Vector{ErpData} ? _average_conditions(erp_data) : erp_data, sample_selection = time_sel)
            elseif btn.label[] == "Topoplot 3D"
                plot_topography_3d(erp_data, sample_selection = time_sel)
            elseif btn.label[] == "Topoplot 3D (average)"
                plot_topography_3d(erp_data isa Vector{ErpData} ? _average_conditions(erp_data) : erp_data, sample_selection = time_sel)
            end
        end
    end

    _display_in_screen(menu_fig; size = (300, 300))
end

function _handle_waveform_right_click!(selection_state, mouse_x, data, condition_checked_ref; window_title::String = "")
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_waveform_topo_menu!(selection_state, data, condition_checked_ref; window_title = window_title)
    end
end

function _setup_waveform_control_panel!(
    on_update!::Function,
    fig::Figure,
    dat_subset::AbstractVector,
    baseline_interval::Interval,
    condition_checked_ref::Ref{Union{Vector{Observable{Bool}},Nothing}},
    legend_refs;
    panel_title::String = "Control Panel",
)
    control_fig = Ref{Union{Figure,Nothing}}(nothing)
    control_screen = Ref{Any}(nothing)

    start_val, stop_val = _extract_baseline_values(baseline_interval)
    baseline_start_obs = Observable(isnothing(start_val) ? "" : string(start_val))
    baseline_stop_obs = Observable(isnothing(stop_val) ? "" : string(stop_val))
    condition_checked = [Observable(true) for _ in dat_subset]
    condition_checked_ref[] = condition_checked
    previous_baseline = Ref{Union{Tuple{Float64,Float64},Nothing}}(nothing)

    function update_plot!()
        s_val, e_val = _parse_baseline_values(baseline_start_obs[], baseline_stop_obs[])
        baseline_interval_new = (!isnothing(s_val) && !isnothing(e_val)) ? (s_val, e_val) : nothing
        baseline_changed = baseline_interval_new != previous_baseline[] && !isnothing(baseline_interval_new)
        baseline_changed && (previous_baseline[] = baseline_interval_new)
        condition_mask = Bool[checked[] for checked in condition_checked]
        on_update!(baseline_changed, baseline_interval_new, condition_mask)
    end

    for checked in condition_checked
        on(_ -> update_plot!(), checked)
    end

    on(events(fig).keyboardbutton) do event
        event.action == Keyboard.press || return
        if event.key == Keyboard.c
            if !isnothing(control_screen[]) && hasmethod(isopen, Tuple{typeof(control_screen[])}) && isopen(control_screen[])
                return
            end

            control_fig[] = _eegfun_figure(title = panel_title, size = (320, 420))
            panel_layout = GridLayout(control_fig[][1, 1], tellwidth = false, rowgap = 10)

            Label(panel_layout[1, 1], "Baseline Interval", fontsize = 14, font = :bold)
            baseline_layout = GridLayout(panel_layout[2, 1], tellwidth = false, colgap = 10)
            _create_baseline_textbox(baseline_layout, 1, "Start (s):", baseline_start_obs, " ", 100)
            _create_baseline_textbox(baseline_layout, 2, "End (s):", baseline_stop_obs, " ", 100)

            apply_btn = Button(panel_layout[3, 1], label = "Apply Baseline", width = 200)
            on(_ -> update_plot!(), apply_btn.clicks)

            Label(panel_layout[4, 1], "Conditions", fontsize = 14, font = :bold)
            conditions_layout = GridLayout(panel_layout[5, 1], tellwidth = false, rowgap = 5)
            for (i, dat) in enumerate(dat_subset)
                cb = Checkbox(conditions_layout[i, 1], checked = condition_checked[i][])
                Label(conditions_layout[i, 2], dat.condition_name)
                connect!(condition_checked[i], cb.checked)
            end

            leg_btn = Button(panel_layout[6, 1], label = "Toggle Legend (L)", width = 200)
            on(_ -> _toggle_legends!(legend_refs), leg_btn.clicks)

            control_screen[] = _display_in_screen(control_fig[])
        elseif event.key == Keyboard.l
            _toggle_legends!(legend_refs)
        end
    end
end
