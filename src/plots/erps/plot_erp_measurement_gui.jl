"""
    plot_erp_measurement_gui(; filter_lp = 30.0, channel = :Cz, analysis_type = "mean_amplitude",
                               analysis_interval = (0.3, 0.5), baseline_interval = (-0.2, 0.0), kwargs...)
    plot_erp_measurement_gui(filepath::String; kwargs...)
    plot_erp_measurement_gui(erp::ErpData; kwargs...)
    plot_erp_measurement_gui(erps::AbstractVector{<:ErpData}; kwargs...)

Launch an interactive GUI for exploring ERP measurements with live visual feedback.

# Arguments
- `filepath::String`: Path to a `.jld2` file containing ERP data
- `erp::ErpData` or `erps::AbstractVector{<:ErpData}`: Single ERP or multiple conditions to overlay

# Keyword Arguments
- `filter_lp::Union{Real,Nothing}`: Low-pass filter cutoff in Hz applied when loading default demo data (default: `30.0`; `nothing` disables)
- `channel::Union{Symbol,Vector{Symbol},Nothing}`: Initial channel(s) to display (default: `:Cz` for demo, first channel otherwise)
- `analysis_type::String`: Initial measurement type (default: `"mean_amplitude"`)
- `analysis_interval::Union{Tuple{Real,Real},Nothing}`: Initial measurement interval
- `baseline_interval::Union{Tuple{Real,Real},Nothing}`: Initial baseline interval
- `display_plot::Bool`: Display the figure window (default: `true`)

# Examples
```julia
plot_erp_measurement_gui()
plot_erp_measurement_gui("participant_1_erps_final.jld2")
plot_erp_measurement_gui(erp, analysis_type = "max_peak_latency", analysis_interval = (0.3, 0.5))
```
"""
function plot_erp_measurement_gui(;
    filter_lp::Union{Real,Nothing} = 30.0,
    channel::Union{Symbol,Vector{Symbol},Nothing} = :Cz,
    analysis_type::String = "mean_amplitude",
    analysis_interval::Union{Tuple{Real,Real},Nothing} = (0.3, 0.5),
    baseline_interval::Union{Tuple{Real,Real},Nothing} = (-0.2, 0.0),
    kwargs...,
)
    data = read_data(example_path("data/julia/erps/example1_erps.jld2"))
    !isnothing(filter_lp) && (data = lowpass_filter(data, filter_lp))
    return plot_erp_measurement_gui(
        data;
        channel = channel,
        analysis_type = analysis_type,
        analysis_interval = analysis_interval,
        baseline_interval = baseline_interval,
        kwargs...,
    )
end

function plot_erp_measurement_gui(filepath::String; kwargs...)
    data = read_data(filepath)
    isnothing(data) && @minimal_error "No data found in file: $filepath"
    return plot_erp_measurement_gui(data; kwargs...)
end

function plot_erp_measurement_gui(
    erp::ErpData;
    channel::Union{Symbol,Vector{Symbol},Nothing} = nothing,
    analysis_type::String = "mean_amplitude",
    analysis_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    display_plot::Bool = true,
)
    return plot_erp_measurement_gui([erp]; channel, analysis_type, analysis_interval, baseline_interval, display_plot)
end

function plot_erp_measurement_gui(
    erp_vec::AbstractVector{<:ErpData};
    channel::Union{Symbol,Vector{Symbol},Nothing} = nothing,
    analysis_type::String = "mean_amplitude",
    analysis_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    baseline_interval::Union{Tuple{Real,Real},Nothing} = nothing,
    display_plot::Bool = true,
)
    first_erp = erp_vec[1]
    all_channels = setdiff(propertynames(first_erp.data), meta_labels(first_erp))
    isempty(all_channels) && @minimal_error "No channels found in data!"

    initial_channels = isnothing(channel) ? [all_channels[1]] : (channel isa Symbol ? [channel] : channel)
    for ch in initial_channels
        ch in all_channels || @minimal_error "Channel $ch not found in data. Available: $(all_channels)"
    end

    time_data = time_vector(first_erp)
    time_min, time_max = extrema(time_data)
    isnothing(analysis_interval) && (analysis_interval = (time_min, time_max))

    selected_channels = Observable{Vector{Symbol}}(initial_channels)
    selected_type = Observable(analysis_type)
    analysis_interval_obs = Observable(analysis_interval)
    baseline_interval_obs = Observable{Interval}(baseline_interval)
    show_markers_obs = Observable(false)

    fig = _eegfun_figure(size = (1280, 960))
    _set_window_title("ERP Measurement Tool")

    main_grid = fig[1, 1] = GridLayout()
    controls_grid = main_grid[1, 1] = GridLayout(valign = :top, tellheight = false)
    plot_area_grid = main_grid[1, 2] = GridLayout()
    colsize!(main_grid, 1, Relative(0.2))
    colsize!(main_grid, 2, Relative(0.8))
    rowsize!(main_grid, 1, Auto())

    Label(controls_grid[1, 1], "Measurement Controls", fontsize = 16, font = :bold, halign = :left)
    Label(controls_grid[2, 1], "Channel:", halign = :left)
    channel_menu = Menu(controls_grid[3, 1], options = zip(string.(all_channels), all_channels), default = string(initial_channels[1]))

    Label(controls_grid[4, 1], "Measurement Type:", halign = :left)
    type_idx = findfirst(p -> p.second == analysis_type, MEASUREMENT_TYPE_LABELS)
    default_type_label = isnothing(type_idx) ? "Mean Amplitude" : MEASUREMENT_TYPE_LABELS[type_idx].first
    type_menu = Menu(
        controls_grid[5, 1],
        options = zip([p.first for p in MEASUREMENT_TYPE_LABELS], MEASUREMENT_TYPE_LABELS),
        default = default_type_label,
    )

    Label(controls_grid[6, 1], "Measurement Interval:", halign = :left)
    meas_window_slider = IntervalSlider(controls_grid[7, 1], range = time_min:0.005:time_max, startvalues = analysis_interval)
    meas_window_label = Label(controls_grid[8, 1], @sprintf("%.3f s - %.3f s", analysis_interval...), halign = :left)

    Label(controls_grid[9, 1], "Baseline Interval:", halign = :left)
    slider_bl = isnothing(baseline_interval) ? (time_min, time_max) : baseline_interval
    baseline_interval_slider = IntervalSlider(controls_grid[10, 1], range = time_min:0.005:time_max, startvalues = slider_bl)
    baseline_interval_label = Label(
        controls_grid[11, 1],
        isnothing(baseline_interval) ? "None" : @sprintf("%.3f s - %.3f s", baseline_interval...),
        halign = :left,
    )

    peak_lbl1 = Label(controls_grid[12, 1], "Peak Local Window:", halign = :left)
    local_interval_slider = Slider(controls_grid[13, 1], range = 1:15, startvalue = 3)
    peak_lbl2 = Label(controls_grid[14, 1], lift(v -> "$v samples", local_interval_slider.value), halign = :left)

    frac_lbl1 = Label(controls_grid[15, 1], "Fraction:", halign = :left)
    fraction_slider = Slider(controls_grid[16, 1], range = 0.05:0.05:0.95, startvalue = 0.5)
    frac_lbl2 = Label(controls_grid[17, 1], lift(v -> @sprintf("%.2f", v), fraction_slider.value), halign = :left)

    function update_param_sliders_visibility!(atype)
        is_peak = atype in (
            "max_peak_amplitude",
            "min_peak_amplitude",
            "max_peak_latency",
            "min_peak_latency",
            "peak_to_peak_amplitude",
            "peak_to_peak_latency",
            "fractional_peak_latency",
        )
        is_frac = atype in ("fractional_area_latency", "fractional_peak_latency")

        peak_lbl1.blockscene.visible[] = is_peak
        local_interval_slider.blockscene.visible[] = is_peak
        peak_lbl2.blockscene.visible[] = is_peak
        frac_lbl1.blockscene.visible[] = is_frac
        fraction_slider.blockscene.visible[] = is_frac
        frac_lbl2.blockscene.visible[] = is_frac
    end
    update_param_sliders_visibility!(analysis_type)

    Label(controls_grid[18, 1], "Show Result Markers:", halign = :left)
    show_markers_toggle = Toggle(controls_grid[19, 1], active = false)

    Label(controls_grid[20, 1], "Result:", fontsize = 14, font = :bold, halign = :right, justification = :right)
    result_label = Label(controls_grid[21, 1], "—", fontsize = 16, halign = :right, justification = :right)

    ax_topo = Axis(controls_grid[22, 1], aspect = DataAspect())
    plot_layout_2d!(fig, ax_topo, first_erp.layout, point_plot = true, label_plot = false, head_linewidth = 1)
    hidedecorations!(ax_topo)
    hidespines!(ax_topo)
    deregister_interaction!(ax_topo, :rectanglezoom)
    deregister_interaction!(ax_topo, :scrollzoom)

    ch_selection_active = Ref(false)
    ch_selection_start = Ref(Point2f(0, 0))
    ch_selection_rect = poly!(ax_topo, Point2f[(0, 0), (0, 0), (0, 0), (0, 0)], color = (:blue, 0.3), strokecolor = :blue, visible = false)

    on(events(ax_topo).mousebutton, priority = 50) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            ax_pos = ax_topo.scene.viewport[]
            mouse_pos = events(fig).mouseposition[]
            if !(
                ax_pos.origin[1] <= mouse_pos[1] <= ax_pos.origin[1] + ax_pos.widths[1] &&
                ax_pos.origin[2] <= mouse_pos[2] <= ax_pos.origin[2] + ax_pos.widths[2]
            )
                return Consume(false)
            end

            ch_selection_active[] = true
            pt0 = mouseposition(ax_topo)
            ch_selection_start[] = pt0
            ch_selection_rect[1] = Point2f[pt0, pt0, pt0, pt0]
            ch_selection_rect.visible[] = true
            return Consume(true)
        elseif event.action == Mouse.release && ch_selection_active[]
            ch_selection_active[] = false
            current_pos = mouseposition(ax_topo)
            ch_selection_rect.visible[] = false

            x1, x2 = minmax(ch_selection_start[][1], current_pos[1])
            y1, y2 = minmax(ch_selection_start[][2], current_pos[2])

            layout_df = first_erp.layout.data
            lx = layout_df[!, :x2]
            ly = layout_df[!, :y2]
            ll = layout_df[!, :label]

            selected_this_time = Symbol[]
            for i in eachindex(ll)
                if x1 <= lx[i] <= x2 && y1 <= ly[i] <= y2
                    ch_sym = Symbol(ll[i])
                    ch_sym in all_channels && push!(selected_this_time, ch_sym)
                end
            end

            if isempty(selected_this_time) && abs(x2 - x1) < 0.05 && abs(y2 - y1) < 0.05
                min_dist = Inf
                min_idx = 0
                for i in eachindex(ll)
                    d2 = (lx[i] - current_pos[1])^2 + (ly[i] - current_pos[2])^2
                    if d2 < min_dist
                        min_dist = d2
                        min_idx = i
                    end
                end
                if min_idx > 0 && min_dist < 0.04
                    ch_sym = Symbol(ll[min_idx])
                    ch_sym in all_channels && push!(selected_this_time, ch_sym)
                end
            end

            if !isempty(selected_this_time)
                ks = events(fig).keyboardstate
                mod_held =
                    Makie.Keyboard.left_shift in ks ||
                    Makie.Keyboard.right_shift in ks ||
                    Makie.Keyboard.left_control in ks ||
                    Makie.Keyboard.right_control in ks

                if mod_held
                    new_sel = copy(selected_channels[])
                    for c in selected_this_time
                        c in new_sel ? filter!(!=(c), new_sel) : push!(new_sel, c)
                    end
                    selected_channels[] = isempty(new_sel) ? [all_channels[1]] : new_sel
                else
                    selected_channels[] = selected_this_time
                end

                if length(selected_channels[]) == 1
                    idx = findfirst(==(selected_channels[][1]), all_channels)
                    channel_menu.i_selected[] = isnothing(idx) ? 0 : idx
                else
                    channel_menu.i_selected[] = 0
                end
            end
            return Consume(true)
        end
        return Consume(false)
    end

    on(events(ax_topo).mouseposition) do _
        if ch_selection_active[]
            cp = mouseposition(ax_topo)
            sp = ch_selection_start[]
            ch_selection_rect[1] = Point2f[(sp[1], sp[2]), (cp[1], sp[2]), (cp[1], cp[2]), (sp[1], cp[2])]
        end
    end

    active_pos = lift(selected_channels) do chs
        ll = first_erp.layout.data[!, :label]
        lx = first_erp.layout.data[!, :x2]
        ly = first_erp.layout.data[!, :y2]
        pts = Point2f[]
        for ch in chs
            idx = findfirst(==(Symbol(ch)), ll)
            !isnothing(idx) && push!(pts, Point2f(lx[idx], ly[idx]))
        end
        return pts
    end
    scatter!(ax_topo, active_pos, color = :red, markersize = 14)
    rowgap!(controls_grid, 10)

    ax = Axis(
        plot_area_grid[1, 1],
        xlabel = "Time (s)",
        ylabel = "μV",
        title = "$(_print_vector(initial_channels)): $(analysis_type)",
        xgridvisible = false,
        ygridvisible = false,
        xminorgridvisible = false,
        yminorgridvisible = false,
    )

    analysis_band_plots = AbstractPlot[]
    baseline_band_plots = AbstractPlot[]

    function _draw_interval_band!(plots_list, interval, band_color, edge_color)
        for p in plots_list
            delete!(ax, p)
        end
        empty!(plots_list)
        isnothing(interval) && return
        bh = 0.01 * ax.finallimits[].widths[2]
        p_band = poly!(ax, Point2f[(interval[1], -bh), (interval[2], -bh), (interval[2], bh), (interval[1], bh)], color = (band_color, 0.3))
        p_left = vlines!(ax, interval[1], color = (edge_color, 0.2), linewidth = 1, linestyle = :dash)
        p_right = vlines!(ax, interval[2], color = (edge_color, 0.2), linewidth = 1, linestyle = :dash)
        push!(plots_list, p_band, p_left, p_right)
    end

    update_analysis_band!(mw) = _draw_interval_band!(analysis_band_plots, mw, :blue, :blue)
    update_baseline_band!(bw) = _draw_interval_band!(baseline_band_plots, bw, :gray, :gray)

    marker_plots = AbstractPlot[]
    function update_result_markers!(results, show)
        for p in marker_plots
            delete!(ax, p)
        end
        empty!(marker_plots)
        (!show || isempty(results)) && return

        colors = length(erp_vec) > 1 ? Makie.cgrad(DEFAULT_COLORMAP, length(erp_vec), categorical = true) : [:black]
        for (idx, result) in enumerate(results)
            (haskey(result, :error) || !haskey(result, :value) || isnothing(result.value) || isnan(result.value)) && continue
            c = (colors[idx], 0.8)

            p = if haskey(result, :point)
                _is_latency_measurement(selected_type[]) ? vlines!(ax, result.point[1]; color = c, linewidth = 2) :
                hlines!(ax, result.point[2]; color = c, linewidth = 2)
            elseif haskey(result, :points)
                lines!(ax, [Point2f(pt) for pt in result.points]; color = c, linestyle = :dash, linewidth = 2)
            elseif haskey(result, :span)
                lines!(ax, Point2f[(result.span[1], result.value), (result.span[2], result.value)]; color = c, linewidth = 4)
            elseif haskey(result, :curve)
                y_band = result.curve
                y0 = zeros(Float64, length(y_band))
                y_lower, y_upper = if selected_type[] == "positive_area"
                    (y0, max.(y_band, 0.0))
                elseif selected_type[] == "negative_area"
                    (min.(y_band, 0.0), y0)
                else
                    (y0, y_band)
                end
                band!(ax, result.x, y_lower, y_upper; color = (colors[idx], 0.3))
            else
                _is_latency_measurement(selected_type[]) ? vlines!(ax, result.value, color = c, linestyle = :dot, linewidth = 2) :
                hlines!(ax, result.value, color = c, linestyle = :dot, linewidth = 2)
            end
            push!(marker_plots, p)
        end
    end

    cached_ylims = Ref{Tuple{Float64,Float64}}((0.0, 1.0))
    cached_ch = Ref{Vector{Symbol}}(Symbol[])
    current_legend = Ref{Union{Legend,Nothing}}(nothing)

    function update_erp_plot!()
        empty!(ax)
        if !isnothing(current_legend[])
            delete!(current_legend[])
            current_legend[] = nothing
        end

        chs = selected_channels[]
        bi = baseline_interval_obs[]

        if chs != cached_ch[]
            max_range = 0.0
            for erp in erp_vec
                raw_data = length(chs) == 1 ? erp.data[!, chs[1]] : _colmeans(erp.data, chs)
                r_min, r_max = extrema(raw_data)
                max_range = max(max_range, r_max - r_min)
            end
            max_range == 0.0 && (max_range = 1.0)
            cached_ylims[] = (-1.05 * max_range, 1.05 * max_range)
            cached_ch[] = chs
        end

        if length(chs) > 1
            lbl_sym = Symbol("ROI_$(length(chs))_channels")
            plot_vec = channel_average(erp_vec, channel_selections = [channels(chs)], output_labels = [lbl_sym], reduce = true)
            plot_ch_selection = channels([lbl_sym])
        else
            plot_vec = erp_vec
            plot_ch_selection = channels(chs)
        end

        plot_erp!(fig, ax, plot_vec, channel_selection = plot_ch_selection, baseline_interval = bi, legend = false, legend_position = :rt)

        ylims!(ax, cached_ylims[][1], cached_ylims[][2])
        update_analysis_band!(analysis_interval_obs[])
        update_baseline_band!(baseline_interval_obs[])
        _set_window_title("ERP Measurement Tool")

        if length(erp_vec) > 1
            line_elements = filter(p -> p isa Lines, ax.scene.plots)
            labels = [erp.condition_name for erp in erp_vec]
            current_legend[] = Legend(plot_area_grid[1, 2], line_elements[1:length(erp_vec)], labels, halign = :right, valign = :top)
        end
    end

    update_erp_plot!()

    drag_target = Observable{Symbol}(:none)
    drag_suppress_slider = Ref(false)
    drag_offset = Ref(0.0)

    function _find_drag_target(mouse_x)
        tolerance = 0.02 * (time_max - time_min)
        ai = analysis_interval_obs[]
        bi = baseline_interval_obs[]

        abs(mouse_x - ai[1]) < tolerance && return :analysis_left
        abs(mouse_x - ai[2]) < tolerance && return :analysis_right
        if !isnothing(bi)
            abs(mouse_x - bi[1]) < tolerance && return :baseline_left
            abs(mouse_x - bi[2]) < tolerance && return :baseline_right
        end
        ai[1] < mouse_x < ai[2] && return :analysis_middle
        !isnothing(bi) && bi[1] < mouse_x < bi[2] && return :baseline_middle
        return :none
    end

    is_recalculating = Ref(false)
    latest_results = Ref{Vector{NamedTuple}}(NamedTuple[])

    function recalculate_and_update!()
        is_recalculating[] && return
        is_recalculating[] = true
        try
            chs = selected_channels[]
            type_str = selected_type[]
            mw = analysis_interval_obs[]
            bw = baseline_interval_obs[]
            frac = fraction_slider.value[]

            measurement_kwargs = Dict{Symbol,Any}(k => v[1] for (k, v) in ERP_MEASUREMENTS_KWARGS)
            measurement_kwargs[:local_interval] = local_interval_slider.value[]
            measurement_kwargs[:fractional_area_fraction] = frac
            measurement_kwargs[:fractional_peak_fraction] = frac

            results = NamedTuple[]
            for erp in erp_vec
                try
                    res = _compute_gui_measurement(erp, chs, type_str, mw, bw, measurement_kwargs)
                    push!(results, (condition = erp.condition_name, res...))
                catch e
                    push!(results, (condition = erp.condition_name, value = NaN, error = string(e)))
                end
            end
            latest_results[] = results

            if isempty(results)
                result_label.text = "—"
            else
                fmt_val(val) =
                    if isnothing(val) || isnan(val)
                        "   N/A  "
                    elseif _is_latency_measurement(type_str)
                        @sprintf("% .3f s", val)
                    elseif type_str in ("rectified_area", "integral", "positive_area", "negative_area")
                        @sprintf("% .3f μV·s", val)
                    else
                        @sprintf("% .3f μV", val)
                    end

                if length(results) == 1
                    r = results[1]
                    result_label.text = haskey(r, :error) ? "Error: $(r.error)" : fmt_val(r.value)
                else
                    result_label.text = join(
                        [haskey(r, :error) ? "$(r.condition): $(r.error)" : "$(r.condition): $(fmt_val(r.value))" for r in results],
                        "\n",
                    )
                end
            end

            update_result_markers!(results, show_markers_obs[])
        finally
            is_recalculating[] = false
        end
    end

    on(events(ax).mousebutton, priority = 100) do event
        event.button == Mouse.left || return Consume(false)
        if event.action == Mouse.press
            mouse_pos = events(fig).mouseposition[]
            area = ax.scene.viewport[]
            if !(
                area.origin[1] <= mouse_pos[1] <= area.origin[1] + area.widths[1] &&
                area.origin[2] <= mouse_pos[2] <= area.origin[2] + area.widths[2]
            )
                return Consume(false)
            end

            world_x = mouseposition(ax)[1]
            target = _find_drag_target(world_x)
            if target != :none
                drag_target[] = target
                if target == :analysis_middle
                    drag_offset[] = world_x - analysis_interval_obs[][1]
                elseif target == :baseline_middle
                    drag_offset[] = world_x - baseline_interval_obs[][1]
                end
                return Consume(true)
            end
        elseif event.action == Mouse.release && drag_target[] != :none
            released_target = drag_target[]
            drag_target[] = :none
            released_target in (:baseline_left, :baseline_right, :baseline_middle) && update_erp_plot!()
            recalculate_and_update!()
            return Consume(true)
        end
        return Consume(false)
    end

    function _apply_interval_drag!(obs, slider, lbl, side::Symbol, world_x::Real)
        iv = obs[]
        new_left, new_right = if side === :left
            (min(world_x, iv[2] - 0.005), iv[2])
        elseif side === :right
            (iv[1], max(world_x, iv[1] + 0.005))
        else
            w = iv[2] - iv[1]
            nl = clamp(world_x - drag_offset[], time_min, time_max - w)
            (nl, nl + w)
        end
        obs[] = (new_left, new_right)
        set_close_to!(slider, new_left, new_right)
        lbl.text = @sprintf("%.3f s - %.3f s", new_left, new_right)
    end

    on(events(fig).mouseposition) do _
        target = drag_target[]
        target == :none && return

        world_x = clamp(mouseposition(ax)[1], time_min, time_max)
        drag_suppress_slider[] = true
        if target in (:analysis_left, :analysis_right, :analysis_middle)
            side = target === :analysis_left ? :left : (target === :analysis_right ? :right : :middle)
            _apply_interval_drag!(analysis_interval_obs, meas_window_slider, meas_window_label, side, world_x)
        else
            side = target === :baseline_left ? :left : (target === :baseline_right ? :right : :middle)
            _apply_interval_drag!(baseline_interval_obs, baseline_interval_slider, baseline_interval_label, side, world_x)
            update_erp_plot!()
        end
        recalculate_and_update!()
        drag_suppress_slider[] = false
    end

    on(channel_menu.selection) do ch
        if !isnothing(ch)
            sym_ch = Symbol(ch)
            if length(selected_channels[]) > 1 || selected_channels[][1] != sym_ch
                selected_channels[] = [sym_ch]
            end
        end
    end

    on(selected_channels) do chs
        update_erp_plot!()
        ax.title = "$(_print_vector(chs)): $(selected_type[])"
        recalculate_and_update!()
    end

    on(type_menu.selection) do type_pair
        type_str = type_pair isa Pair ? type_pair[2] : type_pair
        selected_type[] = type_str
        ax.title = "$(_print_vector(selected_channels[])): $(type_str)"
        update_param_sliders_visibility!(type_str)
        recalculate_and_update!()
    end

    on(meas_window_slider.interval) do interval
        drag_suppress_slider[] && return
        analysis_interval_obs[] = (interval[1], interval[2])
        meas_window_label.text = @sprintf("%.3f s - %.3f s", interval[1], interval[2])
        recalculate_and_update!()
    end

    on(baseline_interval_slider.interval) do interval
        drag_suppress_slider[] && return
        baseline_interval_obs[] = (interval[1], interval[2])
        baseline_interval_label.text = @sprintf("%.3f s - %.3f s", interval[1], interval[2])
        update_erp_plot!()
        recalculate_and_update!()
    end

    on(_ -> recalculate_and_update!(), local_interval_slider.value)
    on(_ -> recalculate_and_update!(), fraction_slider.value)
    on(update_analysis_band!, analysis_interval_obs)
    on(update_baseline_band!, baseline_interval_obs)

    on(show_markers_toggle.active) do enabled
        show_markers_obs[] = enabled
        update_result_markers!(latest_results[], enabled)
    end

    recalculate_and_update!()

    screen = display_plot ? _display_figure(fig) : nothing
    return (fig = fig, axes = [ax], screen = screen)
end

_is_latency_measurement(atype::AbstractString) = endswith(atype, "latency")

function _compute_gui_measurement(
    erp::ErpData,
    channels::Vector{Symbol},
    analysis_type,
    analysis_interval,
    baseline_interval,
    measurement_kwargs,
)
    erp_data = isnothing(baseline_interval) ? erp : baseline(erp, baseline_interval)
    time_data = time_vector(erp_data)
    channel_data = length(channels) == 1 ? erp_data.data[!, channels[1]] : _colmeans(erp_data.data, channels)

    sliced = _select_interval_slice(channel_data, time_data, analysis_interval[1], analysis_interval[2])
    isnothing(sliced) && return (value = NaN,)
    selected_data, selected_times = sliced

    dummy_channel_sym = length(channels) == 1 ? channels[1] : Symbol("ROI_$(length(channels))_channels")
    value = _compute_measurement(selected_data, selected_times, analysis_type, measurement_kwargs, dummy_channel_sym)
    isnothing(value) && return (value = NaN,)

    if analysis_type == "mean_amplitude"
        return (value = value, span = (analysis_interval[1], analysis_interval[2]))
    elseif analysis_type in ("max_peak_amplitude", "min_peak_amplitude", "max_peak_latency", "min_peak_latency")
        peak_type = startswith(analysis_type, "max") ? :max : :min
        peak_val, peak_idx = _compute_peak_measurement(selected_data, peak_type, measurement_kwargs[:local_interval], dummy_channel_sym)
        return (!isnothing(peak_val) && !isnothing(peak_idx)) ? (value = value, point = (selected_times[peak_idx], peak_val)) :
               (value = value,)
    elseif analysis_type in ("peak_to_peak_amplitude", "peak_to_peak_latency")
        li = measurement_kwargs[:local_interval]
        max_val, max_idx = _compute_peak_measurement(selected_data, :max, li, dummy_channel_sym)
        min_val, min_idx = _compute_peak_measurement(selected_data, :min, li, dummy_channel_sym)
        return (!isnothing(max_val) && !isnothing(min_val)) ?
               (value = value, points = [(selected_times[max_idx], max_val), (selected_times[min_idx], min_val)]) : (value = value,)
    elseif analysis_type in ("rectified_area", "integral", "positive_area", "negative_area")
        return (value = value, curve = collect(selected_data), x = collect(selected_times))
    elseif analysis_type in ("fractional_area_latency", "fractional_peak_latency")
        idx = argmin(abs.(time_data .- value))
        return (value = value, point = (value, channel_data[idx]))
    end

    return (value = value,)
end
