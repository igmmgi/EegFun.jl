# Interactive ICA component scalp topographies and scrolling activation browser.

######################################
######### Types & Constants ##########
######################################

const _CHANNEL_OVERLAY_COLORS = [:grey, :dodgerblue, :darkorange, :mediumpurple, :forestgreen, :crimson, :goldenrod, :deeppink]

mutable struct IcaOverlayState
    names::Vector{Symbol}
    plots::Vector{Vector{Any}}
end

mutable struct IcaArtifactState
    data::Union{Nothing,ArtifactComponents}
    toggles::Union{Nothing,Dict{Symbol,Observable{Bool}}}
    textbox::Union{Nothing,Textbox}
end

mutable struct IcaComponentState
    dat::ContinuousData
    ica::InfoIca
    component_data::Matrix{Float64}
    n_visible_components::Int
    window_size::Int
    comp_start::Observable{Int}
    xrange::Observable{UnitRange{Int}}
    ylims::Observable{Tuple{Float64,Float64}}
    channel_yscale::Observable{Float64}
    use_global_scale::Observable{Bool}
    invert_scale::Observable{Bool}
    components::Vector{Int}
    plot_kwargs::Dict{Symbol,Any}
    fig::Union{Nothing,Figure}
    axs::Vector{Axis}
    channel_axs::Vector{Union{Axis,Nothing}}
    topo_axs::Vector{Axis}
    lines_obs::Vector{Observable{Vector{Float64}}}
    overlay::IcaOverlayState
    artifacts::IcaArtifactState

    function IcaComponentState(
        dat,
        ica,
        component_selection,
        n_visible_components,
        window_size;
        method = :thin_plate,
        gridscale = 100,
        kwargs...,
    )
        dat_matrix = _prepare_ica_plot_matrix(dat, ica)
        component_data = ica.unmixing * dat_matrix
        total_components, n_samps = size(component_data)

        comp_start = Observable(1)
        use_global_scale = Observable(Bool(get(kwargs, :use_global_scale, false)))
        invert_scale = Observable(false)

        time_zero_idx = clamp(searchsortedfirst(dat.data.time, 0.0), 1, length(dat.data.time))
        start_idx = max(1, time_zero_idx - div(window_size, 2))
        end_idx = min(n_samps, start_idx + window_size - 1)
        end_idx == n_samps && (start_idx = max(1, end_idx - window_size + 1))
        xrange = Observable(start_idx:end_idx)

        comps_to_use = findall(component_selection(1:total_components))

        max_abs = 0.0
        has_val = false
        if !isempty(comps_to_use) && all(<=(total_components), comps_to_use)
            @inbounds for s = start_idx:end_idx, c in comps_to_use
                v = abs(component_data[c, s])
                v > max_abs && (max_abs = v)
                has_val = true
            end
        end
        initial_range = (has_val && max_abs > 0.0) ? max_abs : 1.0

        ylims = Observable((-initial_range, initial_range))
        channel_yscale = Observable(1.0)

        plot_kwargs = Dict{Symbol,Any}(pairs(kwargs))
        delete!(plot_kwargs, :use_global_scale)
        plot_kwargs[:method] = method
        plot_kwargs[:gridscale] = gridscale

        new(
            dat,
            ica,
            component_data,
            n_visible_components,
            window_size,
            comp_start,
            xrange,
            ylims,
            channel_yscale,
            use_global_scale,
            invert_scale,
            comps_to_use,
            plot_kwargs,
            nothing,
            Axis[],
            Union{Axis,Nothing}[],
            Axis[],
            Observable{Vector{Float64}}[],
            IcaOverlayState(Symbol[], Vector{Any}[]),
            IcaArtifactState(nothing, nothing, nothing),
        )
    end
end

######################################
########### User Functions ###########
######################################

"""
    plot_topography(ica::InfoIca; component_selection = components(), kwargs...) -> NamedTuple

Plot ICA component scalp topographies in a grid layout.

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))
"""
function plot_topography(ica::InfoIca; component_selection = components(), kwargs...)
    merged_kwargs = Dict{Symbol,Any}(:colorbar_label => "a.u.", pairs(kwargs)...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, merged_kwargs)

    dims = plot_kwargs[:dims]
    display_plot = plot_kwargs[:display_plot]
    show_title = plot_kwargs[:show_title]
    colorbar_plot = plot_kwargs[:colorbar_plot]

    ica = _ica_with_plot_coords(ica)

    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    comps = get_selected_components(ica, component_selection)
    if isempty(comps)
        return (fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...), axes = Axis[])
    end

    fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    isnothing(dims) && (dims = _best_rect(length(comps)))
    (length(dims) == 2 && all(dims .> 0)) || @minimal_error "Invalid dimensions: $dims. Expected [rows, cols] with positive values."
    total_cells = dims[1] * dims[2]
    total_cells < length(comps) && @minimal_error "Grid dimensions $dims provide $total_cells cells but need $(length(comps))."

    all_data, x_bounds, y_bounds = _prepare_ica_topo_data(ica, comps, plot_kwargs[:method], plot_kwargs[:gridscale])
    use_global = Bool(plot_kwargs[:use_global_scale])
    levels_result =
        isnothing(plot_kwargs[:colorrange]) ? _calculate_ica_topo_levels(all_data, use_global, plot_kwargs[:num_levels]) : nothing

    axes = Vector{Axis}(undef, length(comps))
    for i in eachindex(comps)
        base_row, base_col = divrem(i - 1, dims[2]) .+ (1, 1)
        default_title = @sprintf("IC %d (%.1f%%)", comps[i], ica.variance[comps[i]] * 100)
        ax_title = if !show_title
            ""
        elseif !isnothing(plot_kwargs[:plot_title]) && length(comps) == 1
            String(plot_kwargs[:plot_title])
        else
            default_title
        end

        ax = if colorbar_plot
            sg = GridLayout(fig[base_row, base_col])
            Axis(sg[1, 1], aspect = DataAspect(), title = ax_title)
        else
            Axis(fig[base_row, base_col], aspect = DataAspect(), title = ax_title)
        end
        deregister_interaction!(ax, :rectanglezoom)
        deregister_interaction!(ax, :scrollzoom)
        deregister_interaction!(ax, :dragpan)
        axes[i] = ax

        comp_levels = isnothing(levels_result) ? nothing : (use_global ? levels_result : levels_result[i])
        _plot_ica_topography_impl!(
            fig,
            ax,
            ica,
            comps[i],
            plot_kwargs,
            colorbar_plot;
            precomputed_topo = (all_data[i], x_bounds, y_bounds),
            pre_calculated_levels = comp_levels,
        )
    end

    zoom_step = plot_kwargs[:zoom_step]
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key in (Keyboard.up, Keyboard.down)
            fig_axes = filter(a -> a isa Axis, fig.content)
            if event.key == Keyboard.up
                _topo_scale_up!.(fig_axes, zoom_step)
            else
                _topo_scale_down!.(fig_axes, zoom_step)
            end
        end
    end

    display_plot && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

"""
    plot_ica_component_activation(dat::ContinuousData, ica::InfoIca; artifacts = nothing, kwargs...) -> NamedTuple

Interactive viewer for ICA component time-course activations and scalp topographies.

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))
"""
function plot_ica_component_activation(dat::ContinuousData, ica::InfoIca; artifacts::Union{Nothing,ArtifactComponents} = nothing, kwargs...)
    _set_window_title(_generate_window_title(dat))

    merged_kwargs = Dict{Symbol,Any}(:colorbar_label => "a.u.", pairs(kwargs)...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, merged_kwargs)

    ica = _ica_with_plot_coords(ica)

    component_selection = pop!(plot_kwargs, :component_selection)
    selected_components = findall(component_selection(1:length(ica.ica_label)))
    isempty(selected_components) && @minimal_error("No ICA components matched component_selection")
    n_visible_components = min(length(selected_components), 10)
    window_size = min(2000, n_samples(dat))
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    display_plot = pop!(plot_kwargs, :display_plot)

    plot_kwargs[:point_plot] = false
    plot_kwargs[:label_plot] = false

    state = IcaComponentState(
        dat,
        ica,
        component_selection,
        n_visible_components,
        window_size;
        method = method,
        gridscale = gridscale,
        plot_kwargs...,
    )

    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    _create_component_activation_plots!(fig, state)
    _add_navigation_controls!(fig, state)
    _add_navigation_sliders!(fig, state)
    _add_channel_menu!(fig, state)
    !isnothing(artifacts) && _add_artifact_category_checkboxes!(fig, state, artifacts)
    _setup_keyboard_interactions!(fig, state)

    colsize!(fig.layout, 1, Relative(0.15))
    colsize!(fig.layout, 2, Relative(0.85))
    rowgap!(fig.layout, 5)
    rowgap!(fig.layout, state.n_visible_components, 80)

    screen = display_plot ? _display_figure(fig) : nothing
    _set_window_title("Makie")
    return (fig = fig, axes = state.axs, screen = screen)
end

######################################
######### Internal Functions #########
######################################

function _ica_level_range(min_val::Float64, max_val::Float64, num_levels::Int)
    num_levels <= 0 && @minimal_error("num_levels must be positive, got $num_levels")
    if min_val == max_val
        min_val -= 0.1
        max_val += 0.1
    end
    return range(min_val, max_val, length = num_levels)
end

function _calculate_topo_levels(
    data::AbstractMatrix{<:Real};
    use_global_scale = false,
    global_min = nothing,
    global_max = nothing,
    num_levels = 20,
)
    lo, hi =
        (use_global_scale && !isnothing(global_min) && !isnothing(global_max)) ? (Float64(global_min), Float64(global_max)) :
        _symmetric_colorrange(data)
    return _ica_level_range(lo, hi, num_levels)
end

function _calculate_ica_topo_levels(
    all_data::Vector{Matrix{Float64}},
    use_global_scale::Bool,
    num_levels::Int;
    global_min = nothing,
    global_max = nothing,
)
    if use_global_scale
        any(mat -> any(!isnan, mat), all_data) || @minimal_error "No valid (non-NaN) data. Cannot calculate topographic levels."
        lo, hi = if !isnothing(global_min) && !isnothing(global_max)
            Float64(global_min), Float64(global_max)
        else
            max_abs = maximum(mat -> _symmetric_colorrange(mat)[2], all_data)
            (-max_abs, max_abs)
        end
        return _ica_level_range(lo, hi, num_levels)
    end
    return [_calculate_topo_levels(mat; num_levels = num_levels) for mat in all_data]
end

function _resolve_ica_levels(data::AbstractMatrix{<:Real}, num_levels::Int; pre_calculated = nothing, colorrange = nothing)
    !isnothing(pre_calculated) && return pre_calculated
    !isnothing(colorrange) && return _calculate_topo_levels(
        data;
        use_global_scale = true,
        global_min = colorrange[1],
        global_max = colorrange[2],
        num_levels = num_levels,
    )
    return _calculate_topo_levels(data; num_levels = num_levels)
end

function _prepare_ica_topo_data(ica::InfoIca, comp_idx::Int, method::Symbol, gridscale::Int)
    return _topo_interpolate(@view(ica.mixing[:, comp_idx]), ica.layout, gridscale; method = method)
end

function _prepare_ica_topo_data(ica::InfoIca, comp_indices::AbstractVector{Int}, method::Symbol, gridscale::Int)
    isempty(comp_indices) && return Vector{Matrix{Float64}}(), (-1.0, 1.0), (-1.0, 1.0)
    if length(comp_indices) == 1
        data, x_bounds, y_bounds = _prepare_ica_topo_data(ica, comp_indices[1], method, gridscale)
        return [data], x_bounds, y_bounds
    end
    return _topo_interpolate(ica.mixing[:, comp_indices], ica.layout, gridscale; method = method)
end

function _plot_ica_topography_impl!(
    fig::Figure,
    ax::Axis,
    ica::InfoIca,
    component::Int,
    plot_kwargs::Dict{Symbol,Any},
    colorbar_plot::Bool;
    precomputed_topo = nothing,
    pre_calculated_levels = nothing,
)
    data, x_bounds, y_bounds =
        isnothing(precomputed_topo) ? _prepare_ica_topo_data(ica, component, plot_kwargs[:method], plot_kwargs[:gridscale]) :
        precomputed_topo
    levels =
        _resolve_ica_levels(data, plot_kwargs[:num_levels]; pre_calculated = pre_calculated_levels, colorrange = plot_kwargs[:colorrange])

    co, _, _ = _render_topo_surface!(
        fig,
        ax,
        Float64[],
        ica.layout;
        plot_kwargs...,
        data_interp = data,
        x_bounds = x_bounds,
        y_bounds = y_bounds,
        levels = levels,
    )

    colorbar_plot_numbers = plot_kwargs[:colorbar_plot_numbers]
    if colorbar_plot && (isempty(colorbar_plot_numbers) || component in colorbar_plot_numbers)
        cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
        sg = ax.layoutobservables.gridcontent[].parent
        Colorbar(sg[cb_pos...], co; plot_kwargs[:colorbar_kwargs]..., tellwidth = true, tellheight = false)
    end
    return co
end

function _plot_topography!(fig::Figure, ax::Axis, ica::InfoIca, component::Int; precomputed_topo = nothing, kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    ica = _ica_with_plot_coords(ica)
    return _plot_ica_topography_impl!(
        fig,
        ax,
        ica,
        component,
        plot_kwargs,
        plot_kwargs[:colorbar_plot];
        precomputed_topo = precomputed_topo,
    )
end

function _fill_ica_plot_row!(dat_matrix::Matrix{Float64}, c::Int, col::AbstractVector{<:Real}, mc::Float64, inv_scale::Float64)
    @inbounds for s in eachindex(col)
        dat_matrix[c, s] = (Float64(col[s]) - mc) * inv_scale
    end
    return nothing
end

function _prepare_ica_plot_matrix(dat::ContinuousData, ica::InfoIca)::Matrix{Float64}
    labels = ica.layout.data.label
    n_ch = length(labels)
    n_samps = nrow(dat.data)
    dat_matrix = Matrix{Float64}(undef, n_ch, n_samps)
    inv_scale = 1.0 / ica.scale
    @inbounds for c = 1:n_ch
        _fill_ica_plot_row!(dat_matrix, c, dat.data[!, labels[c]], ica.mean[c], inv_scale)
    end
    return dat_matrix
end

function _plot_ica_topo_in_viewer!(
    fig::Figure,
    topo_ax::Axis,
    ica::InfoIca,
    comp_idx::Int;
    gridscale = 100,
    num_levels = 20,
    method = :spherical_spline,
    pre_calculated_levels = nothing,
    precomputed_topo = nothing,
    kwargs...,
)
    topo_data, x_bounds, y_bounds =
        isnothing(precomputed_topo) ? _prepare_ica_topo_data(ica, comp_idx, method, gridscale) : precomputed_topo
    levels =
        _resolve_ica_levels(topo_data, num_levels; pre_calculated = pre_calculated_levels, colorrange = get(kwargs, :colorrange, nothing))
    co, _, _ = _render_topo_surface!(
        fig,
        topo_ax,
        Float64[],
        ica.layout;
        kwargs...,
        gridscale = gridscale,
        data_interp = topo_data,
        x_bounds = x_bounds,
        y_bounds = y_bounds,
        levels = levels,
    )
    return co
end

function _ensure_channel_overlay_axis!(state::IcaComponentState, i::Int)::Axis
    existing = state.channel_axs[i]
    !isnothing(existing) && return existing

    fig = state.fig
    ax_channel = Axis(
        fig[i, 2],
        yaxisposition = :right,
        xaxisposition = :top,
        yticklabelsvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        xticksvisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xminorgridvisible = false,
        yminorgridvisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        leftspinevisible = false,
    )
    linkxaxes!(state.axs[i], ax_channel)
    xlims!(ax_channel, (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])]))
    state.channel_axs[i] = ax_channel
    return ax_channel
end

function _create_component_activation_plots!(fig, state::IcaComponentState)
    state.fig = fig
    n_total = size(state.component_data, 1)
    num_plots = state.n_visible_components

    valid_comps = Int[]
    comp_to_k = Dict{Int,Int}()
    for i = 1:num_plots
        c = _get_component_index(state, i)
        if c <= n_total && !haskey(comp_to_k, c)
            push!(valid_comps, c)
            comp_to_k[c] = length(valid_comps)
        end
    end
    all_topo_data, x_bounds, y_bounds =
        _prepare_ica_topo_data(state.ica, valid_comps, state.plot_kwargs[:method], state.plot_kwargs[:gridscale])
    has_explicit_cr = !isnothing(get(state.plot_kwargs, :colorrange, nothing))
    levels_result =
        (!isempty(all_topo_data) && !has_explicit_cr) ?
        _calculate_ica_topo_levels(all_topo_data, state.use_global_scale[], get(state.plot_kwargs, :num_levels, 20)) : nothing

    init_xlims = (state.dat.data.time[first(state.xrange[])], state.dat.data.time[last(state.xrange[])])

    for i = 1:num_plots
        topo_ax = Axis(fig[i, 1], aspect = DataAspect(), titlevisible = false)
        deregister_interaction!(topo_ax, :rectanglezoom)
        deregister_interaction!(topo_ax, :scrollzoom)
        deregister_interaction!(topo_ax, :dragpan)
        push!(state.topo_axs, topo_ax)

        comp_idx = _get_component_index(state, i)
        is_last = i == state.n_visible_components

        ax = Axis(
            fig[i, 2],
            ylabel = @sprintf("IC %d", comp_idx),
            yaxisposition = :left,
            yticklabelsvisible = false,
            yticksvisible = true,
            xticklabelsvisible = is_last,
            xticksvisible = is_last,
            xgridvisible = false,
            ygridvisible = false,
            xminorgridvisible = false,
            yminorgridvisible = false,
            ylabelpadding = 0.0,
            yticklabelpad = 0.0,
            yticklabelspace = 0.0,
            yautolimitmargin = (0, 0),
        )
        deregister_interaction!(ax, :rectanglezoom)
        deregister_interaction!(ax, :dragpan)
        push!(state.axs, ax)
        push!(state.channel_axs, nothing)
        push!(state.overlay.plots, Any[])

        comp_data = comp_idx <= n_total ? state.component_data[comp_idx, :] : zeros(Float64, size(state.component_data, 2))
        lines_obs = Observable(comp_data)
        push!(state.lines_obs, lines_obs)

        lines!(
            ax,
            @lift(@view(state.dat.data.time[$(state.xrange)])),
            @lift(@view($(lines_obs)[$(state.xrange)])),
            color = _resolve_theme_linecolor(ax, nothing),
        )
        xlims!(ax, init_xlims)

        if comp_idx <= n_total
            k = comp_to_k[comp_idx]
            level_to_use = isnothing(levels_result) ? nothing : (state.use_global_scale[] ? levels_result : levels_result[k])
            _plot_ica_topo_in_viewer!(
                fig,
                topo_ax,
                state.ica,
                comp_idx;
                use_global_scale = state.use_global_scale[],
                pre_calculated_levels = level_to_use,
                precomputed_topo = (all_topo_data[k], x_bounds, y_bounds),
                state.plot_kwargs...,
            )
            topo_ax.ylabel = @sprintf("IC %d\n(%.1f%%)", comp_idx, state.ica.variance[comp_idx] * 100)
            topo_ax.ylabelvisible = true
        else
            empty!(topo_ax)
            topo_ax.ylabel = @sprintf("Invalid IC %d", comp_idx)
        end

        plot_index = i
        on(events(topo_ax).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                mouse_pos = mouseposition(topo_ax)
                head_radius = get(state.plot_kwargs, :head_radius, 1.0)
                if hypot(mouse_pos[1], mouse_pos[2]) <= head_radius
                    current_comp_idx = _get_component_index(state, plot_index)
                    plot_ica_component_spectrum(state.dat, state.ica, component_selection = components(current_comp_idx))
                end
            end
        end
    end
end

function _add_navigation_controls!(fig, state)
    topo_nav = GridLayout(fig[state.n_visible_components+1, 1], tellheight = false)
    colsize!(topo_nav, 1, 40)

    prev_topo = Button(topo_nav[1, 2], label = "◄ Previous", tellheight = false)
    next_topo = Button(topo_nav[1, 3], label = "Next ►", tellheight = false)

    Label(topo_nav[2, 2], "Components:", tellheight = false, halign = :right)
    text_input = Textbox(topo_nav[2, 3], placeholder = "e.g. 1,3-5,8", tellheight = false)
    apply_button = Button(topo_nav[2, 4], label = "Apply", tellheight = false)
    state.artifacts.textbox = text_input

    global_scale_check = Checkbox(topo_nav[3, 2], checked = state.use_global_scale[], tellheight = false)
    Label(topo_nav[3, 3], "Use Global Scale", tellwidth = false, tellheight = false)

    invert_scale_check = Checkbox(topo_nav[4, 2], checked = state.invert_scale[], tellheight = false)
    Label(topo_nav[4, 3], "Invert Scale", tellwidth = false, tellheight = false)

    colgap!(topo_nav, 2, 10)
    colgap!(topo_nav, 3, 5)
    rowgap!(topo_nav, 1, 35)
    rowgap!(topo_nav, 2, 45)
    rowgap!(topo_nav, 3, 35)

    on(global_scale_check.checked) do checked
        state.use_global_scale[] = checked
        _update_components!(state)
    end
    on(invert_scale_check.checked) do checked
        state.invert_scale[] = checked
        _update_components!(state)
    end
    on(prev_topo.clicks) do _
        state.comp_start[] = max(1, state.comp_start[] - state.n_visible_components)
        _update_components!(state)
    end
    on(next_topo.clicks) do _
        max_start = max(1, length(state.components) - state.n_visible_components + 1)
        state.comp_start[] = min(max_start, state.comp_start[] + state.n_visible_components)
        _update_components!(state)
    end
    on(apply_button.clicks) do _
        text_value = text_input.displayed_string[]
        if !isempty(text_value) && !startswith(text_value, "e.g.")
            comps = _parse_string_to_ints(text_value, size(state.component_data, 1))
            if !isempty(comps)
                @info "Creating new plot with components: $comps"
                plot_ica_component_activation(state.dat, state.ica, component_selection = components(comps))
            end
        end
    end
end

function _add_navigation_sliders!(fig, state)
    slider_row = state.n_visible_components + 3
    step_size = max(1, div(length(state.dat.data.time), 100))
    x_slider = Slider(
        fig[slider_row, 2],
        range = 1:step_size:length(state.dat.data.time),
        startvalue = first(state.xrange[]),
        tellwidth = false,
        tellheight = false,
        width = Auto(),
    )
    on(x_slider.value) do x
        _update_view_range!(state, Int(round(x)))
    end
end

function _update_view_range!(state, start_pos)
    n_t = length(state.dat.data.time)
    start_pos + state.window_size > n_t && (start_pos = max(1, n_t - state.window_size + 1))
    end_pos = min(n_t, start_pos + state.window_size - 1)
    state.xrange[] = start_pos:end_pos

    new_xlims = (state.dat.data.time[start_pos], state.dat.data.time[end_pos])
    for ax in state.axs
        xlims!(ax, new_xlims)
    end
    for ax in state.channel_axs
        !isnothing(ax) && xlims!(ax, new_xlims)
    end
end

function _show_channel_popup(state)
    all_channels = Symbol.(names(state.dat.data))
    _show_checkbox_dialog(
        "Select Additional Channels",
        all_channels,
        9,
        (ch, _) -> ch in state.overlay.names,
        [("None (Clear)", _ -> false)],
        cbs -> _update_channel_overlays!(state, all_channels[findall(cb -> cb.checked[], cbs)]);
        apply_width = 150,
    )
end

function _add_channel_menu!(fig, state)
    btn = Button(fig[state.n_visible_components+3, 1], label = "Additional Channels", fontsize = 18, tellheight = false)
    on(btn.clicks) do _
        _show_channel_popup(state)
    end
    return btn
end

function _clear_channel_overlays!(state)
    for i in eachindex(state.overlay.plots)
        ax_ch = state.channel_axs[i]
        if !isnothing(ax_ch)
            for p in state.overlay.plots[i]
                delete!(ax_ch, p)
            end
        end
        empty!(state.overlay.plots[i])
    end
    empty!(state.overlay.names)
end

function _update_channel_overlays!(state, selected_channels::Vector{Symbol})
    _clear_channel_overlays!(state)
    isempty(selected_channels) && return
    append!(state.overlay.names, selected_channels)

    n_panels = min(state.n_visible_components, length(state.channel_axs))
    for (ch_idx, ch_name) in enumerate(selected_channels)
        color = _CHANNEL_OVERLAY_COLORS[mod1(ch_idx, length(_CHANNEL_OVERLAY_COLORS))]
        col_data = state.dat.data[!, ch_name]
        is_bool = eltype(col_data) == Bool
        true_times = is_bool ? state.dat.data.time[findall(col_data)] : Float64[]
        ch_vec = is_bool ? Float64[] : Vector{Float64}(col_data)

        for i = 1:n_panels
            ax_ch = _ensure_channel_overlay_axis!(state, i)
            if is_bool
                if !isempty(true_times)
                    p = vlines!(ax_ch, true_times, color = (color, 0.5), linewidth = 1)
                    push!(state.overlay.plots[i], p)
                end
            else
                p = lines!(
                    ax_ch,
                    @lift(@view(state.dat.data.time[$(state.xrange)])),
                    @lift(@view(ch_vec[$(state.xrange)]) .* $(state.channel_yscale)),
                    color = color,
                )
                push!(state.overlay.plots[i], p)
            end
        end
    end
end

function _get_artifact_categories(artifacts::ArtifactComponents)
    categories = Tuple{Symbol,String,Vector{Int}}[]
    for prop in propertynames(artifacts)
        val = getproperty(artifacts, prop)
        if val isa AbstractDict
            for (k, v) in val
                push!(categories, (Symbol(k), string(k), collect(Int, v)))
            end
        elseif val isa AbstractVector
            lbl = join(titlecase.(split(string(prop), "_")), " ")
            push!(categories, (prop, lbl, collect(Int, val)))
        end
    end
    return categories
end

function _add_artifact_category_checkboxes!(fig, state, artifacts)
    state.artifacts.data = artifacts
    state.artifacts.toggles = Dict{Symbol,Observable{Bool}}()
    topo_nav = content(fig[state.n_visible_components+1, 1])

    categories = _get_artifact_categories(artifacts)
    for (idx, (key, label, comps)) in enumerate(categories)
        row = 4 + idx
        comp_str = isempty(comps) ? "none" : join(string.(comps), ", ")
        state.artifacts.toggles[key] = Observable(false)
        checkbox = Checkbox(topo_nav[row, 2], checked = false, tellheight = false)
        Label(topo_nav[row, 3], "$label ($comp_str)", tellwidth = false, tellheight = false)
        on(checkbox.checked) do checked
            state.artifacts.toggles[key][] = checked
            _update_artifact_textbox!(state)
        end
    end
    rowgap!(topo_nav, 4, 45)
    return nothing
end

function _update_artifact_textbox!(state)
    isnothing(state.artifacts.textbox) && return
    comps = Int[]
    for (key, _, source) in _get_artifact_categories(state.artifacts.data)
        state.artifacts.toggles[key][] && append!(comps, source)
    end
    state.artifacts.textbox.displayed_string[] = join(string.(sort!(unique!(comps))), ",")
    return nothing
end

function _setup_keyboard_interactions!(fig, state)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            _show_plot_help(:ica)
        elseif event.action in (Keyboard.press, Keyboard.repeat)
            if event.key == Keyboard.left || event.key == Keyboard.right
                current_range = state.xrange[]
                data_length = size(state.component_data, 2)
                window_size = state.window_size
                step_start = event.key == Keyboard.left ? (first(current_range) - window_size) : (first(current_range) + window_size)
                new_start = clamp(step_start, 1, max(1, data_length - window_size + 1))
                _update_view_range!(state, new_start)
            elseif event.key == Keyboard.up || event.key == Keyboard.down
                shift_pressed = (Keyboard.left_shift in events(fig).keyboardstate) || (Keyboard.right_shift in events(fig).keyboardstate)
                if !shift_pressed
                    zoom_step = get(state.plot_kwargs, :zoom_step, 0.2)
                    if event.key == Keyboard.up
                        _ymore!.(state.axs, zoom_step)
                    else
                        _yless!.(state.axs, zoom_step)
                    end
                    state.ylims[] = state.axs[1].yaxis.attributes.limits[]
                    num_plots = state.n_visible_components
                    for i = 1:num_plots
                        comp_idx = _get_component_index(state, i)
                        if comp_idx <= size(state.component_data, 1)
                            cdata = state.component_data[comp_idx, :]
                            state.lines_obs[i][] = state.invert_scale[] ? -cdata : cdata
                        end
                    end
                else
                    state.channel_yscale[] = event.key == Keyboard.up ? state.channel_yscale[] * 1.1 : state.channel_yscale[] / 1.1
                end
            elseif event.key == Keyboard.page_up || event.key == Keyboard.page_down
                if length(state.components) > state.n_visible_components
                    current_start = state.comp_start[]
                    max_start = max(1, length(state.components) - state.n_visible_components + 1)
                    new_start =
                        event.key == Keyboard.page_up ? max(1, current_start - state.n_visible_components) :
                        min(max_start, current_start + state.n_visible_components)
                    if new_start != current_start
                        state.comp_start[] = new_start
                        _update_components!(state)
                    end
                end
            end
        end
    end
end

function _update_components!(state)
    n_total = size(state.component_data, 1)
    num_plots = state.n_visible_components

    valid_comp_indices = Int[]
    plot_indices = Int[]
    for i = 1:num_plots
        comp_idx = _get_component_index(state, i)
        if comp_idx <= n_total
            push!(valid_comp_indices, comp_idx)
            push!(plot_indices, i)
        end
    end

    all_data, x_bounds, y_bounds =
        _prepare_ica_topo_data(state.ica, valid_comp_indices, state.plot_kwargs[:method], state.plot_kwargs[:gridscale])
    levels_result = _calculate_ica_topo_levels(all_data, state.use_global_scale[], state.plot_kwargs[:num_levels])

    for (k, i) in enumerate(plot_indices)
        comp_idx = valid_comp_indices[k]
        cdata = state.component_data[comp_idx, :]
        if i <= length(state.lines_obs)
            state.lines_obs[i][] = state.invert_scale[] ? -cdata : cdata
        end
        if i <= length(state.topo_axs)
            topo_ax = state.topo_axs[i]
            level_to_use = state.use_global_scale[] ? levels_result : levels_result[k]
            empty!(topo_ax)
            _plot_ica_topo_in_viewer!(
                state.fig,
                topo_ax,
                state.ica,
                comp_idx;
                pre_calculated_levels = level_to_use,
                precomputed_topo = (all_data[k], x_bounds, y_bounds),
                state.plot_kwargs...,
            )
            topo_ax.ylabel = @sprintf("IC %d\n(%.1f%%)", comp_idx, state.ica.variance[comp_idx] * 100)
            topo_ax.ylabelvisible = true
        end
        i <= length(state.axs) && (state.axs[i].ylabel = @sprintf("IC %d", comp_idx))
    end

    valid_plot_set = Set(plot_indices)
    for i = 1:num_plots
        if !(i in valid_plot_set)
            comp_idx = _get_component_index(state, i)
            i <= length(state.lines_obs) && (state.lines_obs[i][] = zeros(Float64, size(state.component_data, 2)))
            if i <= length(state.topo_axs)
                empty!(state.topo_axs[i])
                state.topo_axs[i].ylabel = @sprintf("Invalid IC %d", comp_idx)
                state.topo_axs[i].ylabelvisible = true
            end
            i <= length(state.axs) && (state.axs[i].ylabel = "")
        end
    end
end

function _get_component_index(state, i::Int)
    idx = state.comp_start[] + i - 1
    return idx <= length(state.components) ? state.components[idx] : typemax(Int)
end
