const PLOT_EPOCHS_REJECTION_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :display_plot => (true, "Whether to display the plot"),
    :dims => ((4, 6), "Grid dimensions as (rows, cols) for epoch display"),
    :spine_width => (2, "Width of axis spines"),
    :theme_fontsize => (nothing, "Font size for theme"),
    :linewidth => (nothing, "Line width for epoch traces"),
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    :colormap => (nothing, "Colormap for channel traces"),
    :good_epoch_color => (:green, "Color for good epoch spines"),
    :bad_epoch_color => (:red, "Color for bad epoch spines"),
    :global_ylim => (false, "Whether to use global Y-axis limits across all epochs"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple, or nothing for automatic"),
    :xlim => (nothing, "X-axis limits as (min, max) tuple, or nothing for automatic"),
)

"""
    EpochRejectionState

Stores the state of the epoch rejection interface.
"""
mutable struct EpochRejectionState
    epoch_data::EpochData
    selected_channels::Vector{Symbol}
    selected_channels_set::Set{Symbol}
    rejected::Vector{Bool}
    current_page::Observable{Int}
    epochs_per_page::Int
    n_total_epochs::Int
    n_pages::Int
    checkboxes::Vector{Toggle}
    axes::Vector{Axis}
    fig::Figure
    show_bad_channels_only::Observable{Bool}
    artifact_info::Union{Nothing,EpochRejectionInfo}
    colors::Vector{Any}
    global_ylim::Union{Nothing,Tuple{Float64,Float64}}
    is_updating_toggles::Base.RefValue{Bool}
end

function Base.getproperty(state::EpochRejectionState, sym::Symbol)
    sym === :epoch_axes && return getfield(state, :axes)
    return getfield(state, sym)
end

Base.broadcastable(state::EpochRejectionState) = Ref(state)
Base.iterate(state::EpochRejectionState, idx::Int = 1) = idx == 1 ? (state.fig, 2) : idx == 2 ? (state.axes, 3) : nothing

function Base.show(io::IO, state::EpochRejectionState)
    rejected_count = sum(state.rejected)
    println(io, "EpochRejectionState:")
    println(io, "  Total epochs: $(state.n_total_epochs)")
    println(io, "  Rejected epochs: $rejected_count ($(round(100 * rejected_count / state.n_total_epochs, digits=1))%)")
    println(io, "  Channels displayed: $(length(state.selected_channels))")
end

@inline function _set_axis_spine_color!(ax::Axis, color)
    ax.leftspinecolor = color
    ax.rightspinecolor = color
    ax.bottomspinecolor = color
    ax.topspinecolor = color
    return nothing
end

"""
    detect_bad_epochs_interactive(dat::EpochData;
                                  channel_selection::Function = channels(),
                                  artifact_info::Union{Nothing,EpochRejectionInfo} = nothing,
                                  kwargs...)::EpochRejectionState

Create an interactive interface for visually rejecting epochs.

Displays epochs in a grid layout with toggles for marking each epoch as bad (checked)
or good (unchecked) and navigation buttons for scrolling through pages of epochs.

# Arguments
- `dat::EpochData`: Epoched EEG data to review
- `channel_selection::Function`: Channels to display (default: all channels)
- `artifact_info::Union{Nothing,EpochRejectionInfo}`: Optional artifact detection info for bad channel filtering
- `kwargs`: Additional keyword arguments

$(_generate_kwargs_doc(PLOT_EPOCHS_REJECTION_KWARGS))

# Returns
- `EpochRejectionState`: State object containing rejection decisions

# Examples
```julia
state = detect_bad_epochs_interactive(epochs)
rejected_indices = get_rejected_epochs(state)
cleaned = reject_epochs(epochs, state)
```
"""
function detect_bad_epochs_interactive(
    dat::EpochData;
    channel_selection = channels(),
    artifact_info::Union{Nothing,EpochRejectionInfo} = nothing,
    kwargs...,
)::EpochRejectionState
    @info "Starting interactive epoch rejection interface"
    _set_window_title(_generate_window_title(dat))

    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_REJECTION_KWARGS, _normalize_plot_kwargs(kwargs))
    isempty(dat.data) && @minimal_error("Cannot create rejection interface for empty EpochData")
    dims = plot_kwargs[:dims]
    dims[1] * dims[2] <= 1 && @minimal_error("dims must be positive, got $dims")

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for display")

    n_total_epochs = n_epochs(dat)
    epochs_per_page = dims[1] * dims[2]
    n_pages = cld(n_total_epochs, epochs_per_page)

    @info "Displaying $(length(selected_channels)) channels"
    @info "Total epochs: $n_total_epochs across $n_pages pages"

    rejected = fill(false, n_total_epochs)
    if !isnothing(artifact_info)
        if hasproperty(artifact_info, :rejected) && !isnothing(artifact_info.rejected)
            for r in artifact_info.rejected
                1 <= r.epoch <= n_total_epochs && (rejected[r.epoch] = true)
            end
        elseif hasproperty(artifact_info, :rejected_epochs) && !isnothing(artifact_info.rejected_epochs)
            for ep in artifact_info.rejected_epochs
                1 <= ep <= n_total_epochs && (rejected[ep] = true)
            end
        end
    end

    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; figure_padding = 50, fontsize_kw...)

    selected_channels_set = Set(selected_channels)
    colors = collect(Any, Makie.to_colormap(_resolve_theme_colormap(fig, plot_kwargs[:colormap])))
    global_ylim = plot_kwargs[:global_ylim] ? _calculate_global_ylim(dat, selected_channels) : nothing

    state = EpochRejectionState(
        dat,
        selected_channels,
        selected_channels_set,
        rejected,
        Observable(1),
        epochs_per_page,
        n_total_epochs,
        n_pages,
        Toggle[],
        Axis[],
        fig,
        Observable(false),
        artifact_info,
        colors,
        global_ylim,
        Ref(false),
    )

    _create_rejection_interface!(fig, state, dims, artifact_info, plot_kwargs)
    plot_kwargs[:display_plot] && _display_figure(fig)

    @info "Interface ready. Review epochs and mark bad ones by checking the boxes."
    @info "Checked = REJECT, Unchecked = KEEP"
    return state
end

function detect_bad_epochs_interactive(
    epochs::AbstractVector{<:EpochData};
    channel_selection = channels(),
    artifact_info::Union{Nothing,AbstractVector{<:EpochRejectionInfo}} = nothing,
    kwargs...,
)
    return [
        detect_bad_epochs_interactive(
            epochs[i];
            channel_selection = channel_selection,
            artifact_info = isnothing(artifact_info) ? nothing : artifact_info[i],
            kwargs...,
        ) for i in eachindex(epochs)
    ]
end

function _create_rejection_interface!(
    fig::Figure,
    state::EpochRejectionState,
    dims::Tuple{Int,Int},
    artifact_info::Union{Nothing,EpochRejectionInfo},
    plot_kwargs::Dict{Symbol,Any},
)
    rows, cols = dims
    @info "Creating $rows x $cols grid of epoch plots"

    empty!(state.axes)
    empty!(state.checkboxes)

    root = GridLayout(fig[1, 1])
    axes_gl = GridLayout(root[1, 1])

    for row = 1:rows, col = 1:cols
        idx = (row - 1) * cols + col
        cell_gl = GridLayout(axes_gl[row, col])
        ax = Axis(cell_gl[1, 1])
        push!(state.axes, ax)

        t = Toggle(cell_gl[1, 1]; halign = :left, valign = :top, tellwidth = false, tellheight = false, length = 60, markersize = 30)
        ax.spinewidth = plot_kwargs[:spine_width]
        ax.xgridvisible = plot_kwargs[:xgrid]
        ax.ygridvisible = plot_kwargs[:ygrid]
        ax.xminorgridvisible = plot_kwargs[:xminorgrid]
        ax.yminorgridvisible = plot_kwargs[:yminorgrid]
        colsize!(cell_gl, 1, Relative(1))

        on(t.active) do is_checked
            _set_axis_spine_color!(ax, is_checked ? plot_kwargs[:bad_epoch_color] : plot_kwargs[:good_epoch_color])
            state.is_updating_toggles[] && return
            epoch_idx = (state.current_page[] - 1) * state.epochs_per_page + idx
            1 <= epoch_idx <= state.n_total_epochs && (state.rejected[epoch_idx] = is_checked)
        end
        push!(state.checkboxes, t)
    end

    nav_gl = GridLayout(root[2, 1])
    btn_prev = Button(nav_gl[1, 1], label = "◀ Prev", width = 100, height = 40)
    on(btn_prev.clicks) do _
        if state.current_page[] > 1
            state.current_page[] -= 1
            _update_epoch_display!(state, artifact_info, plot_kwargs)
        end
    end
    btn_next = Button(nav_gl[1, 2], label = "Next ▶", width = 100, height = 40)
    on(btn_next.clicks) do _
        if state.current_page[] < state.n_pages
            state.current_page[] += 1
            _update_epoch_display!(state, artifact_info, plot_kwargs)
        end
    end

    bad_channels_toggle = Toggle(nav_gl[1, 3])
    Label(nav_gl[1, 4], "Show bad channels only")
    on(bad_channels_toggle.active) do active
        if active && (isnothing(artifact_info) || isempty(artifact_info.rejected))
            @info "Bad channel filtering requires artifact_info with rejected epochs"
        else
            state.show_bad_channels_only[] = active
            _update_epoch_display!(state, artifact_info, plot_kwargs)
        end
    end

    global_ylim_toggle = Toggle(nav_gl[1, 5])
    Label(nav_gl[1, 6], "Global Y limits")
    on(global_ylim_toggle.active) do active
        plot_kwargs[:global_ylim] = active
        if active && isnothing(state.global_ylim)
            state.global_ylim = _calculate_global_ylim(state.epoch_data, state.selected_channels)
        end
        _update_epoch_display!(state, artifact_info, plot_kwargs)
    end

    colsize!(root, 1, Relative(1))
    _update_epoch_display!(state, artifact_info, plot_kwargs)
end

function _epoch_channels_to_show(state::EpochRejectionState, epoch_idx::Int)
    state.show_bad_channels_only[] || return state.selected_channels
    info = state.artifact_info
    isnothing(info) && return Symbol[]
    return [ch for ch in state.selected_channels if any(r -> r.epoch == epoch_idx && r.channel == ch, info.rejected)]
end

function _update_epoch_display!(state::EpochRejectionState, artifact_info::Union{Nothing,EpochRejectionInfo}, plot_kwargs::Dict{Symbol,Any})
    start_idx = (state.current_page[] - 1) * state.epochs_per_page + 1

    state.is_updating_toggles[] = true
    try
        for (i, ax) in enumerate(state.axes)
            epoch_idx = start_idx + i - 1
            empty!(ax)
            if epoch_idx <= state.n_total_epochs
                channels_to_show = _epoch_channels_to_show(state, epoch_idx)
                _plot_single_epoch!(ax, state, epoch_idx, channels_to_show, plot_kwargs)
                ax.title = "Epoch $epoch_idx: $(_print_vector(channels_to_show, n_ends = 3))"

                if !isnothing(plot_kwargs[:ylim])
                    ylims!(ax, plot_kwargs[:ylim])
                elseif plot_kwargs[:global_ylim] && !isnothing(state.global_ylim)
                    ylims!(ax, state.global_ylim)
                else
                    ax.limits = (ax.limits[][1], nothing)
                end
                !isnothing(plot_kwargs[:xlim]) && xlims!(ax, plot_kwargs[:xlim])

                is_rej = state.rejected[epoch_idx]
                i <= length(state.checkboxes) && (state.checkboxes[i].active[] = is_rej)
                _set_axis_spine_color!(ax, is_rej ? plot_kwargs[:bad_epoch_color] : plot_kwargs[:good_epoch_color])
            else
                ax.title = ""
                i <= length(state.checkboxes) && (state.checkboxes[i].active[] = false)
                _set_axis_spine_color!(ax, :transparent)
            end
        end
    finally
        state.is_updating_toggles[] = false
    end
    return nothing
end

function _plot_single_epoch!(
    ax::Axis,
    state::EpochRejectionState,
    epoch_idx::Int,
    channels_to_plot::Vector{Symbol},
    plot_kwargs::Dict{Symbol,Any},
)
    epoch = state.epoch_data.data[epoch_idx]
    t = epoch[!, :time]
    lw = _resolve_theme_linewidth(ax, plot_kwargs[:linewidth], 1)
    n_colors = length(state.colors)

    for (ch_idx, ch) in enumerate(channels_to_plot)
        lines!(ax, t, epoch[!, ch], color = state.colors[mod1(ch_idx, n_colors)], linewidth = lw)
    end
    _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
end

@inline _col_minmax(col::AbstractVector) = (Float64(minimum(col)), Float64(maximum(col)))

function _calculate_global_ylim(dat::EpochData, selected_channels::Vector{Symbol})::Tuple{Float64,Float64}
    global_min = Inf
    global_max = -Inf
    for epoch in dat
        for ch in selected_channels
            if hasproperty(epoch, ch)
                lo, hi = _col_minmax(epoch[!, ch])
                global_min = min(global_min, lo)
                global_max = max(global_max, hi)
            end
        end
    end
    (isinf(global_min) || isinf(global_max)) && return (-1.0, 1.0)
    padding = (global_max - global_min) * 0.05
    return (global_min - padding, global_max + padding)
end

"""
    get_rejected_epochs(state::EpochRejectionState)::Vector{Int}

Get indices of rejected epochs from the rejection state.

# Examples
```julia
state = detect_bad_epochs_interactive(epochs)
rejected_indices = get_rejected_epochs(state)
```
"""
get_rejected_epochs(state::EpochRejectionState)::Vector{Int} = findall(state.rejected)

function _to_rejection_info(state::EpochRejectionState)::EpochRejectionInfo
    dat = state.epoch_data
    rejected_epoch_indices = findall(state.rejected)
    layout_channels = state.selected_channels

    rejections = Vector{Rejection}(undef, length(rejected_epoch_indices) * length(layout_channels))
    k = 0
    @inbounds for epoch_idx in rejected_epoch_indices, ch in layout_channels
        k += 1
        rejections[k] = Rejection(ch, epoch_idx)
    end

    info = EpochInfo(dat.condition, dat.condition_name, state.n_total_epochs)
    return EpochRejectionInfo("interactive_rejection", info, length(rejections), 0, nothing, 0, nothing, rejections, nothing, nothing)
end

_to_rejection_info(states::Vector{EpochRejectionState})::Vector{EpochRejectionInfo} = _to_rejection_info.(states)

reject_epochs(dat::EpochData, state::EpochRejectionState)::EpochData = reject_epochs(dat, _to_rejection_info(state))
reject_epochs(dat::Vector{EpochData}, states::Vector{EpochRejectionState})::Vector{EpochData} =
    reject_epochs.(dat, _to_rejection_info.(states))
