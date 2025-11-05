"""
Interactive visual epoch rejection interface using Makie.

This module provides an interactive GUI for manually rejecting epochs based on
visual inspection. It displays epochs in a grid layout with checkboxes for
marking epochs as good/bad and navigation buttons for scrolling through data.
"""

#=============================================================================
    DEFAULT KEYWORD ARGUMENTS
=============================================================================#
const PLOT_EPOCHS_REJECTION_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Grid configuration
    :dims => ((4, 6), "Grid dimensions as (rows, cols) for epoch display"),

    # Plot styling
    :spine_width => (2, "Width of axis spines"),
    :theme_fontsize => (20, "Font size for theme"),

    # Line styling
    :linewidth => (1, "Line width for epoch traces"),
    :origin_linewidth => (0.5, "Line width for origin line at y=0"),
    :origin_linestyle => (:dash, "Line style for origin line"),
    :origin_color => (:black, "Color for origin line"),

    # Origin lines
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Color scheme
    :colormap => (:GnBu, "Colormap for channel traces"),
    :good_epoch_color => (:green, "Color for good epoch spines"),
    :bad_epoch_color => (:red, "Color for bad epoch spines"),

    # Y-axis limits
    :global_ylim => (false, "Whether to use global Y-axis limits across all epochs"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple, or nothing for automatic"),
    :xlim => (nothing, "X-axis limits as (min, max) tuple, or nothing for automatic"),
)

#=============================================================================
    EPOCH REJECTION STATE
=============================================================================#

"""
    EpochRejectionState

Stores the state of the epoch rejection interface.
"""
mutable struct EpochRejectionState
    epoch_data::EpochData
    selected_channels::Vector{Symbol}
    selected_channels_set::Set{Symbol}  # Precomputed set for faster intersection
    rejected::Vector{Bool}  # One per epoch
    current_page::Observable{Int}
    epochs_per_page::Int
    n_total_epochs::Int
    n_pages::Int
    checkboxes::Vector{Toggle}
    epoch_axes::Vector{Axis}
    fig::Figure
    show_bad_channels_only::Observable{Bool}  # Filter to show only bad channels
    artifact_info::Union{Nothing,EpochRejectionInfo}  # Store artifact info for per-epoch filtering
    colors::Vector{Any}  # Cached color array
    global_ylim::Union{Nothing,Tuple{Float64,Float64}}  # Global Y-axis limits for all epochs
end


#=============================================================================
    MAIN PLOTTING FUNCTION
=============================================================================#

"""
    detect_bad_epochs_interactive(dat::EpochData;
                             channel_selection::Function = channels(),
                             epochs_per_page::Int = 12,
                             grid_size::Tuple{Int,Int} = (3, 4))::EpochRejectionState

Create an interactive interface for visually rejecting epochs.

This function creates a Makie window displaying epochs in a grid layout with
checkboxes for marking each epoch as good (✓) or bad (✗). Navigation buttons
allow scrolling through pages of epochs.

# Arguments
- `dat::EpochData`: Epoched EEG data to review
- `channel_selection::Function`: Channels to display (default: all channels)
- `artifact_info::Union{Nothing,EpochRejectionInfo}`: Optional artifact detection info for bad channel filtering
- `kwargs`: Additional keyword arguments

$(generate_kwargs_doc(PLOT_EPOCHS_REJECTION_KWARGS))

# Returns
- `EpochRejectionState`: State object containing rejection decisions

# Usage
After reviewing all epochs and marking bad ones, extract the results:
```julia
# Create interactive interface
state = detect_bad_epochs_interactive(epochs)

# After closing the window or finishing review:
rejected_indices = findall(state.rejected)
clean_epochs = epochs.data[.!state.rejected]
```

# Examples
```julia
using eegfun, JLD2

# Load epochs
epochs = load("participant_1_epochs.jld2", "epochs")

# Launch interactive rejection interface
state = detect_bad_epochs_interactive(epochs)

# Review epochs by:
# 1. Checking/unchecking boxes for each epoch
# 2. Using Previous/Next buttons to navigate
# 3. Using First/Last buttons to jump

# After review, get results
rejected_indices = findall(state.rejected)
println("Rejected epochs: \$rejected_indices")

# Create cleaned dataset
clean_data = EpochData(
    epochs.data[.!state.rejected],
    epochs.layout,
    epochs.sample_rate,
    epochs.analysis_info
)

# Save
jldsave("participant_1_epochs_cleaned.jld2"; data = clean_data)
```

# Interactive Controls
- **Checkboxes**: Toggle to mark epoch as bad (checked) or good (unchecked)
- **Previous/Next**: Navigate one page at a time
- **First/Last**: Jump to first or last page
- **Page indicator**: Shows current page number
- **Show bad channels only**: Filter to display only channels that are marked as bad in any epoch (only available when artifact_info is provided)

# Notes
- Epochs are displayed with all selected channels overlaid (or only bad channels if filter is enabled)
- Checked boxes indicate epochs to REJECT
- Unchecked boxes indicate epochs to KEEP
- Changes are saved in the state object
- Close the window when done reviewing
- Bad channels filter only appears when artifact_info is provided and contains rejected channels
"""
function detect_bad_epochs_interactive(
    dat::EpochData;
    channel_selection::Function = channels(),
    artifact_info::Union{Nothing,EpochRejectionInfo} = nothing,
    kwargs...,
)::EpochRejectionState

    @info "Starting interactive epoch rejection interface"

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_EPOCHS_REJECTION_KWARGS, kwargs)

    # Validate inputs
    _validate_rejection_gui_inputs(dat, plot_kwargs[:dims])

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for display")
    end

    n_total_epochs = n_epochs(dat)
    epochs_per_page = plot_kwargs[:dims][1] * plot_kwargs[:dims][2]
    n_pages = floor(Int, n_total_epochs / epochs_per_page)

    @info "Displaying $(length(selected_channels)) channels"
    @info "Total epochs: $n_total_epochs across $n_pages pages"

    # Initialize rejection state (all epochs start as not rejected)
    rejected = fill(false, n_total_epochs)

    # Create figure sized to fit typical screens
    fig = Figure(figure_padding = 50)

    # Performance optimizations
    selected_channels_set = Set(selected_channels)
    colors = Makie.to_colormap(plot_kwargs[:colormap])
    global_ylim = plot_kwargs[:global_ylim] ? _calculate_global_ylim(dat, selected_channels) : nothing

    # Create state object
    current_page = Observable(1)
    show_bad_channels_only = Observable(false)
    state = EpochRejectionState(
        dat,
        selected_channels,
        selected_channels_set,
        rejected,
        current_page,
        epochs_per_page,
        n_total_epochs,
        n_pages,
        Toggle[],
        Axis[],
        fig,
        show_bad_channels_only,
        artifact_info,
        colors,
        global_ylim,
    )

    # Create UI layout
    _create_rejection_interface!(fig, state, plot_kwargs[:dims], artifact_info, plot_kwargs)

    # Display figure
    if plot_kwargs[:display_plot]
        display(fig)
    end

    @info "Interface ready. Review epochs and mark bad ones by checking the boxes."
    @info "Checked = REJECT, Unchecked = KEEP"

    return state
end


#=============================================================================
    UI CREATION FUNCTIONS
=============================================================================#

"""
Create the main interface layout with epoch grid, checkboxes, and navigation.
"""
function _create_rejection_interface!(
    fig::Figure,
    state::EpochRejectionState,
    dims::Tuple{Int,Int},
    artifact_info::Union{Nothing,EpochRejectionInfo},
    plot_kwargs::Dict{Symbol,Any},
)
    rows, cols = dims

    @info "Creating $rows x $cols grid of epoch plots"

    # Clear any existing UI handles
    empty!(state.epoch_axes)
    empty!(state.checkboxes)

    # Root layout with two rows: axes grid (1) and nav bar (2)
    root = GridLayout(fig[1, 1])

    # Strict sublayout for axes: rows x cols
    axes_gl = GridLayout(root[1, 1])
    # Build per-cell sublayouts: axis on top, checkbox below
    for row = 1:rows
        for col = 1:cols
            idx = (row - 1) * cols + col
            cell_gl = GridLayout(axes_gl[row, col])
            ax = Axis(cell_gl[1, 1])
            push!(state.epoch_axes, ax)
            # Place a small toggle inside the plot area, top-left
            t = Toggle(
                cell_gl[1, 1];
                halign = :left,
                valign = :top,
                tellwidth = false,
                tellheight = false,
                length = 60,
                markersize = 30,
            )

            ax.spinewidth = plot_kwargs[:spine_width]

            # Configure grid
            ax.xgridvisible = plot_kwargs[:xgrid]
            ax.ygridvisible = plot_kwargs[:ygrid]
            ax.xminorgridvisible = plot_kwargs[:xminorgrid]
            ax.yminorgridvisible = plot_kwargs[:yminorgrid]

            # Set font sizes (only need to do this once)
            ax.titlesize = plot_kwargs[:theme_fontsize]
            ax.xlabelsize = plot_kwargs[:theme_fontsize]
            ax.ylabelsize = plot_kwargs[:theme_fontsize]
            ax.xticklabelsize = plot_kwargs[:theme_fontsize]
            ax.yticklabelsize = plot_kwargs[:theme_fontsize]

            on(t.active) do active
                color = active ? plot_kwargs[:bad_epoch_color] : plot_kwargs[:good_epoch_color]
                for spline in (:leftspinecolor, :rightspinecolor, :bottomspinecolor, :topspinecolor)
                    setproperty!(ax, spline, color)
                end
            end

            colsize!(cell_gl, 1, Relative(1))
            on(t.active) do is_checked
                page = state.current_page[]
                epoch_idx = (page - 1) * state.epochs_per_page + idx
                if epoch_idx <= state.n_total_epochs
                    state.rejected[epoch_idx] = is_checked
                end
            end
            push!(state.checkboxes, t)
        end
    end

    # Navigation sublayout
    nav_gl = GridLayout(root[2, 1])

    # Navigation buttons
    btn_prev = Button(nav_gl[1, 1], label = "◀ Prev", width = 100, height = 40)
    on(btn_prev.clicks) do _
        if state.current_page[] > 1
            state.current_page[] -= 1
            _update_epoch_display!(state, artifact_info, plot_kwargs)
        end
    end
    btn_next = Button(nav_gl[1, 2], label = "Next ▶", width = 100, height = 40)
    on(btn_next.clicks) do _
        if state.current_page[] <= state.n_pages
            state.current_page[] += 1
            _update_epoch_display!(state, artifact_info, plot_kwargs)
        end
    end

    # Bad channels filter checkbox (always create to maintain layout)
    bad_channels_toggle = Toggle(nav_gl[1, 3])
    bad_channels_label = Label(nav_gl[1, 4], "Show bad channels only")

    on(bad_channels_toggle.active) do active
        # Only show info message when trying to turn ON without artifact_info
        if active && (isnothing(artifact_info) || isempty(artifact_info.rejected))
            @info "Bad channel filtering requires artifact_info with rejected epochs"
        else # Normal operation: update state and display
            state.show_bad_channels_only[] = active
            _update_epoch_display!(state, artifact_info, plot_kwargs)
        end
    end

    # Global Y limits checkbox
    global_ylim_toggle = Toggle(nav_gl[1, 5])
    global_ylim_label = Label(nav_gl[1, 6], "Global Y limits")
    on(global_ylim_toggle.active) do active
        plot_kwargs[:global_ylim] = active
        # Recalculate global limits if enabling
        if active && isnothing(state.global_ylim)
            state.global_ylim = _calculate_global_ylim(state.epoch_data, state.selected_channels)
        end
        _update_epoch_display!(state, artifact_info, plot_kwargs)
    end

    # # Size root rows/cols now
    colsize!(root, 1, Relative(1))

    _update_epoch_display!(state, artifact_info, plot_kwargs)
end


"""
Update the display when page changes.
"""
function _update_epoch_display!(
    state::EpochRejectionState,
    artifact_info::Union{Nothing,EpochRejectionInfo},
    plot_kwargs::Dict{Symbol,Any},
)
    page = state.current_page[]
    start_idx = (page - 1) * state.epochs_per_page + 1

    for (i, ax) in enumerate(state.epoch_axes)
        epoch_idx = start_idx + i - 1
        empty!(ax)
        if epoch_idx <= state.n_total_epochs
            _plot_single_epoch!(ax, state, epoch_idx, plot_kwargs)

            # Update title to show filtered channels
            channels_to_show = if state.show_bad_channels_only[]
                bad_channels_for_epoch = _get_bad_channels_for_epoch(state.artifact_info, epoch_idx)
                collect(intersect(state.selected_channels_set, bad_channels_for_epoch))
            else
                state.selected_channels
            end

            ax.title = "Epoch $epoch_idx: $(print_vector(channels_to_show, n_ends = 3))"

            # Apply Y limits: explicit ylim takes precedence, then global_ylim, then automatic
            if !isnothing(plot_kwargs[:ylim])
                ylims!(ax, plot_kwargs[:ylim])
            elseif plot_kwargs[:global_ylim] && !isnothing(state.global_ylim)
                ylims!(ax, state.global_ylim)
            else
                # Reset to automatic scaling by clearing the limits
                ax.limits = (ax.limits[][1], nothing)
            end
            
            # Apply X limits if provided
            if !isnothing(plot_kwargs[:xlim])
                xlims!(ax, plot_kwargs[:xlim])
            end
            if i <= length(state.checkboxes)
                # Check if epoch is already rejected by user or by automatic detection
                auto_rejected = if isnothing(artifact_info)
                    false
                else
                    epoch_idx ∈ unique_epochs(artifact_info.rejected)
                end
                state.checkboxes[i].active[] = state.rejected[epoch_idx] || auto_rejected
            end
        else
            ax.title = ""
            state.checkboxes[i].active[] = false
        end
    end
    return nothing
end


"""
Plot a single epoch in the given axis.
"""
function _plot_single_epoch!(ax::Axis, state::EpochRejectionState, epoch_idx::Int, plot_kwargs::Dict{Symbol,Any})
    epoch = state.epoch_data.data[epoch_idx]
    t = epoch.time

    # Filter channels based on bad channels checkbox
    channels_to_plot = if state.show_bad_channels_only[]
        # Only show channels that are bad in this specific epoch
        bad_channels_for_epoch = _get_bad_channels_for_epoch(state.artifact_info, epoch_idx)
        collect(intersect(state.selected_channels_set, bad_channels_for_epoch))
    else
        # Show all selected channels
        state.selected_channels
    end

    for (ch_idx, ch) in enumerate(channels_to_plot)
        color = state.colors[mod1(ch_idx, length(state.colors))]
        lines!(ax, t, epoch[!, ch], color = color, linewidth = plot_kwargs[:linewidth])
    end

    # Add origin lines if enabled
    if plot_kwargs[:add_xy_origin]
        hlines!(
            ax,
            [0.0],
            color = plot_kwargs[:origin_color],
            linewidth = plot_kwargs[:origin_linewidth],
            linestyle = plot_kwargs[:origin_linestyle],
        )
        vlines!(
            ax,
            [0.0],
            color = plot_kwargs[:origin_color],
            linewidth = plot_kwargs[:origin_linewidth],
            linestyle = plot_kwargs[:origin_linestyle],
        )
    end
end

function Base.show(io::IO, state::EpochRejectionState)
    rejected_count = sum(state.rejected)
    println(io, "EpochRejectionState:")
    println(io, "  Total epochs: $(state.n_total_epochs)")
    println(io, "  Rejected epochs: $rejected_count ($(round(100 * rejected_count / state.n_total_epochs, digits=1))%)")
    println(io, "  Channels displayed: $(length(state.selected_channels))")
end


#=============================================================================
    INPUT VALIDATION
=============================================================================#

"""
Validate inputs for the rejection GUI.
"""
function _validate_rejection_gui_inputs(dat::EpochData, dims::Tuple{Int,Int})
    if isempty(dat.data)
        @minimal_error_throw("Cannot create rejection interface for empty EpochData")
    end
    if dims[1] * dims[2] <= 1
        @minimal_error_throw("dims must be positive, got $dims")
    end
end


#=============================================================================
    HELPER FUNCTIONS
=============================================================================#


"""
    _get_bad_channels_for_epoch(artifact_info::Union{Nothing,EpochRejectionInfo}, epoch_idx::Int)::Set{Symbol}

Get the set of bad channels for a specific epoch by filtering artifact info.
"""
function _get_bad_channels_for_epoch(artifact_info::Union{Nothing,EpochRejectionInfo}, epoch_idx::Int)::Set{Symbol}
    if isnothing(artifact_info)
        return Set{Symbol}()
    end

    bad_channels = [r.label for r in artifact_info.rejected if r.epoch == epoch_idx]
    return Set(bad_channels)
end

"""
    _calculate_global_ylim(dat::EpochData, selected_channels::Vector{Symbol})::Tuple{Float64,Float64}

Calculate global Y-axis limits across all epochs and selected channels.
"""
function _calculate_global_ylim(dat::EpochData, selected_channels::Vector{Symbol})::Tuple{Float64,Float64}
    min_vals = Float64[]
    max_vals = Float64[]

    for epoch in dat.data
        for ch in selected_channels
            if hasproperty(epoch, ch)
                ch_data = epoch[!, ch]
                push!(min_vals, minimum(ch_data))
                push!(max_vals, maximum(ch_data))
            end
        end
    end

    global_min = minimum(min_vals)
    global_max = maximum(max_vals)

    # Add some padding (5% on each side)
    padding = (global_max - global_min) * 0.05
    return (global_min - padding, global_max + padding)
end



#=============================================================================
    HELPER FUNCTIONS FOR EXTRACTING RESULTS
=============================================================================#

"""
    get_rejected_epochs(state::EpochRejectionState)::Vector{Int}

Get indices of rejected epochs from the rejection state.

# Examples
```julia
state = detect_bad_epochs_interactive(epochs)
# ... after review ...
rejected_indices = get_rejected_epochs(state)
```
"""
function get_rejected_epochs(state::EpochRejectionState)::Vector{Int}
    return findall(state.rejected)
end
