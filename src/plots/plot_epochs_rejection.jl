"""
Interactive visual epoch rejection interface using Makie.

This module provides an interactive GUI for manually rejecting epochs based on
visual inspection. It displays epochs in a grid layout with checkboxes for
marking epochs as good/bad and navigation buttons for scrolling through data.
"""

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
end


#=============================================================================
    MAIN PLOTTING FUNCTION
=============================================================================#

"""
    reject_epochs_interactive(dat::EpochData;
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
- `grid_size::Tuple{Int,Int}`: Grid dimensions as (rows, cols) (default: (3, 4) for 12 epochs)

# Returns
- `EpochRejectionState`: State object containing rejection decisions

# Usage
After reviewing all epochs and marking bad ones, extract the results:
```julia
# Create interactive interface
state = reject_epochs_interactive(epochs)

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
state = reject_epochs_interactive(epochs)

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
save("participant_1_epochs_cleaned.jld2", "epochs", clean_data)
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
function reject_epochs_interactive(
    dat::EpochData;
    channel_selection::Function = channels(),
    grid_size::Tuple{Int,Int} = (4, 6),
    artifact_info::Union{Nothing,EpochRejectionInfo} = nothing,
)::EpochRejectionState

    @info "Starting interactive epoch rejection interface"

    # Validate inputs
    _validate_rejection_gui_inputs(dat, grid_size)
    
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for display")
    end
    
    n_total_epochs = n_epochs(dat)
    println("n_total_epochs: $n_total_epochs")
    n_pages = floor(Int, n_total_epochs / (grid_size[1] * grid_size[2]))
    println("n_pages: $n_pages")
    
    @info "Displaying $(length(selected_channels)) channels"
    @info "Total epochs: $n_total_epochs across $n_pages pages"
    
    # Initialize rejection state (all epochs start as not rejected)
    rejected = fill(false, n_total_epochs)
    
    # Create figure sized to fit typical screens
    fig = Figure(figure_padding = 50)
    
    # Create state object
    current_page = Observable(1)
    show_bad_channels_only = Observable(false)
    state = EpochRejectionState(
        dat,
        selected_channels,
        rejected,
        current_page,
        grid_size[1] * grid_size[2],
        n_total_epochs,
        n_pages,
        Toggle[],
        Axis[],
        fig,
        show_bad_channels_only,
        artifact_info
    )
    
    # Create UI layout
    _create_rejection_interface!(fig, state, grid_size, artifact_info)
    
    # Display figure
    display(fig)
    
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
function _create_rejection_interface!(fig::Figure, state::EpochRejectionState, grid_size::Tuple{Int,Int}, artifact_info::Union{Nothing,EpochRejectionInfo})
    rows, cols = grid_size
    
    @info "Creating $rows x $cols grid of epoch plots"
    
    # Clear any existing UI handles
    empty!(state.epoch_axes)
    empty!(state.checkboxes)
    
    # Root layout with two rows: axes grid (1) and nav bar (2)
    root = GridLayout(fig[1, 1])

    # Strict sublayout for axes: rows x cols
    axes_gl = GridLayout(root[1, 1])
    # Build per-cell sublayouts: axis on top, checkbox below
    for row in 1:rows
        for col in 1:cols
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
                markersize=30
            )

            ax.spinewidth = 2
            on(t.active) do active
                color = active ? :red : :green
                for spline in (:leftspinecolor, :rightspinecolor, :bottomspinecolor, :topspinecolor)
                    setproperty!(ax, spline, color)
                end
                # if active 
                #     ax.leftspinecolor = :red
                #     ax.rightspinecolor = :red
                #     ax.bottomspinecolor = :red
                #     ax.topspinecolor = :red
                # else
                #     ax.leftspinecolor = :green
                #     ax.rightspinecolor = :green
                #     ax.bottomspinecolor = :green
                #     ax.topspinecolor = :green
                # end
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
    btn_prev = Button(nav_gl[1, 1], label = "◀ Prev")
    on(btn_prev.clicks) do _
        if state.current_page[] > 1
            state.current_page[] -= 1
            _update_epoch_display!(state, artifact_info)
        end
    end
    btn_next = Button(nav_gl[1, 2], label = "Next ▶")
    on(btn_next.clicks) do _
        if state.current_page[] <= state.n_pages
            state.current_page[] += 1
            _update_epoch_display!(state, artifact_info)
        end
    end
    
    # Bad channels filter checkbox
    if !isnothing(artifact_info) && !isempty(artifact_info.rejected_epochs)
        bad_channels_toggle = Toggle(nav_gl[1, 3])
        bad_channels_label = Label(nav_gl[1, 4], "Show bad channels only")
        on(bad_channels_toggle.active) do active
            state.show_bad_channels_only[] = active
            _update_epoch_display!(state, artifact_info)
        end
    end

    # # Size root rows/cols now
    colsize!(root, 1, Relative(1))
    
    _update_epoch_display!(state, artifact_info)
end


"""
Update the display when page changes.
"""
function _update_epoch_display!(state::EpochRejectionState, artifact_info::Union{Nothing,EpochRejectionInfo})
    page = state.current_page[]
    start_idx = (page - 1) * state.epochs_per_page + 1
    
    for (i, ax) in enumerate(state.epoch_axes)
        epoch_idx = start_idx + i - 1
        empty!(ax)
        if epoch_idx <= state.n_total_epochs
            _plot_single_epoch!(ax, state, epoch_idx)
            
            # Update title to show filtered channels
            channels_to_show = if state.show_bad_channels_only[]
                bad_channels_for_epoch = _get_bad_channels_for_epoch(state.artifact_info, epoch_idx)
                intersect(state.selected_channels, bad_channels_for_epoch)
            else
                state.selected_channels
            end
            println("channels_to_show: $channels_to_show")
            
            ax.title = "Epoch $epoch_idx: $(print_vector(channels_to_show, n_ends = 3))"
            ax.titlesize = 22
            if i <= length(state.checkboxes)
                state.checkboxes[i].active[] = state.rejected[epoch_idx] || epoch_idx ∈ unique_epochs(artifact_info.rejected_epochs)
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
function _plot_single_epoch!(ax::Axis, state::EpochRejectionState, epoch_idx::Int)
    epoch = state.epoch_data.data[epoch_idx]
    t = epoch.time
    colors = Makie.wong_colors()
    
    # Filter channels based on bad channels checkbox
    channels_to_plot = if state.show_bad_channels_only[]
        # Only show channels that are bad in this specific epoch
        bad_channels_for_epoch = _get_bad_channels_for_epoch(state.artifact_info, epoch_idx)
        intersect(state.selected_channels, bad_channels_for_epoch)
    else
        # Show all selected channels
        state.selected_channels
    end
    
    for (ch_idx, ch) in enumerate(channels_to_plot)
        color = colors[mod1(ch_idx, length(colors))]
        lines!(ax, t, epoch[!, ch], color = color, linewidth = 1)
    end
    hlines!(ax, [0.0], color = :black, linewidth = 0.5, linestyle = :dash)
end

function Base.show(io::IO, state::EpochRejectionState)
    rejected_count = sum(state.rejected)
    kept_count = state.n_total_epochs - rejected_count
    
    println(io, "EpochRejectionState:")
    println(io, "  Total epochs: $(state.n_total_epochs)")
    println(io, "  Rejected epochs: $rejected_count ($(round(100 * rejected_count / state.n_total_epochs, digits=1))%)")
    println(io, "  Kept epochs: $kept_count ($(round(100 * kept_count / state.n_total_epochs, digits=1))%)")
    println(io, "  Current page: $(state.current_page[])/$(state.n_pages)")
    println(io, "  Epochs per page: $(state.epochs_per_page)")
    println(io, "  Channels displayed: $(length(state.selected_channels))")
end


#=============================================================================
    INPUT VALIDATION
=============================================================================#

"""
Validate inputs for the rejection GUI.
"""
function _validate_rejection_gui_inputs(dat::EpochData, grid_size::Tuple{Int,Int})
    if isempty(dat.data)
        @minimal_error_throw("Cannot create rejection interface for empty EpochData")
    end
    if grid_size[1] * grid_size[2] <= 1
        @minimal_error_throw("grid_size must be positive, got $grid_size")
    end
end


#=============================================================================
    HELPER FUNCTIONS
=============================================================================#

"""
    _get_bad_channels_for_epoch(artifact_info::Union{Nothing,EpochRejectionInfo}, epoch_idx::Int)::Vector{Symbol}

Get the list of bad channels for a specific epoch from artifact detection info.
"""
function _get_bad_channels_for_epoch(artifact_info::Union{Nothing,EpochRejectionInfo}, epoch_idx::Int)::Vector{Symbol}
    if isnothing(artifact_info)
        return Symbol[]
    end
    
    # Filter rejections for this specific epoch
    bad_channels = Symbol[]
    for rejection in artifact_info.rejected_epochs
        if rejection.epoch == epoch_idx
            push!(bad_channels, rejection.label)
        end
    end
    
    return unique(bad_channels)
end


#=============================================================================
    HELPER FUNCTIONS FOR EXTRACTING RESULTS
=============================================================================#

"""
    get_rejected_epochs(state::EpochRejectionState)::Vector{Int}

Get indices of rejected epochs from the rejection state.

# Examples
```julia
state = reject_epochs_interactive(epochs)
# ... after review ...
rejected_indices = get_rejected_epochs(state)
```
"""
function get_rejected_epochs(state::EpochRejectionState)::Vector{Int}
    return findall(state.rejected)
end


"""
    save_rejection_decisions(state::EpochRejectionState, filename::String)

Save rejection decisions to a text file.

# Arguments
- `state::EpochRejectionState`: State from interactive rejection
- `filename::String`: Output filename

# Examples
```julia
state = reject_epochs_interactive(epochs)
# ... after review ...
save_rejection_decisions(state, "participant_1_manual_rejection.txt")
```
"""
function save_rejection_decisions(state::EpochRejectionState, filename::String)
    rejected_indices = get_rejected_epochs(state)
    
    open(filename, "w") do io
        println(io, "="^70)
        println(io, "MANUAL EPOCH REJECTION DECISIONS")
        println(io, "="^70)
        println(io, "")
        println(io, "Total epochs: $(state.n_total_epochs)")
        println(io, "Rejected epochs: $(length(rejected_indices))")
        println(io, "Kept epochs: $(state.n_total_epochs - length(rejected_indices))")
        println(io, "Rejection rate: $(round(100 * length(rejected_indices) / state.n_total_epochs, digits=1))%")
        println(io, "")
        println(io, "REJECTED EPOCH INDICES:")
        println(io, "-"^70)
        if isempty(rejected_indices)
            println(io, "None")
        else
            println(io, rejected_indices)
        end
        println(io, "")
        println(io, "KEPT EPOCH INDICES:")
        println(io, "-"^70)
        kept_indices = findall(.!state.rejected)
        println(io, kept_indices)
        println(io, "")
        println(io, "="^70)
    end
    
    @info "Rejection decisions saved to: $filename"
end

