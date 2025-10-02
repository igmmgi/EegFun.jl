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
- `epochs_per_page::Int`: Number of epochs to show per page (default: 12)
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

# Notes
- Epochs are displayed with all selected channels overlaid
- Checked boxes indicate epochs to REJECT
- Unchecked boxes indicate epochs to KEEP
- Changes are saved in the state object
- Close the window when done reviewing
"""
function reject_epochs_interactive(
    dat::EpochData;
    channel_selection::Function = channels(),
    epochs_per_page::Int = 12,
    grid_size::Tuple{Int,Int} = (3, 4),
)::EpochRejectionState

    @info "Starting interactive epoch rejection interface"
    
    # Validate inputs
    _validate_rejection_gui_inputs(dat, epochs_per_page, grid_size)
    
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for display")
    end
    
    n_total_epochs = length(dat.data)
    n_pages = ceil(Int, n_total_epochs / epochs_per_page)
    
    @info "Displaying $(length(selected_channels)) channels"
    @info "Total epochs: $n_total_epochs across $n_pages pages"
    
    # Initialize rejection state (all epochs start as not rejected)
    rejected = fill(false, n_total_epochs)
    
    # Create figure
    fig = Figure(size = (1400, 900))
    
    # Create state object
    current_page = Observable(1)
    state = EpochRejectionState(
        dat,
        selected_channels,
        rejected,
        current_page,
        epochs_per_page,
        n_total_epochs,
        n_pages,
        Toggle[],
        Axis[],
        fig
    )
    
    # Create UI layout
    _create_rejection_interface!(fig, state, grid_size)
    
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
function _create_rejection_interface!(fig::Figure, state::EpochRejectionState, grid_size::Tuple{Int,Int})
    rows, cols = grid_size
    
    # Title
    Label(fig[1, 1:cols], "Interactive Epoch Rejection", 
          fontsize = 24, font = :bold, halign = :center)
    
    # Info label (updates with page info)
    info_label = Label(fig[2, 1:cols], "", fontsize = 18, halign = :center)
    
    # Update info label when page changes
    on(state.current_page) do page
        rejected_count = sum(state.rejected)
        info_label.text[] = "Page $page/$(state.n_pages) | " *
                           "Rejected: $rejected_count/$(state.n_total_epochs) epochs | " *
                           "✓ = REJECT, ✗ = KEEP"
    end
    notify(state.current_page)  # Trigger initial update
    
    # Create grid of epoch plots with checkboxes
    for row in 1:rows
        for col in 1:cols
            epoch_idx_in_page = (row - 1) * cols + col
            
            # Calculate grid position (epochs start at row 3)
            grid_row = 3 + (row - 1) * 2  # 2 rows per epoch (plot + checkbox)
            grid_col = col
            
            # Create axis for epoch plot
            ax = Axis(fig[grid_row, grid_col],
                     xlabel = "Time (s)",
                     ylabel = "μV",
                     xlabelsize = 12,
                     ylabelsize = 12,
                     xticklabelsize = 10,
                     yticklabelsize = 10)
            
            push!(state.epoch_axes, ax)
            
            # Create checkbox below the plot
            checkbox = Toggle(fig[grid_row + 1, grid_col], active = false)
            push!(state.checkboxes, checkbox)
            
            # Checkbox label
            checkbox_label = Label(fig[grid_row + 1, grid_col], "",
                                  fontsize = 14, halign = :center)
            
            # Setup checkbox observer to update rejection state
            on(checkbox.active) do is_checked
                # Calculate actual epoch index
                page = state.current_page[]
                epoch_idx = (page - 1) * state.epochs_per_page + epoch_idx_in_page
                
                if epoch_idx <= state.n_total_epochs
                    state.rejected[epoch_idx] = is_checked
                    
                    # Update label text
                    if is_checked
                        checkbox_label.text[] = "Epoch $epoch_idx: REJECT ✗"
                        checkbox_label.color[] = :red
                    else
                        checkbox_label.text[] = "Epoch $epoch_idx: Keep ✓"
                        checkbox_label.color[] = :green
                    end
                    
                    # Update info label
                    rejected_count = sum(state.rejected)
                    info_label.text[] = "Page $(state.current_page[])/$(state.n_pages) | " *
                                       "Rejected: $rejected_count/$(state.n_total_epochs) epochs | " *
                                       "✓ = REJECT, ✗ = KEEP"
                end
            end
        end
    end
    
    # Navigation buttons at the bottom
    nav_row = 3 + rows * 2 + 1
    
    # First button
    btn_first = Button(fig[nav_row, 1], label = "|◀ First", fontsize = 16)
    on(btn_first.clicks) do _
        state.current_page[] = 1
        _update_epoch_display!(state)
    end
    
    # Previous button
    btn_prev = Button(fig[nav_row, 2], label = "◀ Previous", fontsize = 16)
    on(btn_prev.clicks) do _
        if state.current_page[] > 1
            state.current_page[] -= 1
            _update_epoch_display!(state)
        end
    end
    
    # Page indicator
    Label(fig[nav_row, 3:4], "", fontsize = 16, halign = :center)
    
    # Next button
    btn_next = Button(fig[nav_row, 5], label = "Next ▶", fontsize = 16)
    on(btn_next.clicks) do _
        if state.current_page[] < state.n_pages
            state.current_page[] += 1
            _update_epoch_display!(state)
        end
    end
    
    # Last button
    btn_last = Button(fig[nav_row, 6], label = "Last ▶|", fontsize = 16)
    on(btn_last.clicks) do _
        state.current_page[] = state.n_pages
        _update_epoch_display!(state)
    end
    
    # Initial display
    _update_epoch_display!(state)
end


"""
Update the display when page changes.
"""
function _update_epoch_display!(state::EpochRejectionState)
    page = state.current_page[]
    start_idx = (page - 1) * state.epochs_per_page + 1
    end_idx = min(start_idx + state.epochs_per_page - 1, state.n_total_epochs)
    
    for (i, ax) in enumerate(state.epoch_axes)
        epoch_idx = start_idx + i - 1
        
        # Clear previous content
        empty!(ax)
        
        if epoch_idx <= state.n_total_epochs
            # Plot this epoch
            _plot_single_epoch!(ax, state, epoch_idx)
            
            # Update checkbox state
            state.checkboxes[i].active[] = state.rejected[epoch_idx]
            
            # Set title
            ax.title = "Epoch $epoch_idx"
            ax.titlesize = 14
        else
            # Hide empty slots
            ax.title = ""
            hidedecorations!(ax)
            hidespines!(ax)
            state.checkboxes[i].active[] = false
        end
    end
end


"""
Plot a single epoch in the given axis.
"""
function _plot_single_epoch!(ax::Axis, state::EpochRejectionState, epoch_idx::Int)
    epoch = state.epoch_data.data[epoch_idx]
    
    # Plot each channel
    colors = Makie.wong_colors()
    for (ch_idx, ch) in enumerate(state.selected_channels)
        color = colors[mod1(ch_idx, length(colors))]
        lines!(ax, epoch.time, epoch[!, ch], 
              color = color, 
              linewidth = 1,
              label = string(ch))
    end
    
    # Add zero lines
    hlines!(ax, [0.0], color = :black, linewidth = 0.5, linestyle = :dash)
    vlines!(ax, [0.0], color = :black, linewidth = 0.5, linestyle = :dash)
    
    # Only show legend if few channels
    if length(state.selected_channels) <= 5
        axislegend(ax, position = :rt, labelsize = 10)
    end
end


"""
Validate inputs for the rejection GUI.
"""
function _validate_rejection_gui_inputs(dat::EpochData, epochs_per_page::Int, grid_size::Tuple{Int,Int})
    if isempty(dat.data)
        @minimal_error_throw("Cannot create rejection interface for empty EpochData")
    end
    
    if epochs_per_page < 1
        @minimal_error_throw("epochs_per_page must be positive, got $epochs_per_page")
    end
    
    rows, cols = grid_size
    if rows * cols != epochs_per_page
        @minimal_error_throw("grid_size ($rows × $cols = $(rows*cols)) must equal epochs_per_page ($epochs_per_page)")
    end
    
    if rows < 1 || cols < 1
        @minimal_error_throw("grid_size must have positive dimensions, got $grid_size")
    end
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
    get_clean_epochs(state::EpochRejectionState)::EpochData

Create a new EpochData object with rejected epochs removed.

# Examples
```julia
state = reject_epochs_interactive(epochs)
# ... after review ...
clean_epochs = get_clean_epochs(state)
save("epochs_cleaned.jld2", "epochs", clean_epochs)
```
"""
function get_clean_epochs(state::EpochRejectionState)::EpochData
    kept_epochs = state.epoch_data.data[.!state.rejected]
    
    if isempty(kept_epochs)
        @minimal_warning "All epochs were rejected! Returning empty EpochData."
    end
    
    return EpochData(
        kept_epochs,
        state.epoch_data.layout,
        state.epoch_data.sample_rate,
        state.epoch_data.analysis_info
    )
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


"""
    Base.show(io::IO, state::EpochRejectionState)

Display epoch rejection state information.
"""
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

