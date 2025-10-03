"""
    AutoRejectParams

Parameters for automated artifact rejection

# Fields
- `n_interpolate::Vector{Int}`: Numbers of channels to interpolate to try
- `consensus::Vector{Float64}`: Consensus percentages to try
- `thresh_method::Symbol`: Method for threshold calculation (:mad or :std)
- `thresh_range::Vector{Float64}`: Range of thresholds to try
"""
struct AutoRejectParams
    n_interpolate::Vector{Int}
    consensus::Vector{Float64}
    thresh_method::Symbol
    thresh_range::Vector{Float64}
end

function AutoRejectParams(;
    n_interpolate = 1:4,
    consensus = 0.1:0.1:1.0,
    thresh_method = :mad,
    thresh_range = 1.5:0.5:4.0,
)
    AutoRejectParams(n_interpolate, consensus, thresh_method, thresh_range)
end

"""
    find_bad_channels(epochs::EpochData, params::AutoRejectParams)

Find optimal parameters and identify bad channels using cross-validation.

Returns:
- Dictionary of bad channels per epoch
- Optimal threshold parameters
"""
function find_bad_channels(epochs::EpochData, params::AutoRejectParams)
    n_epochs = length(epochs.data)
    n_channels = length(epochs.layout.data.label)

    # Pre-allocate results
    cv_scores = zeros(length(params.n_interpolate), length(params.consensus), length(params.thresh_range))

    # Calculate peak-to-peak amplitude for each channel/epoch
    p2p = zeros(n_channels, n_epochs)
    for (e_idx, epoch) in enumerate(epochs.data)
        for (c_idx, chan) in enumerate(epochs.layout.data.label)
            p2p[c_idx, e_idx] = maximum(epoch[!, chan]) - minimum(epoch[!, chan])
        end
    end

    # Cross-validation to find optimal parameters
    for (i_idx, n_interp) in enumerate(params.n_interpolate)
        for (c_idx, cons) in enumerate(params.consensus)
            for (t_idx, thresh) in enumerate(params.thresh_range)
                # Calculate threshold for each channel
                if params.thresh_method == :mad
                    channel_thresh = median(p2p, dims = 2) .+ thresh .* mads(p2p, dims = 2)
                else # :std
                    channel_thresh = mean(p2p, dims = 2) .+ thresh .* std(p2p, dims = 2)
                end
                # Find violations
                violations = p2p .> channel_thresh
                # Calculate consensus score
                consensus_score = mean(violations, dims = 1)[:]
                bad_epochs = consensus_score .> cons
                # Score this parameter combination
                cv_scores[i_idx, c_idx, t_idx] = mean(bad_epochs)
            end
        end
    end

    # Find optimal parameters
    best_params = argmin(cv_scores)
    opt_n_interp = params.n_interpolate[best_params[1]]
    opt_consensus = params.consensus[best_params[2]]
    opt_thresh = params.thresh_range[best_params[3]]

    # Find bad channels using optimal parameters
    bad_channels = Dict{Int,Vector{Symbol}}()
    for epoch_idx = 1:n_epochs
        epoch_p2p = p2p[:, epoch_idx]
        if params.thresh_method == :mad
            thresh = median(epoch_p2p) + opt_thresh * mad(epoch_p2p)
        else
            thresh = mean(epoch_p2p) + opt_thresh * std(epoch_p2p)
        end
        bad_idx = findall(epoch_p2p .> thresh)
        bad_channels[epoch_idx] = epochs.layout.data.label[bad_idx]
    end

    return bad_channels, (opt_n_interp, opt_consensus, opt_thresh)
end

"""
    repair_epochs!(epochs::EpochData, bad_channels::Dict, neighbours::Dict)

Repair epochs by interpolating bad channels using neighbor data.
"""
function repair_epochs!(epochs::EpochData, bad_channels::Dict, neighbours::Dict)
    for (epoch_idx, bad_chans) in bad_channels
        # Skip if no bad channels
        isempty(bad_chans) && continue

        # Get current epoch
        epoch = epochs.data[epoch_idx]

        # Interpolate each bad channel
        for bad_chan in bad_chans
            # Get neighbor data
            neighbor_labels = neighbours[Symbol(bad_chan)].electrodes
            good_neighbors = setdiff(neighbor_labels, bad_chans)

            # Skip if not enough good neighbors
            # length(good_neighbors) < 3 && continue

            # Get weights based on distance
            weights =
                [1/neighbours[Symbol(bad_chan)].distances[findfirst(==(n), neighbor_labels)] for n in good_neighbors]
            weights ./= sum(weights)

            # Interpolate
            interpolated = zeros(nrow(epoch))
            for (w, n) in zip(weights, good_neighbors)
                interpolated .+= w .* epoch[!, n]
            end

            # Update data
            epoch[!, bad_chan] = interpolated
        end
    end
end

"""
    plot_interpolation_comparison(epochs::EpochData, bad_channels::Dict, epoch_idx::Int)

Visualize the effect of different n_interpolate values.

Returns a figure showing:
1. Original data
2. Bad channels marked
3. Interpolated result
"""
function plot_interpolation_comparison(epochs::EpochData, bad_channels::Dict, epoch_idx::Int)
    fig = Figure()

    # Create grid layout and make it fill the figure
    gl = fig[1, 1] = GridLayout()
    colsize!(gl, 1, Relative(1.0))

    # Create observable for selected channel label
    selected_label = Observable("")
    label_text = gl[1, 1] = Label(fig, lift(x -> "Channel: $x", selected_label))

    # Get data range for consistent y-limits
    epoch_data = epochs.data[epoch_idx]
    all_data = Matrix(epoch_data[!, epochs.layout.data.label])
    y_min, y_max = extrema(all_data)
    y_range = y_max - y_min
    limits = (y_min - 0.1*y_range, y_max + 0.1*y_range)

    # Create observable for selected channel
    selected_channel = Observable("")

    # Original data - all lightgrey
    ax1 = Axis(gl[2, 1], title = "Original", limits = (nothing, limits))
    for chan in epochs.layout.data.label
        lines!(
            ax1,
            epoch_data.time,
            epoch_data[!, chan],
            color = :lightgrey,
            linewidth = 1,
            alpha = 0.5,
            inspectable = true,
        )
    end

    # Mark bad channels - bad in red, others lightgrey
    ax2 = Axis(gl[3, 1], title = "Bad Channels Marked", limits = (nothing, limits))
    for chan in epochs.layout.data.label
        color = chan in bad_channels[epoch_idx] ? :red : :lightgrey
        lines!(ax2, epoch_data.time, epoch_data[!, chan], color = color, linewidth = 1, alpha = 0.5, inspectable = true)
    end

    # After interpolation
    ax3 = Axis(gl[4, 1], title = "After Interpolation", limits = (nothing, limits))
    epochs_copy = copy(epochs)
    get_layout_neighbours_xyz!(epochs.layout, 1.5)

    for chan in epochs.layout.data.label
        color = chan in bad_channels[epoch_idx] ? :red : :lightgrey
        lines!(
            ax3,
            epochs_copy.data[epoch_idx].time,
            epochs_copy.data[epoch_idx][!, chan],
            color = color,
            linewidth = 1,
            alpha = 0.5,
            inspectable = true,
        )
    end

    # Add click interaction to highlight selected channel
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            plt = pick(fig)
            if plt !== nothing && plt[1] isa Lines
                # Get the index of the clicked line
                for (i, chan) in enumerate(epochs.layout.data.label)
                    if plt[1] == ax1.scene.plots[i] || plt[1] == ax2.scene.plots[i] || plt[1] == ax3.scene.plots[i]
                        selected_channel[] = string(chan)
                        selected_label[] = string(chan)
                        break
                    end
                end

                # Top plot - all lightgrey, selected channel thicker
                empty!(ax1)
                for chan in epochs.layout.data.label
                    linewidth = chan == selected_channel[] ? 2 : 1
                    lines!(
                        ax1,
                        epoch_data.time,
                        epoch_data[!, chan],
                        color = :lightgrey,
                        linewidth = linewidth,
                        alpha = 0.5,
                        inspectable = true,
                    )
                end

                # Middle and bottom plots - show bad channels in red
                for ax in [ax2, ax3]
                    data = ax == ax2 ? epoch_data : epochs_copy.data[epoch_idx]
                    empty!(ax)
                    for chan in epochs.layout.data.label
                        println(selected_channel[])
                        color = if chan == selected_channel[]
                            :black
                        elseif chan in bad_channels[epoch_idx]
                            :red
                        else
                            :lightgrey
                        end
                        linewidth = chan == selected_channel[] ? 2 : 1
                        lines!(
                            ax,
                            data.time,
                            data[!, chan],
                            color = color,
                            linewidth = linewidth,
                            alpha = 0.5,
                            inspectable = true,
                        )
                    end
                end
            end
        end
    end

    # Set row sizes to give more space to plots
    rowsize!(gl, 1, 30)  # Label row
    rowsize!(gl, 2, Relative(0.3)) # Plot rows
    rowsize!(gl, 3, Relative(0.3))
    rowsize!(gl, 4, Relative(0.3))

    rowgap!(gl, 10)  # Add gap between rows

    # Link all axes
    linkaxes!(ax1, ax2, ax3)

    return fig
end

"""
    autoreject(epochs::EpochData; params::AutoRejectParams = AutoRejectParams())

Main function that:
1. Finds optimal rejection parameters
2. Identifies bad channels
3. Repairs epochs through interpolation

Returns cleaned EpochData
"""
function autoreject(epochs::EpochData; params::AutoRejectParams = AutoRejectParams())
    # Find bad channels
    bad_chans, opt_params = find_bad_channels(epochs, params)

    # Get neighbors for interpolation
    neighbours = get_electrode_neighbours_xyz(epochs.layout)

    # Repair epochs
    repair_epochs!(epochs, bad_chans, neighbours)

    return epochs
end

"""
    mad(x::AbstractArray; dims=nothing, normalize=true)

Calculate Median Absolute Deviation (MAD) along specified dimensions.

# Arguments
- `x::AbstractArray`: Input array
- `dims`: Dimension(s) to compute MAD over
- `normalize=true`: If true, multiply by 1.4826 for normal distribution consistency

Returns MAD statistic.
"""
function mads(x::AbstractArray; dims = nothing, normalize = true)
    if isnothing(dims)
        med = median(x)
        mad_val = median(abs.(x .- med))
    else
        med = median(x, dims = dims)
        mad_val = median(abs.(x .- med), dims = dims)
    end
    return normalize ? 1.4826 * mad_val : mad_val
end
