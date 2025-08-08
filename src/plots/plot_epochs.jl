function plot_epochs(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(), 
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
    kwargs = Dict())::Tuple{Figure, Union{Axis, Vector{Axis}}}
    
    # Use subset to get the data we want to plot (same pattern as other functions)
    dat_subset = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra
    )
    
    # Get all non-metadata channels from the subset (it already contains only what we want)
    selected_channels = channel_labels(dat_subset)  # Gets EEG channels from layout
    extra_channels = extra_labels(dat_subset)       # Gets extra channels (EOG, etc.)
    all_plot_channels = vcat(selected_channels, extra_channels)
    
    # Validate we have channels to plot
    isempty(all_plot_channels) && throw(ArgumentError("No channels selected for plotting"))
    
    # Info about what we're plotting
    @info "plot_epochs: Plotting $(length(all_plot_channels)) channels across $(length(dat_subset.data)) epochs"

    # Default keyword arguments
    default_kwargs = Dict(
        :average_channels => false,
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => [1, 3],
        :color => [:grey, :black],
        :yreversed => false,
        :layout => nothing,
        :legend => true,
        :legend_label => "",
        :dims => nothing,
        :hidedecorations => false,
        :theme_fontsize => 24,
    )
    kwargs = merge(default_kwargs, kwargs)

    fig = Figure()
    axes = Axis[]  # Keep track of all axes created

    if kwargs[:average_channels]

        # Single plot averaging across channels
        ax = Axis(fig[1, 1])
        push!(axes, ax)
        _plot_epochs!(ax, dat_subset, all_plot_channels, kwargs)

        # Set axis properties
        if length(all_plot_channels) == 1
            _set_axis_properties!(ax, kwargs, "$(all_plot_channels[1])")
        else
            _set_axis_properties!(ax, kwargs, "Avg: $(print_vector(all_plot_channels, max_length = 8, n_ends = 3))")
        end

    else

        # Separate subplot for each channel
        n_channels = length(all_plot_channels)
        if isnothing(kwargs[:layout])
            rows, cols = best_rect(n_channels)
        else
            # Validate custom layout
            layout = kwargs[:layout]
            (length(layout) != 2 || any(x -> x <= 0, layout)) && 
                throw(ArgumentError("layout must be a 2-element vector of positive integers [rows, cols]"))
            rows, cols = layout
            (rows * cols < n_channels) && 
                throw(ArgumentError("layout grid ($(rows)×$(cols)=$(rows*cols)) is too small for $n_channels channels"))
        end

        # Calculate global y-range if ylim is not provided
        ylim = kwargs[:ylim]
        if isnothing(ylim)
            ylim = _calculate_global_yrange(dat_subset, all_plot_channels)
        end

        for (idx, channel) in enumerate(all_plot_channels)
            row = fld(idx-1, cols) + 1
            col = mod(idx-1, cols) + 1
            ax = Axis(fig[row, col])
            push!(axes, ax)
            _plot_epochs!(ax, dat_subset, [channel], kwargs)

            # Set axis properties with ylim
            axis_kwargs = merge(kwargs, Dict(:ylim => ylim))

            # Only add x and y labels to outer left column and bottom row
            if col != 1
                axis_kwargs = merge(axis_kwargs, Dict(:ylabel => ""))
                ax.yticklabelsvisible = false
            end
            if row != rows
                axis_kwargs = merge(axis_kwargs, Dict(:xlabel => ""))
                ax.xticklabelsvisible = false
            end

            _set_axis_properties!(ax, axis_kwargs, "$channel")

        end
    end

    # Theme adjustments
    fontsize_theme = Theme(fontsize = kwargs[:theme_fontsize])
    update_theme!(fontsize_theme)

    display(fig)
    
    # Return fig and all axes (or single axis if only one)
    return fig, length(axes) == 1 ? first(axes) : axes

end

function _plot_epochs!(ax, dat, channels, kwargs)::Nothing
    # Pre-allocate arrays for better performance
    n_timepoints = nrow(dat.data[1])
    n_trials = length(dat.data)
    avg_data = zeros(Float64, n_timepoints)
    
    # Cache time vector and colors for efficiency
    time_vec = dat.data[1][!, :time]
    trial_color = kwargs[:color][1]
    avg_color = kwargs[:color][2]
    trial_linewidth = kwargs[:linewidth][1]
    avg_linewidth = kwargs[:linewidth][2]

    # Process all trials efficiently
    for trial_idx in eachindex(dat.data)
        trial_data = colmeans(dat.data[trial_idx], channels)
        
        # Accumulate for average (more efficient than repeated allocation)
        avg_data .+= trial_data
        
        # Plot individual trial
        lines!(ax, time_vec, trial_data, color = trial_color, linewidth = trial_linewidth)
    end

    # Calculate and plot average
    avg_data ./= n_trials
    lines!(ax, time_vec, avg_data, color = avg_color, linewidth = avg_linewidth)
    
    return nothing
end


function _set_axis_properties!(ax::Axis, kwargs::Dict, default_title::String)::Nothing
    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    ax.title = isnothing(kwargs[:title]) ? default_title : kwargs[:title]
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]
    return nothing
end




"""
    _calculate_global_yrange(dat::EpochData, channels::Vector{Symbol})

Calculate the global y-range across all specified channels.

# Arguments
- `dat::EpochData`: The epoched EEG data.
- `channels::Vector{Symbol}`: List of channels to include.
- `buffer::Float64`: Buffer to add to the range (default: 0.1)

# Returns
- `Tuple{Float64, Float64}`: The minimum and maximum values across all channels.
"""
function _calculate_global_yrange(dat::EpochData, channels::Vector{Symbol}; buffer::Float64 = 0.1)::Tuple{Float64, Float64}
    min_val, max_val = Inf, -Inf
    
    for channel in channels
        for trial in eachindex(dat.data)
            trial_data = dat.data[trial][!, channel]
            min_val = min(min_val, minimum(trial_data))
            max_val = max(max_val, maximum(trial_data))
        end
    end

    # Handle edge case where all values are the same
    if min_val == max_val
        buffer_val = max(abs(min_val) * 0.1, 1.0)  # 10% of value or minimum 1.0
        return (min_val - buffer_val, max_val + buffer_val)
    end

    # Add a proportional buffer to the range
    buffer_val = buffer * (max_val - min_val)
    return (min_val - buffer_val, max_val + buffer_val)
end




