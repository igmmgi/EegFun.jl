"""
    plot_epochs(dat::EpochData, channel_predicate::Function = channels; kwargs=Dict())

Plot epoched EEG data for specified channels.

# Arguments
- `dat::EpochData`: The epoched EEG data to plot.
- `channel_predicate::Function`: Function that returns boolean vector for channel filtering (default: channels - all channels)
- `kwargs`: Additional keyword arguments for customization.

# Returns
- `fig::Figure`: The Figure object containing the plot.
- `ax::Axis`: The Axis object containing the plot.

# Examples
```julia
# Plot all channels
plot_epochs(dat)

# Plot specific channels
plot_epochs(dat, channels([:Fp1, :Fp2]))

# Exclude reference channels
plot_epochs(dat, channels_not([:M1, :M2]))

# Plot frontal channels only
plot_epochs(dat, channels(1:10))

# Custom predicate
plot_epochs(dat, x -> startswith.(string.(x), "F"))

# Plot multiple channels separately
plot_epochs(dat, channels([:Fp1, :Fp2, :Fpz, :C1]), average_channels=false)
```

# Keyword Arguments
- `average_channels::Bool`: Whether to average across channels (default: true)
- `xlim`: X-axis limits
- `ylim`: Y-axis limits
- `title`: Plot title
- `xlabel`: X-axis label (default: "Time (S)")
- `ylabel`: Y-axis label (default: "mV")
- `linewidth`: Line width(s) (default: [1, 3])
- `color`: Line color(s) (default: [:grey, :black])
- `yreversed::Bool`: Whether to reverse Y-axis (default: false)
- `layout`: Subplot layout for multiple channels (default: auto-calculated)
"""
function plot_epochs(dat::EpochData, channel_predicate::Function = channels; kwargs=Dict())
    # Get the channels using the predicate
    selected_channels = channel_predicate(dat.layout.label)

    # Validate inputs
    isempty(selected_channels) && throw(ArgumentError("At least one channel must be specified"))
    invalid_channels = setdiff(selected_channels, dat.layout.label)
    !isempty(invalid_channels) && throw(ArgumentError("Invalid channels: $(join(invalid_channels, ", "))"))

    # Default keyword arguments
    default_kwargs = Dict(
        :average_channels => true,
        :xlim => nothing,
        :ylim => nothing,
        :title => nothing,
        :xlabel => "Time (S)",
        :ylabel => "mV",
        :linewidth => [1, 3],
        :color => [:grey, :black],
        :yreversed => false,
        :layout => nothing,
    )
    kwargs = merge(default_kwargs, kwargs)

    fig = Figure()
    
    if kwargs[:average_channels]

        # Single plot averaging across channels
        ax = Axis(fig[1, 1])
        _plot_epochs!(ax, dat, selected_channels, kwargs)

        # Set axis properties
        if length(selected_channels) == 1
            _set_axis_properties!(ax, kwargs, "$(selected_channels[1])")
        else
            _set_axis_properties!(ax, kwargs, "Avg: $(_print_vector(selected_channels, max_length = 8, n_ends = 3))")
        end

    else

        # Separate subplot for each channel
        n_channels = length(selected_channels)
        rows, cols = isnothing(kwargs[:layout]) ? best_rect(n_channels) : kwargs[:layout]
        
        # Calculate global y-range if ylim is not provided
        ylim = kwargs[:ylim]
        if isnothing(ylim)
            ylim = _calculate_global_yrange(dat, selected_channels)
        end
        
        for (idx, channel) in enumerate(selected_channels)
            row = fld(idx-1, cols) + 1
            col = mod(idx-1, cols) + 1
            ax = Axis(fig[row, col])
            _plot_epochs!(ax, dat, [channel], kwargs)
            
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
    fontsize_theme = Theme(fontsize=24)
    update_theme!(fontsize_theme)

    display(fig)
    return fig, ax

end

# Backward compatibility - keep the old method for existing code
function plot_epochs(dat::EpochData, channels::Vector{Symbol}; kwargs=Dict())
    channel_predicate = x -> x .∈ Ref(channels)
    return plot_epochs(dat, channel_predicate; kwargs = kwargs)
end

"""
    _plot_epochs!(ax, dat, channels, kwargs)

Internal function to plot epoched data for specified channels on a given axis.

# Arguments
- `ax::Axis`: The axis to plot on.
- `dat::EpochData`: The epoched EEG data.
- `channels::Vector{Symbol}`: List of channels to include.
- `kwargs::Dict`: Keyword arguments for customization.
"""
function _plot_epochs!(ax, dat, channels, kwargs)
    avg_data = zeros(nrow(dat.data[1]))
    
    for trial in eachindex(dat.data)
        trial_data = colmeans(dat.data[trial], channels)
        avg_data .+= trial_data
        lines!(ax, dat.data[trial][!, :time], trial_data, 
               color=kwargs[:color][1], linewidth=kwargs[:linewidth][1])
    end
    
    avg_data ./= length(dat.data)
    lines!(ax, dat.data[1][!, :time], avg_data, color=kwargs[:color][2], linewidth=kwargs[:linewidth][2])

end

"""
    _set_axis_properties!(ax, kwargs, default_title)

Internal function to set axis properties.

# Arguments
- `ax::Axis`: The axis to configure.
- `kwargs::Dict`: Keyword arguments for customization.
- `default_title::String`: Default title if no custom title is provided.
"""
function _set_axis_properties!(ax, kwargs, default_title)
    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    ax.title = isnothing(kwargs[:title]) ? default_title : kwargs[:title]
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    ax.yreversed = kwargs[:yreversed]
end

"""
    plot_epochs(dat::EpochData, channel::Symbol; kwargs...)

Convenience method for plotting a single channel.

# Arguments
- `dat::EpochData`: The epoched EEG data to plot.
- `channel::Symbol`: The channel to plot.
- `kwargs...`: Additional keyword arguments for customization.

# Returns
- `fig::Figure`: The Figure object containing the plot.
- `ax::Axis`: The Axis object containing the plot.

# Example
```julia
plot_epochs(epoch, :Fp1)
```
"""
plot_epochs(dat::EpochData, channel::Symbol; kwargs...) = plot_epochs(dat, [channel]; kwargs...)

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
function _calculate_global_yrange(dat::EpochData, channels::Vector{Symbol}; buffer::Float64 = 0.1)

    min_val, max_val = Inf, -Inf
    for channel in channels
        for trial in eachindex(dat.data)
            trial_data = dat.data[trial][!, channel]
            min_val = min(min_val, minimum(trial_data))
            max_val = max(max_val, maximum(trial_data))
        end
    end
    
    # Add a small buffer to the range
    buffer = buffer * (max_val - min_val)
    return (min_val - buffer, max_val + buffer)

end

"""
    plot_epochs_table(epochs::Vector{EpochData})

Create a table showing the number of epochs for each condition.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData objects, one for each condition

# Returns
- `fig::Figure`: The Figure object containing the table
"""
function plot_epochs_table(epochs::Vector{EpochData})
    # Create figure
    fig = Figure()
    
    # Create table data
    conditions = [1, 4, 5, 3]  # The conditions used in your code
    n_epochs = [length(epoch.data) for epoch in epochs]
    
    # Create table
    table = Table(fig[1, 1])
    
    # Add headers
    header = ["Condition", "Number of Epochs"]
    for (i, h) in enumerate(header)
        Label(table[1, i], h, fontsize=16, font=:bold)
    end
    
    # Add data rows
    for (i, (cond, n)) in enumerate(zip(conditions, n_epochs))
        Label(table[i+1, 1], string(cond), fontsize=14)
        Label(table[i+1, 2], string(n), fontsize=14)
    end
    
    # Adjust table properties
    table.halign = :center
    table.valign = :center
    
    return fig
end

"""
    epochs_to_dataframe(epochs::Vector{EpochData})

Create a DataFrame showing the number of epochs for each condition.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData objects, one for each condition

# Returns
- `DataFrame`: DataFrame with condition numbers and epoch counts
"""
function epochs_to_dataframe(epochs::Vector{EpochData})
    conditions = [1, 4, 5, 3]  # The conditions used in your code
    n_epochs = [length(epoch.data) for epoch in epochs]
    
    return DataFrame(
        condition = conditions,
        n_epochs = n_epochs
    )
end
