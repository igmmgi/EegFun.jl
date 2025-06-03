"""
    plot_epochs(dat::EpochData, channels::Vector{Symbol}; kwargs...)

Plot epoched EEG data for specified channels. This function can either average across channels
and plot a single trace or create separate subplots for each channel.

# Arguments
- `dat::EpochData`: The epoched EEG data to plot.
- `channels::Vector{Symbol}`: List of channels to include in the plot.
- `kwargs...`: Additional keyword arguments for customization.

# Keyword Arguments
- `average_channels::Bool=true`: If true, averages across channels for a single plot. If false, creates separate subplots for each channel.
- `xlim::Union{Nothing, Tuple}`: X-axis limits (default: nothing).
- `ylim::Union{Nothing, Tuple}`: Y-axis limits (default: nothing).
- `title::Union{Nothing, String}`: Plot title (default: nothing).
- `xlabel::String`: X-axis label (default: "Time (S)").
- `ylabel::String`: Y-axis label (default: "mV").
- `linewidth::Vector{Int}`: Line widths for individual trials and average (default: [1, 3]).
- `color::Vector{Symbol}`: Colors for individual trials and average (default: [:grey, :black]).
- `yreversed::Bool`: Whether to reverse the Y-axis (default: false).
- `layout::Union{Nothing, Tuple{Int, Int}}`: Custom subplot layout (rows, cols) when `average_channels=false` (default: nothing).

# Returns
- `fig::Figure`: The Figure object containing the plot.
- `ax::Union{Axis, Vector{Axis}}`: The Axis object(s) containing the plot(s).

# Example
```julia
# Plot a single channel
plot_epochs(epoch, :Fp1)

# Plot multiple channels with averaging
plot_epochs(epoch, [:Fp1, :Fp2])

# Plot multiple channels separately
plot_epochs(epoch, [:Fp1, :Fp2, :Fpz, :C1], average_channels=false)
```
"""
function plot_epochs(dat::EpochData, channels::Vector{Symbol}; kwargs=Dict())

    # Validate inputs
    isempty(channels) && throw(ArgumentError("At least one channel must be specified"))
    invalid_channels = setdiff(channels, dat.layout.label)
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
        _plot_epochs!(ax, dat, channels, kwargs)

        # Set axis properties
        if length(channels) == 1
            _set_axis_properties!(ax, kwargs, "$(channels[1])")
        else
            _set_axis_properties!(ax, kwargs, "Avg: $(print_vector_(channels, max_length = 8, n_ends = 3))")
        end

    else

        # Separate subplot for each channel
        n_channels = length(channels)
        rows, cols = isnothing(kwargs[:layout]) ? best_rect(n_channels) : kwargs[:layout]
        
        # Calculate global y-range if ylim is not provided
        ylim = kwargs[:ylim]
        if isnothing(ylim)
            ylim = _calculate_global_yrange(dat, channels)
        end
        
        for (idx, channel) in enumerate(channels)
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
