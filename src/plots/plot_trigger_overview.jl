########################################################
# Trigger overview plotting functions
########################################################

"""
    _trigger_time_count(time, triggers)

Internal function to process trigger data and count occurrences.

# Arguments
- `time`: Vector of time points
- `triggers`: Vector of trigger values

# Returns
- `trigger_times`: Vector of times when triggers occurred
- `trigger_values`: Vector of trigger values at those times
- `trigger_count`: OrderedDict mapping trigger values to their counts
"""
function _trigger_time_count(time, triggers)

    # Since triggers are already cleaned (only onset values), just find non-zero values
    trigger_indices = findall(triggers .!= 0)

    if isempty(trigger_indices)
        return Float64[], Int[], OrderedDict{Int,Int}()
    end

    trigger_values = triggers[trigger_indices]
    trigger_times = time[trigger_indices]
    trigger_count = OrderedDict(i => 0 for i in sort!(collect(Set(trigger_values))))
    for val in trigger_values
        trigger_count[val] += 1
    end

    return trigger_times, trigger_values, trigger_count

end

"""
    plot_trigger_overview(trigger_times, trigger_values, trigger_count)

Plot trigger events as a scatter plot with vertical lines.

# Arguments
- `trigger_times`: Vector of times when triggers occurred
- `trigger_values`: Vector of trigger values at those times
- `trigger_count`: OrderedDict mapping trigger values to their counts

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = true)

    if isempty(trigger_count)
        @warn "No triggers found in the data"
        fig = Figure()
        ax = Axis(fig[1, 1])
        return fig, ax
    end

    fig = Figure()
    ax = Axis(fig[1, 1], yticks = (1:length(trigger_count.keys), string.(trigger_count.keys)))

    # Pre-compute trigger data for each type to avoid repeated filtering
    trigger_data = Dict{Int,Vector{Float64}}()
    for (key, _) in trigger_count
        trigger_data[key] = trigger_times[trigger_values .== key]
    end

    for (unique, (key, value)) in enumerate(trigger_count)
        times = trigger_data[key]
        y_pos = fill(unique, length(times))
        scatter!(ax, times, y_pos, label = "$key: $(string(value))", markersize = DEFAULT_MARKER_SIZE)
        # Add vertical lines
        for (t, y) in zip(times, y_pos)
            lines!(
                ax,
                [t, t],
                [y - DEFAULT_LINE_OFFSET, y + DEFAULT_LINE_OFFSET],
                color = :black,
                linewidth = DEFAULT_LINE_WIDTH_TRIGGER,
            )
        end
    end
    fig[1, 2] = Legend(fig, ax)
    ax.ylabel = "Trigger Value"
    ax.xlabel = "Time (S)"

    if display_plot
        display(fig)
    end

    return fig, ax

end

"""
    plot_trigger_overview(dat::BioSemiBDF.BioSemiData)

Plot trigger events from BioSemi BDF data.

# Arguments
- `dat`: BioSemiData object containing the EEG data

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_trigger_overview(dat::BioSemiBDF.BioSemiData; display_plot = true)
    @info "Plotting trigger (raw) overview for BioSemi data"
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.time, dat.triggers.raw)
    return plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = display_plot)
end

"""
    plot_trigger_overview(dat::ContinuousData)

Plot trigger events from ContinuousData object.

# Arguments
- `dat`: ContinuousData object containing the EEG data

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_trigger_overview(dat::ContinuousData; display_plot = true)
    @info "Plotting trigger (cleaned) overview for ContinuousData"
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.data.time, dat.data.triggers)
    return plot_trigger_overview(trigger_times, trigger_values, trigger_count; display_plot = display_plot)
end
