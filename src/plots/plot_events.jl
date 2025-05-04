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

    trigger_indices = findall(diff(triggers) .>= 1)

    if isempty(trigger_indices)
        return Float64[], Int[], OrderedDict{Int,Int}()
    end

    trigger_values = triggers[trigger_indices.+1]
    trigger_times = time[trigger_indices.+1]
    trigger_count = OrderedDict(i => 0 for i in sort!(collect(Set(trigger_values))))
    for val in trigger_values
        trigger_count[val] += 1
    end

    return trigger_times, trigger_values, trigger_count

end

"""
    plot_events(trigger_times, trigger_values, trigger_count)

Plot trigger events as a scatter plot with vertical lines.

# Arguments
- `trigger_times`: Vector of times when triggers occurred
- `trigger_values`: Vector of trigger values at those times
- `trigger_count`: OrderedDict mapping trigger values to their counts

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_events(trigger_times, trigger_values, trigger_count)

    if isempty(trigger_count)
        @warn "No triggers found in the data"
        fig = Figure()
        ax = Axis(fig[1, 1])
        return fig, ax
    end

    fig = Figure()
    ax = Axis(fig[1, 1], yticks = (1:length(trigger_count.keys), string.(trigger_count.keys)))
    for (unique, (key, value)) in enumerate(trigger_count)
        times = trigger_times[trigger_values.==key]
        y_pos = repeat([unique], length(trigger_values[trigger_values.==key]))
        scatter!(ax, times, y_pos, label = "$key: $(string(value))", markersize = 15)
        # Add vertical lines
        for (t, y) in zip(times, y_pos)
            lines!(ax, [t, t], [y - 0.1, y + 0.1], color = :black, linewidth = 1)
        end
    end
    fig[1, 2] = Legend(fig, ax)
    ax.ylabel = "Trigger Value"
    ax.xlabel = "Time (S)"

    return fig, ax

end

"""
    plot_events(dat::BioSemiBDF.BioSemiData)

Plot trigger events from BioSemi BDF data.

# Arguments
- `dat`: BioSemiData object containing the EEG data

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_events(dat::BioSemiBDF.BioSemiData)
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.time, dat.triggers.raw)
    return plot_events(trigger_times, trigger_values, trigger_count)
end

"""
    plot_events(dat::ContinuousData)

Plot trigger events from ContinuousData object.

# Arguments
- `dat`: ContinuousData object containing the EEG data

# Returns
- `fig`: The Makie Figure object
- `ax`: The Axis object containing the plot
"""
function plot_events(dat::ContinuousData)
    trigger_times, trigger_values, trigger_count = _trigger_time_count(dat.data.time, dat.data.triggers)
    return plot_events(trigger_times, trigger_values, trigger_count)
end

# layout = read_layout("./layouts/biosemi72.csv");
# dat = read_bdf("../Flank_C_3.bdf");
# plot_events(dat)
# dat = create_eeg_dataframe(dat, layout);
# plot_events(dat)