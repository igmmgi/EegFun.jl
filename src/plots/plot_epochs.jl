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
        :linewidth => [1, 2],
        :color => [:grey, :red],
        :yreversed => false,
        # layout can be: false (default auto grid), true (use data's layout semantics), or [rows, cols]
        :layout => false,
        :layout_plot_width => 0.12,
        :layout_plot_height => 0.12,
        :layout_margin => 0.02,
        :layout_show_scale => true,
        :layout_scale_position => [0.95, 0.05],
        :layout_scale_width => 0.14,
        :layout_scale_height => 0.14,
        :legend => true,
        :legend_label => "",
        :dims => nothing,
        :hidedecorations => false,
        :theme_fontsize => 24,
        :plot_avg_trials => true,                # draw ERP average overlay
    )
    kwargs = merge(default_kwargs, kwargs)

    fig = Figure()
    axes = Axis[]  # Keep track of all axes created

    # Precompute ERP once when average overlay requested
    erp_dat = nothing
    if kwargs[:plot_avg_trials]
        erp_dat = average_epochs(dat_subset)
    end

    if kwargs[:average_channels]

        # Single plot averaging across channels
        ax = Axis(fig[1, 1])
        push!(axes, ax)
        _plot_epochs!(ax, dat_subset, all_plot_channels, kwargs)
        kwargs[:plot_avg_trials] && _plot_epochs_from_erp!(ax, erp_dat, all_plot_channels, kwargs)

        # Set axis properties
        if length(all_plot_channels) == 1
            _set_axis_properties!(ax, kwargs, "$(all_plot_channels[1])")
        else
            _set_axis_properties!(ax, kwargs, "Avg: $(print_vector(all_plot_channels, max_length = 8, n_ends = 3))")
        end

    else

        # Separate subplot for each channel
        n_channels = length(all_plot_channels)
        # Handle layout kwarg: Bool or Vector{Int}
        layout_kw = kwargs[:layout]
        if layout_kw === true
            _plot_epochs_layout!(fig, axes, dat_subset, erp_dat, all_plot_channels, kwargs)
        else
            if layout_kw === false || isnothing(layout_kw)
                rows, cols = best_rect(n_channels)
            else
                # Validate custom grid dims
                if length(layout_kw) != 2 || any(x -> x <= 0, layout_kw)
                    throw(ArgumentError("layout must be either a Bool or a 2-element vector [rows, cols] with positive integers"))
                end
                rows, cols = layout_kw
                if rows * cols < n_channels
                    throw(ArgumentError("layout grid ($(rows)×$(cols)=$(rows*cols)) is too small for $n_channels channels"))
                end
            end

            # Calculate y-range if not provided (using shared yrange helper)
            ylim = kwargs[:ylim]
            if isnothing(ylim)
                yr = ylimits(dat_subset; channel_selection = channels(all_plot_channels))
                ylim = (yr[1], yr[2])
            end

            for (idx, channel) in enumerate(all_plot_channels)
                row = fld(idx-1, cols) + 1
                col = mod(idx-1, cols) + 1
                ax = Axis(fig[row, col])
                push!(axes, ax)
                _plot_epochs!(ax, dat_subset, [channel], kwargs)
                kwargs[:plot_avg_trials] && _plot_epochs_from_erp!(ax, erp_dat, [channel], kwargs)

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
    end

    # Link axes (both x and y) when multiple axes are present
    if length(axes) > 1
        Makie.linkaxes!(axes...)
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
    trial_buffer = zeros(Float64, n_samples(dat))

    # Cache time vector and styles
    time_vec = dat.data[1][!, :time]
    trial_color = kwargs[:color][1]
    trial_linewidth = kwargs[:linewidth][1]

    # Local helper to compute row-wise mean across selected channels without allocations
    @inline function _mean_over_channels!(dest::Vector{Float64}, df::DataFrame, cols::Vector{Symbol})
        n = length(dest)
        # Assume consistent sample length across trials; only resize if mismatched
        if n != nrow(df)
            resize!(dest, nrow(df))
            n = length(dest)
        end
        fill!(dest, 0.0)
        @inbounds for j in 1:length(cols)
            col = df[!, cols[j]]
            @inbounds @simd for i in 1:n
                dest[i] += col[i]
            end
        end
        invm = 1.0 / length(cols)
        @inbounds @simd for i in 1:n
            dest[i] *= invm
        end
        return dest
    end

    # Plot only trials (no average here)
    for trial_df in dat.data
        _mean_over_channels!(trial_buffer, trial_df, channels)
        lines!(ax, time_vec, copy(trial_buffer), color = trial_color, linewidth = trial_linewidth)
    end

    return nothing
end


function _plot_epochs_from_erp!(ax, erp_dat::ErpData, channels::Vector{Symbol}, kwargs)::Nothing
    time_vec = erp_dat.data[!, :time]
    avg_color = kwargs[:color][2]
    avg_linewidth = kwargs[:linewidth][2]

    if length(channels) == 1
        ch = channels[1]
        lines!(ax, time_vec, erp_dat.data[!, ch], color = avg_color, linewidth = avg_linewidth)
    else
        vals = similar(time_vec, Float64)
        fill!(vals, 0.0)
        @inbounds for ch in channels
            col = erp_dat.data[!, ch]
            @inbounds @simd for i in eachindex(vals)
                vals[i] += col[i]
            end
        end
        invm = 1.0 / length(channels)
        @inbounds @simd for i in eachindex(vals)
            vals[i] *= invm
        end
        lines!(ax, time_vec, vals, color = avg_color, linewidth = avg_linewidth)
    end
    return nothing
end

function _plot_epochs_layout!(fig::Figure, axes::Vector{Axis}, dat::EpochData, erp_dat, all_plot_channels::Vector{Symbol}, kwargs::Dict)
    # Ensure 2D coordinates exist
    if !all(in.([:x2, :y2], Ref(propertynames(dat.layout.data))))
        polar_to_cartesian_xy!(dat.layout)
    end

    # Determine global y-lims for consistency across small axes
    ylim = kwargs[:ylim]
    if isnothing(ylim)
        yr = ylimits(dat; channel_selection = channels(all_plot_channels))
        ylim = (yr[1], yr[2])
    end

    # Normalize positions to [0,1]
    x2 = dat.layout.data.x2
    y2 = dat.layout.data.y2
    minx, maxx = extrema(x2)
    miny, maxy = extrema(y2)
    xrange = maxx - minx
    yrange = maxy - miny
    xrange = xrange == 0 ? 1.0 : xrange
    yrange = yrange == 0 ? 1.0 : yrange

    plot_w = get(kwargs, :layout_plot_width, 0.12)
    plot_h = get(kwargs, :layout_plot_height, 0.12)
    margin = get(kwargs, :layout_margin, 0.02)

    # Map channel -> position
    pos_map = Dict{Symbol,Tuple{Float64,Float64}}()
    for (lab, x, y) in zip(dat.layout.data.label, x2, y2)
        nx = (x - minx) / xrange
        ny = (y - miny) / yrange
        pos_map[Symbol(lab)] = (nx, ny)
    end

    # Create axes at positions
    for ch in all_plot_channels
        pos = get(pos_map, ch, (0.5, 0.5))
        ax = Axis(
            fig[1, 1],
            width = Relative(plot_w),
            height = Relative(plot_h),
            halign = clamp(pos[1], margin, 1 - margin),
            valign = clamp(pos[2], margin, 1 - margin),
        )
        push!(axes, ax)
        _plot_epochs!(ax, dat, [ch], kwargs)
        erp_dat !== nothing && _plot_epochs_from_erp!(ax, erp_dat, [ch], kwargs)

        # Suppress axis labels on all but the final axis; set only limits and title for now
        axis_kwargs = merge(kwargs, Dict(:ylim => ylim, :xlabel => "", :ylabel => ""))
        _set_axis_properties!(ax, axis_kwargs, string(ch))
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
    end

    # Optional extra scale axis in bottom-right
    if get(kwargs, :layout_show_scale, true)
        sp = get(kwargs, :layout_scale_position, [0.95, 0.05])
        sw = get(kwargs, :layout_scale_width, 0.14)
        sh = get(kwargs, :layout_scale_height, 0.14)
        scale_ax = Axis(
            fig[1, 1],
            width = Relative(sw),
            height = Relative(sh),
            halign = sp[1],
            valign = sp[2],
        )
        push!(axes, scale_ax)
        # No data in this axis; just show labels and limits
        tmin, tmax = (dat.data[1].time[1], dat.data[1].time[end])
        axis_kwargs = merge(kwargs, Dict(:ylim => ylim, :xlim => (tmin, tmax)))
        _set_axis_properties!(scale_ax, axis_kwargs, "")
        scale_ax.xticklabelsvisible = true
        scale_ax.yticklabelsvisible = true
    end
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
