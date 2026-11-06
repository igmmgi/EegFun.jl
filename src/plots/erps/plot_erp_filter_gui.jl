"""
    plot_erp_filter_gui(; channel = :Cz, kwargs...)
    plot_erp_filter_gui(filepath::String; kwargs...)
    plot_erp_filter_gui(erp::Union{ErpData,ContinuousData}; channel = nothing, kwargs...)
    plot_erp_filter_gui(erps::AbstractVector{<:Union{ErpData,ContinuousData}}; channel = nothing, display_plot = true)

Launch an interactive GUI for exploring IIR/FIR low-pass and high-pass filter effects on ERP or continuous waveforms.

# Arguments
- `filepath::String`: Path to a `.jld2` file containing ERP data
- `erp` / `erps`: Single dataset or vector of `ErpData` / `ContinuousData`

# Keyword Arguments
- `channel::Union{Symbol,Nothing}`: Initial channel to display (default: `:Cz` for demo, first channel otherwise)
- `display_plot::Bool`: Display the interactive figure (default: `true`)

# Examples
```julia
plot_erp_filter_gui()
plot_erp_filter_gui(channel = :Cz)
plot_erp_filter_gui("participant_1_erps.jld2", channel = :Pz)
```
"""
function plot_erp_filter_gui(; channel::Union{Symbol,Nothing} = :Cz, kwargs...)
    data = read_data(example_path("data/julia/erps/example1_erps.jld2"))
    return plot_erp_filter_gui(data[1]; channel = channel, kwargs...)
end

function plot_erp_filter_gui(filepath::String; kwargs...)
    data = read_data(filepath)
    isnothing(data) && @minimal_error "No data found in file: $filepath"
    return plot_erp_filter_gui(data; kwargs...)
end

plot_erp_filter_gui(erp::Union{ErpData,ContinuousData}; channel::Union{Symbol,Nothing} = nothing, kwargs...) =
    plot_erp_filter_gui([erp]; channel = channel, kwargs...)

function _add_labeled_toggle!(grid::GridLayout, row::Int, label_text::String = "Enable")
    hbox = GridLayout(grid[row, 1])
    tog = Toggle(hbox[1, 1], active = false, halign = :left)
    Label(hbox[1, 2], label_text, halign = :left)
    colsize!(hbox, 1, Auto())
    colsize!(hbox, 2, Auto())
    return tog
end

function plot_erp_filter_gui(
    erp_vec::AbstractVector{<:Union{ErpData,ContinuousData}};
    channel::Union{Symbol,Nothing} = nothing,
    display_plot::Bool = true,
)
    first_erp = erp_vec[1]
    all_channels = setdiff(propertynames(first_erp.data), meta_labels(first_erp))
    isempty(all_channels) && @minimal_error "No channels found in ERP data"

    initial_channel = isnothing(channel) ? all_channels[1] : channel
    initial_channel in all_channels || @minimal_error "Channel $initial_channel not found in data. Available: $(all_channels)"

    nyquist = first_erp.sample_rate / 2

    selected_channel = Observable(initial_channel)
    lowpass_enabled = Observable(false)
    lowpass_cutoff = Observable(30.0)
    lowpass_order = Observable(2)
    highpass_enabled = Observable(false)
    highpass_cutoff = Observable(0.1)
    highpass_order = Observable(1)
    filter_method_obs = Observable("butterworth")
    filter_func_obs = Observable("filtfilt")
    mirror_enabled = Observable(false)

    _set_window_title("ERP Filter Tool")
    fig = _eegfun_figure(size = (1200, 700))

    main_grid = fig[1, 1] = GridLayout()
    controls_grid = main_grid[1, 1] = GridLayout(valign = :top, tellheight = false)
    plot_area_grid = main_grid[1, 2] = GridLayout()
    colsize!(main_grid, 1, Relative(0.2))
    colsize!(main_grid, 2, Relative(0.8))
    rowsize!(main_grid, 1, Auto())

    Label(controls_grid[1, 1], "Filter Controls", fontsize = 16, font = :bold, halign = :left)
    Label(controls_grid[2, 1], "Channel:", halign = :left)
    channel_menu = Menu(controls_grid[3, 1], options = zip(string.(all_channels), all_channels), default = string(initial_channel))

    Label(controls_grid[4, 1], "Filter Method:", halign = :left)
    filter_methods = ["Butterworth (IIR)" => "butterworth", "FIR (Hamming)" => "fir"]
    method_menu = Menu(controls_grid[5, 1], options = zip([p.first for p in filter_methods], filter_methods), default = "Butterworth (IIR)")

    Label(controls_grid[6, 1], "Filter Function:", halign = :left)
    filter_funcs =
        ["filtfilt (2-pass zero-phase)" => "filtfilt", "filt_zerophase (1-pass FIR)" => "filt_zerophase", "filt (1-pass causal)" => "filt"]
    func_menu =
        Menu(controls_grid[7, 1], options = zip([p.first for p in filter_funcs], filter_funcs), default = "filtfilt (2-pass zero-phase)")

    Label(controls_grid[8, 1], "Mirror (for edge artifacts):", halign = :left)
    mirror_toggle = _add_labeled_toggle!(controls_grid, 9)

    Label(controls_grid[10, 1], "Lowpass Filter:", fontsize = 14, font = :bold, halign = :left)
    lowpass_toggle = _add_labeled_toggle!(controls_grid, 11)
    Label(controls_grid[12, 1], "Cutoff Frequency:", halign = :left)
    lowpass_cutoff_slider = Slider(controls_grid[13, 1], range = 1.0:1.0:(nyquist*0.9), startvalue = 30.0)
    lowpass_cutoff_label = Label(controls_grid[14, 1], "30.0 Hz", halign = :left)
    lp_order_title = Label(controls_grid[15, 1], "Order:", halign = :left)
    lowpass_order_slider = Slider(controls_grid[16, 1], range = 1:10, startvalue = 2)
    lowpass_order_label = Label(controls_grid[17, 1], "2", halign = :left)

    Label(controls_grid[18, 1], "Highpass Filter:", fontsize = 14, font = :bold, halign = :left)
    highpass_toggle = _add_labeled_toggle!(controls_grid, 19)
    Label(controls_grid[20, 1], "Cutoff Frequency:", halign = :left)
    highpass_cutoff_slider = Slider(controls_grid[21, 1], range = 0.01:0.01:10.0, startvalue = 0.1)
    highpass_cutoff_label = Label(controls_grid[22, 1], "0.10 Hz", halign = :left)
    hp_order_title = Label(controls_grid[23, 1], "Order:", halign = :left)
    highpass_order_slider = Slider(controls_grid[24, 1], range = 1:10, startvalue = 1)
    highpass_order_label = Label(controls_grid[25, 1], "1", halign = :left)

    rowgap!(controls_grid, 8)

    n_conditions = length(erp_vec)
    nrows, ncols = _best_rect(n_conditions)
    axes = Vector{Axis}(undef, n_conditions)
    for i = 1:n_conditions
        row_idx = div(i - 1, ncols) + 1
        col_idx = mod(i - 1, ncols) + 1
        axes[i] = Axis(
            plot_area_grid[row_idx, col_idx],
            xlabel = row_idx == nrows ? "Time (s)" : "",
            ylabel = col_idx == 1 ? "μV" : "",
            title = hasproperty(erp_vec[i], :condition_name) ? erp_vec[i].condition_name : "Continuous Data",
        )
    end
    length(axes) > 1 && linkxaxes!(axes...)

    function _fir_info_str(cutoff_hz::Float64, ftype::Symbol, filt_func_sym::Symbol, eff_len::Int, fs::Real)::String
        n_taps = _effective_fir_taps(eff_len, fs, cutoff_hz, ftype, nothing, filt_func_sym)
        if filt_func_sym == :filt
            delay_ms = (n_taps - 1) / (2.0 * fs) * 1000.0
            return @sprintf("FIR: %d taps (+%.1f ms shift)", n_taps, delay_ms)
        elseif filt_func_sym == :filt_zerophase
            return @sprintf("FIR: %d taps (1-pass, 0 ms)", n_taps)
        else
            return @sprintf("FIR: %d taps (2-pass, 0 ms)", n_taps)
        end
    end

    function refresh_param_labels!()
        base_len = nrow(first_erp.data)
        eff_len = (mirror_enabled[] && base_len > 1) ? (3 * base_len - 2) : base_len
        if filter_method_obs[] == "butterworth"
            lp_order_title.text = "Order:"
            lp_order_title.color = :black
            hp_order_title.text = "Order:"
            hp_order_title.color = :black
            suffix = filter_func_obs[] == "filt_zerophase" ? " (2-pass IIR)" : ""
            lowpass_order_label.text = "$(lowpass_order[])$suffix"
            highpass_order_label.text = "$(highpass_order[])$suffix"
        else
            lp_order_title.text = "Order (IIR only):"
            lp_order_title.color = :gray
            hp_order_title.text = "Order (IIR only):"
            hp_order_title.color = :gray
            f_sym = Symbol(filter_func_obs[])
            lowpass_order_label.text = _fir_info_str(lowpass_cutoff[], :lp, f_sym, eff_len, first_erp.sample_rate)
            highpass_order_label.text = _fir_info_str(highpass_cutoff[], :hp, f_sym, eff_len, first_erp.sample_rate)
        end
    end

    function _draw_trace!(ax, dat, ch_sym, color)
        if dat isa ErpData
            kw = color === :black ? (;) : (; color = color)
            plot_erp!(fig, ax, [dat]; channel_selection = channels(ch_sym), legend = false, kw...)
        else
            t = (0:(size(dat.data, 1)-1)) ./ dat.sample_rate
            lines!(ax, t, dat.data[!, ch_sym], color = color)
        end
    end

    function update_plot!()
        refresh_param_labels!()
        ch_sym = selected_channel[]
        for (idx, ax) in enumerate(axes)
            empty!(ax)
            _draw_trace!(ax, erp_vec[idx], ch_sym, :black)

            if lowpass_enabled[] || highpass_enabled[]
                filtered_erp = copy(erp_vec[idx])
                mirror_enabled[] && mirror!(filtered_erp, :both)

                filt_method = filter_method_obs[] == "butterworth" ? :iir : :fir
                filt_func = Symbol(filter_func_obs[])
                filt_func_eff = (filt_method == :iir && filt_func == :filt_zerophase) ? :filtfilt : filt_func

                try
                    if lowpass_enabled[]
                        lowpass_filter!(
                            filtered_erp,
                            lowpass_cutoff[];
                            order = lowpass_order[],
                            filter_method = filt_method,
                            filter_func = filt_func_eff,
                        )
                    end
                    if highpass_enabled[]
                        highpass_filter!(
                            filtered_erp,
                            highpass_cutoff[];
                            order = highpass_order[],
                            filter_method = filt_method,
                            filter_func = filt_func_eff,
                        )
                    end
                    mirror_enabled[] && unmirror!(filtered_erp, :both)

                    _draw_trace!(ax, filtered_erp, ch_sym, :red)
                    line_plots = filter(p -> p isa Lines, ax.scene.plots)
                    if length(line_plots) >= 2
                        axislegend(ax, [line_plots[1], line_plots[2]], ["Original", "Filtered"], position = :rt, framevisible = false)
                    end
                catch e
                    if e isa BoundsError
                        text!(
                            ax,
                            0.5,
                            0.5,
                            text = "Filter error:\nFIR too long for data.\nUse IIR instead.",
                            align = (:center, :center),
                            fontsize = 10,
                            color = :red,
                        )
                    else
                        rethrow(e)
                    end
                end
            end
        end
    end

    update_plot!()

    on(ch -> (selected_channel[] = ch; update_plot!()), channel_menu.selection)
    on(p -> (filter_method_obs[] = p isa Pair ? p[2] : p; update_plot!()), method_menu.selection)
    on(p -> (filter_func_obs[] = p isa Pair ? p[2] : p; update_plot!()), func_menu.selection)
    on(val -> (mirror_enabled[] = val; update_plot!()), mirror_toggle.active)
    on(val -> (lowpass_enabled[] = val; update_plot!()), lowpass_toggle.active)
    on(val -> (lowpass_cutoff[] = val; lowpass_cutoff_label.text = @sprintf("%.1f Hz", val); update_plot!()), lowpass_cutoff_slider.value)
    on(val -> (lowpass_order[] = val; update_plot!()), lowpass_order_slider.value)
    on(val -> (highpass_enabled[] = val; update_plot!()), highpass_toggle.active)
    on(
        val -> (highpass_cutoff[] = val; highpass_cutoff_label.text = @sprintf("%.2f Hz", val); update_plot!()),
        highpass_cutoff_slider.value,
    )
    on(val -> (highpass_order[] = val; update_plot!()), highpass_order_slider.value)

    screen = display_plot ? _display_figure(fig) : nothing
    _set_window_title("ERP Filter Tool")
    return (fig = fig, axes = axes, screen = screen)
end
