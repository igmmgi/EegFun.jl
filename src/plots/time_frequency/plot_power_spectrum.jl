
# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_POWER_SPECTRUM_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Spectral analysis parameters
    :window_size => (1024, "Size of the FFT window for spectral estimation"),
    :overlap => (0.5, "Overlap between windows for Welch's method (0.0 to 1.0)"),
    :max_freq => (200.0, "Maximum frequency to display in Hz"),
    :window_function => (DSP.hanning, "Window function for spectral estimation"),

    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    :show_legend => (true, "Whether to show the legend"),
    :show_freq_bands => (true, "Whether to show frequency band indicators"),

    # Axis styling
    :xlabel => ("Frequency (Hz)", "X-axis label"),
    :ylabel => ("Power Spectral Density (μV²/Hz)", "Y-axis label"),
    :title => ("Power Spectrum", "Plot title"),

    # Scale parameters
    :x_scale => (:linear, "X-axis scale: :linear or :log10"),
    :y_scale => (:linear, "Y-axis scale: :linear or :log10"),
    :unit => (:linear, "Power unit: :linear (μV²/Hz) or :dB (decibels)"),

    # Font sizes
    :title_fontsize => (16, "Font size for title"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :legend_fontsize => (12, "Font size for legend"),

    # Line styling
    :linewidth => (2, "Line width for spectrum lines"),
    :line_alpha => (0.8, "Transparency for spectrum lines"),

    # Frequency band styling
    :freq_band_alpha => (0.3, "Transparency for frequency band indicators"),
    :freq_band_height => (0.1, "Height of frequency band indicators"),

    # Grid styling
    :xgrid => (true, "Whether to show x-axis grid"),
    :ygrid => (true, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
)

"""
    _plot_power_spectrum!(fig, ax, df, channels_to_plot, fs; kwargs...)

Internal function to plot power spectra on an existing axis with interactive controls.

# Arguments
- `fig`: Makie Figure
- `ax`: Makie Axis to plot on
- `df`: DataFrame containing the data
- `channels_to_plot`: Vector of channel symbols to plot
- `fs`: Sampling frequency in Hz

$(generate_kwargs_doc(PLOT_POWER_SPECTRUM_KWARGS))

# Features
- Creates interactive checkboxes for toggling between linear and log scales
- Automatically sets up grid layout with controls in the top-right corner
- Shows frequency band indicators below the main plot when enabled

This is an internal function used by the public plotting functions.

# Example
```julia
# Basic usage with custom styling
fig, ax = plot_channel_spectrum(dat;
    title = "Custom Power Spectrum",
    line_width = 3,
    max_freq = 100.0,
    show_freq_bands = false)

# With custom colors and scales
fig, ax = plot_channel_spectrum(dat;
    x_scale = :log10,
    y_scale = :log10,
    line_alpha = 0.6,
    grid_visible = false)
```
"""
function _plot_power_spectrum!(fig, ax, df::DataFrame, channels_to_plot::Vector{Symbol}, fs::Real; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs)

    # Extract commonly used values
    window_size = plot_kwargs[:window_size]
    overlap = plot_kwargs[:overlap]
    max_freq = plot_kwargs[:max_freq]
    x_scale = plot_kwargs[:x_scale]
    y_scale = plot_kwargs[:y_scale]
    unit = plot_kwargs[:unit]
    window_function = plot_kwargs[:window_function]
    show_freq_bands = plot_kwargs[:show_freq_bands]

    # Validate scale parameters
    valid_scales = [:linear, :log10]
    if !(x_scale in valid_scales)
        @minimal_warning "Invalid x_scale '$x_scale'. Using :linear instead."
        x_scale = :linear
    end
    if !(y_scale in valid_scales)
        @minimal_warning "Invalid y_scale '$y_scale'. Using :linear instead."
        y_scale = :linear
    end

    # Validate unit parameter
    valid_units = [:linear, :dB]
    if !(unit in valid_units)
        @minimal_warning "Invalid unit '$unit'. Using :linear instead."
        unit = :linear
    end

    # Define EEG frequency bands
    freq_bands = Dict(
        "δ" => (0.5, 4.0),    # Delta
        "θ" => (4.0, 8.0),    # Theta
        "α" => (8.0, 13.0),   # Alpha
        "β" => (13.0, 30.0),  # Beta
        "γ" => (30.0, 100.0),  # Gamma
    )

    # Set axis labels and title
    ax.xlabel = plot_kwargs[:xlabel]
    # Update y-axis label based on unit
    ylabel_text = unit == :dB ? "Power Spectral Density (dB)" : plot_kwargs[:ylabel]
    ax.ylabel = ylabel_text
    ax.title = plot_kwargs[:title]
    ax.titlesize = plot_kwargs[:title_fontsize]
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    # Configure grid using the new axis styling function
    _set_axis_grid!(
        ax;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    # Create interactive controls in the figure
    controls_area = fig[1, 2] = GridLayout()

    # Set the column/row proportions 
    colsize!(fig.layout, 1, Relative(0.9))  # Main plot column
    colsize!(fig.layout, 2, Relative(0.1))  # Controls column
    rowsize!(fig.layout, 1, Relative(0.7))  # Main content row

    x_scale_obs = Observable(x_scale)
    y_scale_obs = Observable(y_scale)
    unit_obs = Observable(unit)

    # Add x/y checkboxes for axis types with labels
    x_log_checkbox = Checkbox(controls_area[1, 1], checked = x_scale == :log10)
    Label(controls_area[1, 2], "X: Linear/Log")
    y_log_checkbox = Checkbox(controls_area[2, 1], checked = y_scale == :log10)
    Label(controls_area[2, 2], "Y: Linear/Log")
    unit_checkbox = Checkbox(controls_area[3, 1], checked = unit == :dB)
    Label(controls_area[3, 2], "Unit: μV²/Hz/dB")

    # Update Observables when checkboxes change
    on(x_log_checkbox.checked) do checked
        x_scale_obs[] = checked ? :log10 : :linear
    end
    on(y_log_checkbox.checked) do checked
        y_scale_obs[] = checked ? :log10 : :linear
    end
    on(unit_checkbox.checked) do checked
        unit_obs[] = checked ? :dB : :linear
    end

    # Apply initial scale settings
    if x_scale == :log10
        xlims!(ax, (0.1, max_freq))
        ax.xscale = log10
    else
        xlims!(ax, (0, max_freq))
    end

    # Calculate and plot spectra for all channels in a single loop
    # Adjust window_size if signal is shorter than window_size
    n_samples = nrow(df)
    effective_window_size = min(window_size, n_samples)
    # Ensure window_size is at least 4 samples for meaningful FFT
    effective_window_size = max(4, effective_window_size)

    noverlap = Int(round(effective_window_size * overlap))

    # Store frequency and raw power data for reactive updates
    freq_data = []
    power_raw_data = []
    line_plots = []

    @views for ch in channels_to_plot
        signal = df[!, ch]
        # Use effective_window_size instead of window_size
        pgram = DSP.welch_pgram(signal, effective_window_size, noverlap; fs = fs, window = window_function)
        freqs, psd_raw = DSP.freq(pgram), DSP.power(pgram)

        # Convert power to requested unit
        if unit == :dB
            # Convert to dB: 10 * log10(power / reference)
            # Reference is 1 μV²/Hz (standard in EEG)
            psd = 10.0 .* log10.(max.(psd_raw, 1e-10))  # Avoid log(0) with small threshold
        else
            psd = psd_raw
        end

        # Store raw data for reactive updates
        push!(freq_data, freqs)
        push!(power_raw_data, psd_raw)

        # Create observable for power data
        psd_obs = Observable(psd)

        # Plot this channel's spectrum with observable
        line_plot = lines!(ax, freqs, psd_obs, label = string(ch), linewidth = plot_kwargs[:linewidth], alpha = plot_kwargs[:line_alpha])
        push!(line_plots, (line_plot, psd_obs))
    end

    # Function to get current power data (converted based on unit)
    function get_current_power_data(current_unit)
        if current_unit == :dB
            return [10.0 .* log10.(max.(psd_raw, 1e-10)) for psd_raw in power_raw_data]
        else
            return power_raw_data
        end
    end

    # Apply initial y-axis scale settings (after data is calculated)
    # This will be handled by update_y_scale, but we need to set initial state
    initial_power_data = get_current_power_data(unit)
    all_initial_powers = vcat(initial_power_data...)
    initial_max = maximum(all_initial_powers)
    initial_min = minimum(all_initial_powers)

    if y_scale == :log10
        # Calculate valid limits first
        if unit == :dB
            # For dB, use linear scale (dB is already logarithmic)
            ax.yscale = identity
            ylims!(ax, (initial_min - 5, initial_max + 5))
        else
            # For linear units, use log scale
            positive_powers = filter(x -> x > 0, all_initial_powers)
            if isempty(positive_powers)
                log_min, log_max = 0.001, 1.0
            else
                log_min = max(0.001, minimum(positive_powers))
                log_max = max(log_min * 10, initial_max)
            end
            # Set valid limits BEFORE changing scale to avoid validation errors
            ylims!(ax, (log_min, log_max))
            ax.yscale = log10
        end
    else
        # For linear scale, ensure axis is in linear mode first
        ax.yscale = identity
        if unit == :dB
            ylims!(ax, (initial_min - 5, initial_max + 5))
        else
            ylims!(ax, (0.0, max(0.1, initial_max)))
        end
    end

    # Add frequency band indicators below the x-axis if requested
    if show_freq_bands

        band_ax = Axis(fig[2, 1], height = 30)
        linkxaxes!(ax, band_ax)

        # Set limits and hide decorations
        xlims!(band_ax, 0, max_freq)
        ylims!(band_ax, 0, 1)
        hidedecorations!(band_ax)
        hidespines!(band_ax)

        band_colors = [:lightblue, :lightgreen, :yellow, :orange, :lightcoral]
        for (i, (band_name, (fmin, fmax))) in enumerate(freq_bands)
            if fmax <= max_freq
                # Add colored bar and label
                bar_x = [fmin, fmax]
                bar_y = [0.5, 0.5]
                lines!(band_ax, bar_x, bar_y, color = band_colors[i], linewidth = 8, alpha = plot_kwargs[:freq_band_alpha])
                text!(band_ax, (fmin + fmax) / 2, 0.5, text = band_name, align = (:center, :center), fontsize = 26, color = :black)
            end
        end
    end

    # Add reactive updates for axis scales
    # Function to update x-axis scale
    function update_x_scale(scale)
        if scale == :log10
            xlims!(ax, (0.1, max_freq))
            ax.xscale = log10
        else
            xlims!(ax, (0.1, max_freq))
            ax.xscale = identity
        end
    end

    # Function to update unit
    function update_unit(new_unit)
        # Convert power data
        converted_power_data = get_current_power_data(new_unit)

        # Update line plots
        for (i, (_, psd_obs)) in enumerate(line_plots)
            psd_obs[] = converted_power_data[i]
        end

        # Update y-axis label
        ax.ylabel = new_unit == :dB ? "Power Spectral Density (dB)" : plot_kwargs[:ylabel]

        # Recalculate max power and update limits
        all_powers = vcat(converted_power_data...)
        new_max_power = maximum(all_powers)
        new_min_power = minimum(all_powers)

        # Update y-axis scale and limits based on current y_scale setting
        if y_scale_obs[] == :log10
            if new_unit == :dB
                # For dB, use linear scale (dB is already logarithmic)
                ax.yscale = identity
                ylims!(ax, (new_min_power - 5, new_max_power + 5))
            else
                # For linear units, use log scale
                positive_powers = filter(x -> x > 0, all_powers)
                if isempty(positive_powers)
                    log_min, log_max = 0.001, 1.0
                else
                    log_min = max(0.001, minimum(positive_powers))
                    log_max = max(log_min * 10, new_max_power)
                end
                ylims!(ax, (log_min, log_max))
                ax.yscale = log10
            end
        else
            ax.yscale = identity
            if new_unit == :dB
                ylims!(ax, (new_min_power - 5, new_max_power + 5))
            else
                ylims!(ax, (0.0, max(0.1, new_max_power)))
            end
        end
    end

    # Function to update y-axis scale
    function update_y_scale(scale)
        current_unit = unit_obs[]
        current_power_data = get_current_power_data(current_unit)
        all_powers = vcat(current_power_data...)
        current_max = maximum(all_powers)
        current_min = minimum(all_powers)

        if scale == :log10
            if current_unit == :dB
                # For dB, use linear scale (dB is already logarithmic)
                ax.yscale = identity
                ylims!(ax, (current_min - 5, current_max + 5))
            else
                # Calculate valid limits first
                positive_powers = filter(x -> x > 0, all_powers)
                if isempty(positive_powers)
                    # Fallback: no positive values, use safe defaults
                    log_min = 0.001
                    log_max = 1.0
                else
                    log_min = max(0.001, minimum(positive_powers))
                    log_max = max(log_min * 10, current_max)  # Ensure reasonable range
                end
                # Set valid limits BEFORE changing scale to avoid validation errors
                ylims!(ax, (log_min, log_max))
                ax.yscale = log10
            end
        else
            # For linear scale, change scale FIRST to avoid log validation of linear limits
            ax.yscale = identity
            if current_unit == :dB
                ylims!(ax, (current_min - 5, current_max + 5))
            else
                ylims!(ax, (0.0, max(0.1, current_max)))
            end
        end
    end

    # Set up reactive updates
    on(x_scale_obs) do scale
        update_x_scale(scale)
    end

    on(y_scale_obs) do scale
        update_y_scale(scale)
    end

    on(unit_obs) do new_unit
        update_unit(new_unit)
    end


end


"""
    plot_channel_spectrum(dat; kwargs...)

Plot power spectrum for specified channels from EEG data with interactive controls.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024)
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5)
- `max_freq::Real`: Maximum frequency to display in Hz (default: 200.0)
- `x_scale::Symbol`: Initial scale for x-axis, one of :linear or :log10 (default: :linear)
- `y_scale::Symbol`: Initial scale for y-axis, one of :linear or :log10 (default: :linear)
- `unit::Symbol`: Power unit, one of :linear (μV²/Hz) or :dB (decibels) (default: :linear)
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning)
- `show_freq_bands::Bool`: Whether to show frequency band indicators (default: true)
- `show_legend::Bool`: Whether to show the legend (default: true)
- `display_plot::Bool`: Whether to display the plot (default: true)

# Features
- Interactive checkboxes to toggle between linear and log scales for both axes
- Automatic multi-column legend for many channels
- Frequency band indicators (δ, θ, α, β, γ) below the main plot
- Responsive grid layout with controls in the top-right corner

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot
- `ax::Axis`: The axis containing the plot

# Examples
```julia
# Basic usage with interactive controls
fig, ax = plot_channel_spectrum(dat)

# With specific channels
fig, ax = plot_channel_spectrum(dat, 
    channel_selection = channels([:Fp1, :Fp2]))

# With log scales
fig, ax = plot_channel_spectrum(dat, 
    x_scale = :log10, 
    y_scale = :log10)

# With dB units (EEGLAB-style)
fig, ax = plot_channel_spectrum(dat, unit = :dB)
```
"""

function plot_channel_spectrum(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    interval_selection::TimeInterval = times(),
    channel_selection::Function = channels(),
    kwargs...,
)
    # Generate window title from dataset
    title_str = _generate_window_title(dat)
    set_window_title(title_str)

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs)

    # data selection
    dat_subset = subset(dat, sample_selection = sample_selection, interval_selection = interval_selection, channel_selection = channel_selection)

    # Create figure and main axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    _plot_power_spectrum!(fig, ax, dat_subset.data, dat_subset.layout.data.label, sample_rate(dat_subset); kwargs...)

    # Add legend if requested
    if plot_kwargs[:show_legend]
        n_channels = length(dat_subset.layout.data.label)
        n_cols = n_channels > 10 ? cld(n_channels, 20) : 1
        axislegend(ax, nbanks = n_cols)
    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    set_window_title("Makie")
    return fig, ax

end


"""
    plot_ica_component_spectrum(dat::ContinuousData, ica_result::InfoIca; component_selection::Function = components(), kwargs...)

Plot power spectrum of ICA component(s) with interactive controls for axis scaling.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `ica_result::InfoIca`: The ICA result object.
- `component_selection::Function`: Function that returns boolean vector for component filtering (default: all components).

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Initial scale for x-axis, one of :linear or :log10 (default: :linear).
- `y_scale::Symbol`: Initial scale for y-axis, one of :linear or :log10 (default: :linear).
- `unit::Symbol`: Power unit, one of :linear (μV²/Hz) or :dB (decibels) (default: :linear).
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning).
- `show_freq_bands::Bool`: Whether to show frequency band indicators (default: true).
- `show_legend::Bool`: Whether to show the legend with component variance percentages (default: true).
- `display_plot::Bool`: Whether to display the plot (default: true).

# Features
- Interactive checkboxes to toggle between linear and log scales for both axes
- Automatic multi-column legend for many components with variance percentages
- Frequency band indicators (δ, θ, α, β, γ) below the main plot
- Responsive grid layout with controls in the top-right corner

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.

# Examples
```julia
# Plot all components (default)
plot_ica_component_spectrum(dat, ica_result)

# Plot specific components
plot_ica_component_spectrum(dat, ica_result, component_selection = components([1, 3, 5]))

# Plot components 1-10
plot_ica_component_spectrum(dat, ica_result, component_selection = components(1:10))

# Plot all components except 1 and 2
plot_ica_component_spectrum(dat, ica_result, component_selection = components_not([1, 2]))

# Plot component 1 only
plot_ica_component_spectrum(dat, ica_result, component_selection = components(1))

# With sample selection
plot_ica_component_spectrum(dat, ica_result, component_selection = components(1), sample_selection = samples_not(:is_extreme_value_100))

# With dB units (EEGLAB-style)
plot_ica_component_spectrum(dat, ica_result, component_selection = components(1), unit = :dB)
```
"""

function plot_ica_component_spectrum(
    dat::ContinuousData,
    ica_result::InfoIca;
    sample_selection::Function = samples(),
    interval_selection::TimeInterval = times(),
    component_selection::Function = components(),
    kwargs...,
)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs)
    # Get selected components using the predicate
    selected_components = get_selected_components(ica_result, component_selection)

    # If empty vector, use all components
    if isempty(selected_components)
        selected_components = 1:size(ica_result.unmixing, 1)
    end

    # Debug: verify selected components
    @debug "plot_ica_component_spectrum: selected_components = $selected_components"
    # Apply sample selection
    selected_samples = get_selected_samples(dat, sample_selection)

    # Prepare data matrix for selected samples only
    relevant_cols = vcat(ica_result.layout.data.label)
    dat_matrix = permutedims(Matrix(dat.data[selected_samples, relevant_cols]))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica_result.scale

    # Calculate component activations for selected samples
    components = ica_result.unmixing * dat_matrix
    fs = dat.sample_rate

    # Create a DataFrame to store all component signals
    component_df = DataFrame(:time => dat.data.time[selected_samples])

    # Add each component to the DataFrame with consistent naming
    for comp_idx in selected_components
        component_name = Symbol("IC_$(comp_idx)")
        component_df[!, component_name] = components[comp_idx, :]
    end

    # Create list of component symbols to plot
    components_to_plot = [Symbol("IC_$(comp_idx)") for comp_idx in selected_components]

    # Create figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])

    _plot_power_spectrum!(fig, ax, component_df, components_to_plot, sample_rate(dat); kwargs...)

    # Add custom legend with component variance percentages if requested
    if plot_kwargs[:show_legend]
        # Create component labels with variance percentages
        component_labels = []
        for comp_idx in selected_components
            variance_pct = ica_result.variance[comp_idx] * 100
            push!(component_labels, @sprintf("IC %d (%.1f%%)", comp_idx, variance_pct))
        end

        # Use axislegend for consistency with plot_channel_spectrum
        # Calculate number of columns based on number of components
        ncols = length(selected_components) > 10 ? ceil(Int, length(selected_components) / 10) : 1

        # Place legend on the axis itself
        axislegend(ax, nbanks = ncols)
    end

    if plot_kwargs[:display_plot]
        display_figure(fig)
    end

    return fig, ax
end
