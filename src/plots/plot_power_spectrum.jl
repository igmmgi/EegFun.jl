
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

    # Font sizes
    :title_fontsize => (16, "Font size for title"),
    :label_fontsize => (14, "Font size for axis labels"),
    :tick_fontsize => (12, "Font size for tick labels"),
    :legend_fontsize => (12, "Font size for legend"),

    # Line styling
    :line_width => (2, "Line width for spectrum lines"),
    :line_alpha => (0.8, "Transparency for spectrum lines"),

    # Frequency band styling
    :freq_band_alpha => (0.3, "Transparency for frequency band indicators"),
    :freq_band_height => (0.1, "Height of frequency band indicators"),

    # Grid styling
    :grid_visible => (true, "Whether to show grid"),
    :grid_alpha => (0.3, "Transparency of grid"),
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
    ax.ylabel = plot_kwargs[:ylabel]
    ax.title = plot_kwargs[:title]
    ax.titlesize = plot_kwargs[:title_fontsize]
    ax.xlabelsize = plot_kwargs[:label_fontsize]
    ax.ylabelsize = plot_kwargs[:label_fontsize]
    ax.xticklabelsize = plot_kwargs[:tick_fontsize]
    ax.yticklabelsize = plot_kwargs[:tick_fontsize]

    # Configure grid
    ax.xgridvisible = plot_kwargs[:grid_visible]
    ax.ygridvisible = plot_kwargs[:grid_visible]
    ax.xgridwidth = 1
    ax.ygridwidth = 1
    ax.xgridcolor = (:gray, plot_kwargs[:grid_alpha])
    ax.ygridcolor = (:gray, plot_kwargs[:grid_alpha])

    # Create interactive controls in the figure
    controls_area = fig[1, 2] = GridLayout()

    # Set the column/row proportions 
    colsize!(fig.layout, 1, Relative(0.9))  # Main plot column
    colsize!(fig.layout, 2, Relative(0.1))  # Controls column
    rowsize!(fig.layout, 1, Relative(0.7))  # Main content row

    x_scale_obs = Observable(x_scale)
    y_scale_obs = Observable(y_scale)

    # Add x/y checkboxes for axis types with labels
    x_log_checkbox = Checkbox(controls_area[1, 1], checked = x_scale == :log10)
    Label(controls_area[1, 2], "X: Linear/Log")
    y_log_checkbox = Checkbox(controls_area[2, 1], checked = y_scale == :log10)
    Label(controls_area[2, 2], "Y: Linear/Log")

    # Update Observables when checkboxes change
    on(x_log_checkbox.checked) do checked
        x_scale_obs[] = checked ? :log10 : :linear
    end
    on(y_log_checkbox.checked) do checked
        y_scale_obs[] = checked ? :log10 : :linear
    end

    # Apply initial scale settings
    if x_scale == :log10
        xlims!(ax, (0.1, max_freq))
        ax.xscale = log10
    else
        xlims!(ax, (0, max_freq))
    end

    # Calculate and plot spectra for all channels in a single loop
    noverlap = Int(round(window_size * overlap))
    max_power = 0.0

    # Store frequency and power data for reactive updates
    freq_data = []
    power_data = []

    @views for ch in channels_to_plot
        signal = df[!, ch]
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = fs, window = window_function)
        freqs, psd = DSP.freq(pgram), DSP.power(pgram)

        # Track max power for y-axis limits
        max_power = max(max_power, maximum(psd))

        # Store data for reactive updates
        push!(freq_data, freqs)
        push!(power_data, psd)

        # Plot this channel's spectrum
        lines!(
            ax,
            freqs,
            psd,
            label = string(ch),
            linewidth = plot_kwargs[:line_width],
            alpha = plot_kwargs[:line_alpha],
        )
    end

    # Apply initial y-axis scale settings (after data is calculated)
    if y_scale == :log10
        # Calculate valid limits first
        positive_powers = Base.filter(x -> x > 0, vcat(power_data...))
        if isempty(positive_powers)
            log_min, log_max = 0.001, 1.0
        else
            log_min = max(0.001, minimum(positive_powers))
            log_max = max(log_min * 10, max_power)
        end
        # Set valid limits BEFORE changing scale to avoid validation errors
        ylims!(ax, (log_min, log_max))
        ax.yscale = log10
    else
        # For linear scale, ensure axis is in linear mode first
        ax.yscale = identity
        ylims!(ax, (0.0, max(0.1, max_power)))
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
                lines!(
                    band_ax,
                    bar_x,
                    bar_y,
                    color = band_colors[i],
                    linewidth = 8,
                    alpha = plot_kwargs[:freq_band_alpha],
                )
                text!(
                    band_ax,
                    (fmin + fmax) / 2,
                    0.5,
                    text = band_name,
                    align = (:center, :center),
                    fontsize = 26,
                    color = :black,
                )
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

    # Function to update y-axis scale
    function update_y_scale(scale)
        if scale == :log10
            # Calculate valid limits first
            positive_powers = Base.filter(x -> x > 0, vcat(power_data...))
            if isempty(positive_powers)
                # Fallback: no positive values, use safe defaults
                log_min = 0.001
                log_max = 1.0
            else
                log_min = max(0.001, minimum(positive_powers))
                log_max = max(log_min * 10, max_power)  # Ensure reasonable range
            end
            # Set valid limits BEFORE changing scale to avoid validation errors
            ylims!(ax, (log_min, log_max))
            ax.yscale = log10
        else
            # For linear scale, change scale FIRST to avoid log validation of linear limits
            ax.yscale = identity
            ylims!(ax, (0.0, max(0.1, max_power)))
        end
    end

    # Set up reactive updates
    on(x_scale_obs) do scale
        update_x_scale(scale)
    end

    on(y_scale_obs) do scale
        update_y_scale(scale)
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
```
"""

function plot_channel_spectrum(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    kwargs...,
)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_POWER_SPECTRUM_KWARGS, kwargs)

    # data selection
    dat_subset = subset(dat, sample_selection = sample_selection, channel_selection = channel_selection)

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

    return fig, ax

end


"""
    plot_component_spectrum(ica_result::InfoIca, dat::ContinuousData; component_selection::Function = components(), kwargs...)

Plot power spectrum of ICA component(s) with interactive controls for axis scaling.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.
- `component_selection::Function`: Function that returns boolean vector for component filtering (default: all components).

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Initial scale for x-axis, one of :linear or :log10 (default: :linear).
- `y_scale::Symbol`: Initial scale for y-axis, one of :linear or :log10 (default: :linear).
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
plot_component_spectrum(ica_result, dat)

# Plot specific components
plot_component_spectrum(ica_result, dat, component_selection = components([1, 3, 5]))

# Plot components 1-10
plot_component_spectrum(ica_result, dat, component_selection = components(1:10))

# Plot all components except 1 and 2
plot_component_spectrum(ica_result, dat, component_selection = components_not([1, 2]))

# Plot component 1 only
plot_component_spectrum(ica_result, dat, component_selection = components(1))

# With sample selection
plot_component_spectrum(ica_result, dat, component_selection = components(1), sample_selection = samples_not(:is_extreme_value_100))
```
"""

function plot_component_spectrum(
    ica_result::InfoIca,
    dat::ContinuousData;
    sample_selection::Function = samples(),
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
