"""
    _plot_power_spectrum!(fig, ax, df, channels_to_plot, fs; kwargs...)

Internal function to plot power spectra on an existing axis.

# Arguments
- `fig`: Makie Figure
- `ax`: Makie Axis to plot on
- `df`: DataFrame containing the data
- `channels_to_plot`: Vector of channel symbols to plot
- `fs`: Sampling frequency in Hz

# Keyword Arguments
- `window_size::Int`: Size of the FFT window (default: 1024)
- `overlap::Real`: Overlap between windows (default: 0.5)
- `max_freq::Real`: Maximum frequency to display (default: 200.0)
- `x_scale::Symbol`: X-axis scale, :linear or :log10 (default: :linear)
- `y_scale::Symbol`: Y-axis scale, :linear or :log10 (default: :linear)
- `window_function::Function`: Window function for spectral estimation (default: DSP.hanning)
- `show_freq_bands::Bool`: Whether to show frequency band indicators (default: true)

This is an internal function used by the public plotting functions.
"""
function _plot_power_spectrum!(
    fig,
    ax,
    df::DataFrame,
    channels_to_plot::Vector{Symbol},
    fs::Real,
    controls_area;
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 200.0,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    window_function::Function = DSP.hanning,
    show_freq_bands::Bool = true,
)

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
        "γ" => (30.0, 100.0)  # Gamma
    )

    # Set axis labels and title
    ax.xlabel = "Frequency (Hz)"
    ax.ylabel = "Power Spectral Density (μV²/Hz)"
    ax.title = "Power Spectrum"

    # Create interactive controls
    x_scale_obs = Observable(x_scale)
    y_scale_obs = Observable(y_scale)
    
    # Add checkboxes for axis types with labels
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

      # Apply scale settings
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
        lines!(ax, freqs, psd, label = string(ch))
    end
    
    # Add frequency band indicators below the x-axis if requested
    if show_freq_bands
        # Create a small axis below the main plot for frequency bands
        band_ax = Axis(fig[2, 1], height = 30)
        linkxaxes!(ax, band_ax)
        
        # Set limits and hide decorations
        xlims!(band_ax, 0, max_freq)
        ylims!(band_ax, 0, 1)
        hideydecorations!(band_ax)
        
        # Hide all spines and ticks
        band_ax.leftspinevisible = false
        band_ax.rightspinevisible = false
        band_ax.topspinevisible = false
        band_ax.bottomspinevisible = false
        band_ax.xticksvisible = false
        band_ax.xticklabelsvisible = false
        
        band_colors = [:lightblue, :lightgreen, :yellow, :orange, :lightcoral]
        for (i, (band_name, (fmin, fmax))) in enumerate(freq_bands)
            if fmax <= max_freq
                # Add colored bar
                bar_x = [fmin, fmax]
                bar_y = [0.5, 0.5]
                lines!(band_ax, bar_x, bar_y, color = band_colors[i], linewidth = 8)
                
                # Add band label
                text!(band_ax, (fmin + fmax) / 2, 0.5, text = band_name, 
                      align = (:center, :center), fontsize = 26, color = :black)
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
            ylims!(ax, (0.1, max_power))
            ax.yscale = log10
        else
            ylims!(ax, (0.1, max_power))
            ax.yscale = identity
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
    plot_channel_spectrum(data; kwargs...)
    plot_channel_spectrum(data, fs; kwargs...)

Plot power spectrum for specified channels from EEG data.

# Arguments
- `data`: DataFrame or ContinuousData object containing EEG data
- `fs`: Sampling frequency in Hz (optional, extracted from data if not provided)

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0)
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz)
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024)
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5)
- `max_freq::Real`: Maximum frequency to display in Hz (default: 200.0)
- `x_scale::Symbol`: Scale for x-axis, one of :linear or :log10 (default: :linear)
- `y_scale::Symbol`: Scale for y-axis, one of :linear or :log10 (default: :linear)
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning)
- `show_legend::Bool`: Whether to show the legend (default: true)

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
    show_legend::Bool = true,
    display_plot::Bool = true,
    kwargs...,
)

    # data selection
    dat_subset = subset(dat, sample_selection = sample_selection, channel_selection = channel_selection)

    fig = Figure()
    
    # Create a proper grid layout with interactive controls
    # Main plot takes most of the width and full height
    main_plot_area = fig[1, 1] = GridLayout()
    # Controls take a small portion of width and are positioned at the top
    controls_area = fig[1, 2] = GridLayout()
    
    # Set the column proportions using colsize!
    colsize!(fig.layout, 1, Relative(0.9))  # Main plot column
    colsize!(fig.layout, 2, Relative(0.1))  # Controls column
    
    # Set the row proportions - controls only take top portion
    rowsize!(fig.layout, 1, Relative(0.7))  # Main content row
    
    ax = Axis(main_plot_area[1, 1])

    _plot_power_spectrum!(fig, ax, dat_subset.data, dat_subset.layout.data.label, sample_rate(dat_subset), controls_area; kwargs...)
    
    # Add legend if requested
    if show_legend
        # Automatically arrange legend in multiple columns if many channels
        n_channels = length(dat_subset.layout.data.label)
        n_cols = n_channels > 10 ? cld(n_channels, 20) : 1
        axislegend(ax, nbanks = n_cols)
    end
    if display_plot
        display_figure(fig)
    end

    return fig, ax

end


"""
    plot_component_spectrum(ica_result::InfoIca, dat::ContinuousData; component_selection::Function = components(), kwargs...)

Plot power spectrum of ICA component(s) using a component selection predicate.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.
- `component_selection::Function`: Function that returns boolean vector for component filtering (default: all components).

# Keyword Arguments
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: samples() - all samples)
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Scale for x-axis, one of :linear or :log10 (default: :linear).
- `y_scale::Symbol`: Scale for y-axis, one of :linear or :log10 (default: :linear).
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning).
- `show_legend::Bool`: Whether to show the legend with component variance percentages (default: true).
- `display_plot::Bool`: Whether to display the plot (default: true).

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
    show_legend::Bool = true,
    display_plot::Bool = true,
    kwards...,
)
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
    if show_legend
        # Create component labels with variance percentages
        component_labels = []
        for comp_idx in selected_components
            variance_pct = ica_result.variance[comp_idx] * 100
            push!(component_labels, @sprintf("IC %d (%.1f%%)", comp_idx, variance_pct))
        end
        
        # Create legend entries from the plotted lines
        legend_entries = [LineElement(color = ax.scene.plots[i].color) for i = 1:length(selected_components)]
        
        # Calculate number of columns based on number of components
        ncols = length(selected_components) > 10 ? ceil(Int, length(selected_components) / 10) : 1
        
        # Place legend - use multi-column for many components
        Legend(fig[1, 2], legend_entries, component_labels, "Components", nbanks = ncols)
    end

    if display_plot
        display_figure(fig)
    end

    return fig, ax
end
