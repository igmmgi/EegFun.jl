function _plot_power_spectrum_implementation(
    df::DataFrame,
    channels_to_plot::Vector{Symbol},
    fs::Real;
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 200.0,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    window_function::Function = DSP.hanning,
    show_legend::Bool = true,
    display_plot::Bool = true,
)
    # Check if we have enough data and adjust window size if needed
    if size(df, 1) < window_size
        @warn "Selected region is too short for the specified window size. Adjusting window size..."
        window_size = 2^floor(Int, log2(size(df, 1) / 2))
        if window_size < 32
            @warn "Selected region is too short for meaningful spectral analysis. Minimum window size of 32 samples required."
            return Figure(), Axis()
        end
    end

    # Validate scale parameters
    valid_scales = [:linear, :log10]
    if !(x_scale in valid_scales)
        @warn "Invalid x_scale '$x_scale'. Using :linear instead."
        x_scale = :linear
    end
    if !(y_scale in valid_scales)
        @warn "Invalid y_scale '$y_scale'. Using :linear instead."
        y_scale = :linear
    end

    # Prepare signals dictionary
    signals = Dict{Symbol,Vector{Float64}}()
    for ch in channels_to_plot
        signals[ch] = Vector{Float64}(df[!, ch])
    end

    # Calculate spectra for all requested channels
    noverlap = Int(round(window_size * overlap))
    spectra = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()

    for ch in channels_to_plot
        signal = signals[ch]
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = fs, window = window_function)
        spectra[ch] = (DSP.freq(pgram), DSP.power(pgram))
    end

    # Create figure with a layout for plot and controls
    fig = Figure()
    
    # Create main grid layout with plot area and control area
    gl = fig[1, 1] = GridLayout()
    
    # Create a toggle grid for the controls
    control_grid = fig[2, 1] = GridLayout()
    
    # Initialize with the specified scales
    is_xlog10 = x_scale == :log10
    is_ylog10 = y_scale == :log10
    
    # Create axis with default linear scales first
    ax = Axis(gl[1, 1]; 
        xlabel = "Frequency (Hz)",
        ylabel = "Power Spectral Density (μV²/Hz)",
        title = "Channel Power Spectrum"
    )
    
    # Apply initial scale settings and limits separately, in the correct order
    if is_xlog10
        xlims!(ax, (0.1, max_freq))
        ax.xscale = log10
    else
        xlims!(ax, (0, max_freq))
    end
    
    # Add the same initialization for y-scale
    # Calculate max y-value for limits
    max_power = maximum([maximum(psd) for (_, psd) in values(spectra)])
    if is_ylog10
        ylims!(ax, (1e-10, max_power))
        ax.yscale = log10
    else
        ylims!(ax, (0, max_power))
    end

    # Add scale menus for each axis
    Label(control_grid[1, 1], "X-axis scale:")
    x_scale_options = ["linear", "log10"]
    x_scale_menu = Menu(control_grid[1, 2], 
                        options = x_scale_options,
                        default = x_scale == :log10 ? "log10" : "linear")
    
    Label(control_grid[2, 1], "Y-axis scale:")
    y_scale_options = ["linear", "log10"]
    y_scale_menu = Menu(control_grid[2, 2], 
                        options = y_scale_options,
                        default = y_scale == :log10 ? "log10" : "linear")
    
    # Connect menu selection to axis scales
    on(x_scale_menu.selection) do selection
        if selection == "linear"
            # Apply linear scale
            ax.xscale = identity
            xlims!(ax, (0, max_freq))
        else # "log10"
            # Apply log10 scale
            xlims!(ax, (0.1, max_freq))
            ax.xscale = log10
        end
    end
    
    on(y_scale_menu.selection) do selection
        if selection == "linear"
            # Apply linear scale
            ax.yscale = identity
            ylims!(ax, (0, max_power))
        else # "log10"
            # Apply log10 scale
            ylims!(ax, (1e-10, max_power))
            ax.yscale = log10
        end
    end

    # Plot spectra for all channels
    for ch in channels_to_plot
        freqs, psd = spectra[ch]
        lines!(ax, freqs, psd, label = string(ch))
    end

    # Calculate max value for vertical lines
    max_power = maximum([maximum(psd) for (_, psd) in values(spectra)])

    # Highlight line frequency and harmonics
    for h = 1:3
        freq = line_freq * h
        if freq <= max_freq
            # Add vertical line
            vlines!(
                ax,
                [freq],
                color = :red,
                linestyle = :dash,
            )

            # Add shaded region
            band_x = [freq - freq_bandwidth, freq + freq_bandwidth]
            band_y = [0, max_power]
            poly!(
                ax,
                [
                    Point2f(band_x[1], band_y[1]),
                    Point2f(band_x[2], band_y[1]),
                    Point2f(band_x[2], band_y[2]),
                    Point2f(band_x[1], band_y[2]),
                ],
                color = (:red, 0.1),
            )

        end
    end

    if show_legend
        # Calculate number of columns for the legend based on the number of entries
        ncols = length(channels_to_plot) > 10 ? ceil(Int, length(channels_to_plot) / 10) : 1
        axislegend(ax, position = (1.0, 1.0), nbanks = ncols)
    end

    if display_plot
        display(fig)
    end

    return fig, ax
end



"""
    plot_channel_spectrum(df::DataFrame, fs::Real, channels::Function = channels();
                          kwargs = Dict())

Plot power spectrum for specified channels.

# Arguments
- `df::DataFrame`: DataFrame containing EEG data
- `fs::Real`: Sampling frequency in Hz
- `channels::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `kwargs`: Additional keyword arguments for customization

# Examples
```julia
# Plot all channels
plot_channel_spectrum(df, 1000)

# Plot specific channels
plot_channel_spectrum(df, 1000, channels([:Fp1, :Fp2]))

# Exclude reference channels
plot_channel_spectrum(df, 1000, channels_not([:M1, :M2]))

# Plot frontal channels only
plot_channel_spectrum(df, 1000, channels(1:10))

# Custom predicate
plot_channel_spectrum(df, 1000, x -> startswith.(string.(x), "F"))
```

# Keyword Arguments
- `xlim`: X-axis limits (default: auto-calculated)
- `ylim`: Y-axis limits (default: auto-calculated)
- `xlabel`: X-axis label (default: "Frequency (Hz)")
- `ylabel`: Y-axis label (default: "Power")
- `dims`: Grid dimensions [rows, cols] (default: auto-calculated)
- `hidedecorations`: Whether to hide axis decorations (default: false)
- `theme_fontsize`: Font size for theme (default: 24)
- `yreversed`: Whether to reverse Y-axis (default: false)
- `colormap`: Colormap for the plot (default: :jet)
- `colorrange`: Color range for the plot (default: auto-calculated)
"""
function plot_channel_spectrum(df::DataFrame, fs::Real, channels::Function = channels();
                              kwargs = Dict())
    # Get all column names that are not time, sample, or triggers
    all_columns = filter(col -> !(col in [:time, :sample, :triggers]), propertynames(df))
    
    # Filter channels
    channel_mask = channels(all_columns)
    selected_channels = all_columns[channel_mask]
    
    # Delegate to the implementation function
    return _plot_power_spectrum_implementation(
        df,
        selected_channels,
        fs;
        kwargs...
    )
end

# Backward compatibility - keep the old method for existing code
function plot_channel_spectrum(
    df::DataFrame,
    channel::Union{Symbol,Vector{Symbol}};
    kwargs...
)
    # Process channel input to get vector of channels to plot
    channels_to_plot = channel isa Symbol ? [channel] : channel
    # Check if channels exist
    for ch in channels_to_plot
        if !(ch in propertynames(df))
            error("Channel $ch not found in data")
        end
    end

    # Use the implementation function with filtered data
    return _plot_power_spectrum_implementation(
        df,
        channels_to_plot,
        sample_rate(df);
        kwargs...
    )
end

"""
    plot_channel_spectrum(dat::ContinuousData, channels::Function; kwargs...)

Plot the power spectrum of specified channels from a ContinuousData object.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `channels::Function`: Function that returns boolean vector for channel filtering.
- `kwargs`: Additional keyword arguments passed to plot_channel_spectrum.

# Examples
```julia
# Plot specific channels
plot_channel_spectrum(dat, channels([:Fp1, :Fp2]))

# Exclude reference channels
plot_channel_spectrum(dat, channels_not([:M1, :M2]))

# Plot frontal channels only
plot_channel_spectrum(dat, channels(1:10))
```
"""
function plot_channel_spectrum(dat::ContinuousData, channels::Function; kwargs...)
    return plot_channel_spectrum(dat.data, dat.sample_rate, channels; kwargs...)
end

"""
    plot_selected_spectrum(selected_data::DataFrame, channels::Function = channels();
                          line_freq::Real=50.0,
                          freq_bandwidth::Real=1.0,
                          window_size::Int=1024,
                          overlap::Real=0.5,
                          max_freq::Real=100.0,
                          x_scale::Symbol=:linear,
                          y_scale::Symbol=:linear)

Plot the power spectrum of a selected time region for specific channel(s).

# Arguments
- `selected_data::DataFrame`: The selected time region data.
- `channels::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels).

# Examples
```julia
# Plot all channels in selected region
plot_selected_spectrum(selected_data, channels)

# Plot specific channels
plot_selected_spectrum(selected_data, channels([:Fp1, :Fp2]))

# Exclude reference channels
plot_selected_spectrum(selected_data, channels_not([:M1, :M2]))

# Plot frontal channels only
plot_selected_spectrum(selected_data, channels(1:10))
```

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Scale for x-axis, one of :linear, :log10, or :log (default: :linear).
- `y_scale::Symbol`: Scale for y-axis, one of :linear, :log10, or :log (default: :linear).
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning).
  Common options include DSP.hamming, DSP.blackman, DSP.bartlett, DSP.rect (rectangular window).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.
"""
function plot_selected_spectrum(
    selected_data::DataFrame,
    channels::Function = channels();
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 200.0,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    window_function::Function = DSP.hanning,
)
    # Get sample rate from the data (assuming it's available)
    fs = sample_rate(selected_data)
    
    # Delegate to the main plotting function
    fig, ax = plot_channel_spectrum(
        selected_data,
        fs,
        channels;
        line_freq = line_freq,
        freq_bandwidth = freq_bandwidth,
        window_size = window_size,
        overlap = overlap,
        max_freq = max_freq,
        x_scale = x_scale,
        y_scale = y_scale,
        window_function = window_function,
    )

    return fig, ax
end

# Backward compatibility - keep the old method for existing code
function plot_selected_spectrum(
    selected_data::DataFrame,
    channel::Union{Symbol,Vector{Symbol}};
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 200.0,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    window_function::Function = DSP.hanning,
)
    # Convert to predicate for backward compatibility
    if channel isa Symbol
        channel_predicate = x -> x .== channel
    else
        channel_predicate = x -> x .∈ Ref(channel)
    end
    
    return plot_selected_spectrum(
        selected_data,
        channel_predicate;
        line_freq = line_freq,
        freq_bandwidth = freq_bandwidth,
        window_size = window_size,
        overlap = overlap,
        max_freq = max_freq,
        x_scale = x_scale,
        y_scale = y_scale,
        window_function = window_function,
    )
end






"""
    plot_component_spectrum(ica_result::InfoIca, dat::ContinuousData, comp_idx::Int;
                         line_freq::Real=50.0,
                         freq_bandwidth::Real=1.0,
                         window_size::Int=1024,
                         overlap::Real=0.5,
                         max_freq::Real=100.0,
                         x_scale::Symbol=:linear,
                         y_scale::Symbol=:linear,
                         window_function::Function = DSP.hanning)

Plot the power spectrum of a specific ICA component.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.
- `comp_idx::Int`: Index of the component to plot.

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Scale for x-axis, one of :linear or :log10 (default: :linear).
- `y_scale::Symbol`: Scale for y-axis, one of :linear or :log10 (default: :linear).
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.
"""
function plot_component_spectrum(
    ica_result::InfoIca,
    dat::ContinuousData,
    comp_idx::Int;
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 100.0,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    window_function::Function = DSP.hanning,
)
    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    dat_matrix = permutedims(Matrix(dat.data[:, relevant_cols]))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica_result.scale

    # Calculate component activation
    components = ica_result.unmixing * dat_matrix
    fs = dat.sample_rate

    # Get the component time series and create a temporary DataFrame
    signal = components[comp_idx, :]
    component_name = Symbol("Component_$(comp_idx)")
    time_values = dat.data.time
    component_df = DataFrame(component_name => signal, :time => dat.data.time)

    # Custom title for the component
    variance_pct = ica_result.variance[comp_idx] * 100
    title_text = @sprintf("Power Spectrum of Component %d (%.1f%% variance)", comp_idx, variance_pct)
    
    # Call the implementation function with the component data
    fig, ax = _plot_power_spectrum_implementation(
        component_df,
        [component_name],
        fs,
        line_freq,
        freq_bandwidth,
        window_size,
        overlap,
        max_freq,
        x_scale,
        y_scale,
        window_function,
    )
    
    # Update the plot title with component-specific information
    ax.title = title_text

    return fig, ax
end

"""
    plot_components_spectra(ica_result::InfoIca, dat::ContinuousData, comp_indices::Vector{Int}=Int[];
                          line_freq::Real=50.0,
                          freq_bandwidth::Real=1.0,
                          window_size::Int=1024,
                          overlap::Real=0.5,
                          max_freq::Real=100.0,
                          x_scale::Symbol=:linear,
                          y_scale::Symbol=:linear,
                          window_function::Function=DSP.hanning)

Plot power spectra of multiple ICA components on the same axis.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.
- `comp_indices::Vector{Int}`: Vector of component indices to plot. If empty, plots all components.

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `x_scale::Symbol`: Scale for x-axis, one of :linear or :log10 (default: :linear).
- `y_scale::Symbol`: Scale for y-axis, one of :linear or :log10 (default: :linear).
- `window_function::Function`: Window function to use for spectral estimation (default: DSP.hanning).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.
"""
function plot_components_spectra(
    ica_result::InfoIca,
    dat::ContinuousData,
    comp_indices::Vector{Int}=Int[];
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    window_function::Function=DSP.hanning
)
    # If no components specified, use all components
    if isempty(comp_indices)
        comp_indices = 1:size(ica_result.unmixing, 1)
    end
    
    # Prepare data matrix
    relevant_cols = vcat(ica_result.data_label)
    dat_matrix = permutedims(Matrix(dat.data[:, relevant_cols]))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate component activations
    components = ica_result.unmixing * dat_matrix
    fs = dat.sample_rate
    
    # Create a DataFrame to store all component signals
    component_df = DataFrame(:time => dat.data.time)
    
    # Add each component to the DataFrame with appropriate names
    for comp_idx in comp_indices
        component_name = Symbol("IC_$(comp_idx)")
        component_df[!, component_name] = components[comp_idx, :]
    end
    
    # Create list of component symbols to plot
    components_to_plot = [Symbol("IC_$(comp_idx)") for comp_idx in comp_indices]
    
    # Use the implementation function to create the plot
    fig, ax = _plot_power_spectrum_implementation(
        component_df,
        components_to_plot,
        fs,
        line_freq,
        freq_bandwidth,
        window_size,
        overlap,
        max_freq,
        x_scale,
        y_scale,
        window_function
    )
    
    # Add component variance percentage to the legend
    component_labels = []
    for comp_idx in comp_indices
        variance_pct = ica_result.variance[comp_idx] * 100
        push!(component_labels, @sprintf("IC %d (%.1f%%)", comp_idx, variance_pct))
    end
    
    # Create a new legend with the custom labels
    legend_entries = [LineElement(color=ax.scene.plots[i].color) for i in 1:length(comp_indices)]
    
    # Calculate number of columns based on number of components
    ncols = length(comp_indices) > 10 ? ceil(Int, length(comp_indices) / 10) : 1
    
    # Place legend - use multi-column for many components
    if ncols > 1
        Legend(fig[1, 2], legend_entries, component_labels, "Components", nbanks=ncols)
    else
        Legend(fig[1, 2], legend_entries, component_labels, "Components")
    end
    
    # Set a more descriptive title
    ax.title = "ICA Component Power Spectra"
    
    return fig, ax
end

