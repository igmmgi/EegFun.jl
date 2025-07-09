function _plot_power_spectrum!(fig, ax, df::DataFrame, channels_to_plot::Vector{Symbol}, fs::Real;
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 200.0,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    window_function::Function = DSP.hanning,
)
    # Check if we have enough data and adjust window size if needed
    if size(df, 1) < window_size
        @minimal_warning "Selected region is too short for the specified window size. Adjusting window size..."
        window_size = 2^floor(Int, log2(size(df, 1) / 2))
        if window_size < 32
            @minimal_warning "Selected region is too short for meaningful spectral analysis. Minimum window size of 32 samples required."
            return
        end
    end

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

    # Set axis labels and title
    ax.xlabel = "Frequency (Hz)"
    ax.ylabel = "Power Spectral Density (μV²/Hz)"
    ax.title = "Channel Power Spectrum"
    
    # Apply scale settings
    if x_scale == :log10
        xlims!(ax, (0.1, max_freq))
        ax.xscale = log10
    else
        xlims!(ax, (0, max_freq))
    end
    
    # Calculate max y-value for limits
    max_power = maximum([maximum(psd) for (_, psd) in values(spectra)])
    if y_scale == :log10
        ylims!(ax, (1e-10, max_power))
        ax.yscale = log10
    else
        ylims!(ax, (0, max_power))
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
# DataFrame with auto-extracted sample rate
plot_channel_spectrum(selected_data, channel_selection = channels([:Fp1]))
plot_channel_spectrum(selected_data, channel_selection = channels([:Fp1, :Fp2]))

# DataFrame with explicit sample rate
plot_channel_spectrum(df, 1000, channel_selection = channels([:Fp1]))
plot_channel_spectrum(df, 1000, channel_selection = channels([:Fp1, :Fp2]))

# ContinuousData (auto-extracts sample rate)
plot_channel_spectrum(dat, channel_selection = channels([:Fp1]))
plot_channel_spectrum(dat, channel_selection = channels_not([:M1, :M2]))

# Function-based channel selection
plot_channel_spectrum(dat, channel_selection = channels(x -> startswith.(string.(x), "F")))
plot_channel_spectrum(dat, channel_selection = channels(1:10))

# Sample and channel filtering
plot_channel_spectrum(dat, 
                     sample_selection = samples_not(:is_extreme_value_100),
                     channel_selection = channels_not([:M1, :M2]))

# Custom parameters
plot_channel_spectrum(dat, channel_selection = channels([:Fp1, :Fp2]), 
                     line_freq=60.0, max_freq=100.0, x_scale=:log10)
```
"""
function plot_channel_spectrum(data::DataFrame; 
                              sample_selection::Function = samples(),
                              channel_selection::Function = channels(),
                              kwargs...)
    fs = sample_rate(data)
    return plot_channel_spectrum(data, fs; 
                                sample_selection = sample_selection,
                                channel_selection = channel_selection,
                                kwargs...)
end

function plot_channel_spectrum(data::DataFrame, fs::Real; 
                              sample_selection::Function = samples(),
                              channel_selection::Function = channels(),
                              kwargs...)
    # Get selected channels using the helper function
    all_columns = filter(col -> !(col in [:time, :sample, :triggers]), propertynames(data))
    channel_mask = channel_selection(all_columns)
    selected_channels = all_columns[channel_mask]
    
    # Apply sample selection
    selected_samples = get_selected_samples(data, sample_selection)
    data_subset = data[selected_samples, :]
    
    # Create figure and axis exactly like plot_topoplot
    fig = Figure()
    ax = Axis(fig[1, 1])
    _plot_power_spectrum!(fig, ax, data_subset, selected_channels, fs; kwargs...)
    return fig, ax
end

function plot_channel_spectrum(dat::ContinuousData; 
                              sample_selection::Function = samples(),
                              channel_selection::Function = channels(),
                              display_plot::Bool=true,
                              kwargs...)
    # Get selected channels using the helper function
    selected_channels = get_selected_channels(dat, channel_selection)
    
    # Apply sample selection
    selected_samples = get_selected_samples(dat, sample_selection)
    data_subset = dat.data[selected_samples, :]
    
    # Create figure and axis exactly like plot_topoplot
    fig = Figure()
    ax = Axis(fig[1, 1])
    _plot_power_spectrum!(fig, ax, data_subset, selected_channels, dat.sample_rate; kwargs...)
    if display_plot
        display_figure(fig)
    end
    return fig, ax
end

"""
    plot_ica_component_spectrum(ica_result::InfoIca, dat::ContinuousData, comp_idx::Int; kwargs...)
    plot_ica_component_spectrum(ica_result::InfoIca, dat::ContinuousData, comp_indices::Vector{Int}; kwargs...)
    plot_ica_component_spectrum(ica_result::InfoIca, dat::ContinuousData, comp_indices::AbstractRange{Int}; kwargs...)

Plot power spectrum of ICA component(s).

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.
- `comp_idx::Int`: Single component index to plot.
- `comp_indices::Vector{Int}`: Vector of component indices to plot.
- `comp_indices::AbstractRange{Int}`: Range of component indices to plot.

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

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.

# Examples
```julia
# Plot single component
plot_ica_component_spectrum(ica_result, dat, 1)
plot_ica_component_spectrum(ica_result, dat, 1, sample_selection = samples_not(:is_extreme_value_100))

# Plot multiple components
plot_ica_component_spectrum(ica_result, dat, [1, 3, 5])
plot_ica_component_spectrum(ica_result, dat, 1:10, sample_selection = samples(:epoch_window))

# Plot all components
plot_ica_component_spectrum(ica_result, dat, Int[])
```
"""
function plot_ica_component_spectrum(
    ica_result::InfoIca,
    dat::ContinuousData,
    comp_idx::Int;
    sample_selection::Function = samples(),
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    window_function::Function=DSP.hanning
)
    return plot_ica_component_spectrum(ica_result, dat, [comp_idx];
                                     sample_selection=sample_selection,
                                     line_freq=line_freq,
                                     freq_bandwidth=freq_bandwidth,
                                     window_size=window_size,
                                     overlap=overlap,
                                     max_freq=max_freq,
                                     x_scale=x_scale,
                                     y_scale=y_scale,
                                     window_function=window_function)
end

function plot_ica_component_spectrum(
    ica_result::InfoIca,
    dat::ContinuousData,
    comp_indices::AbstractRange{Int};
    sample_selection::Function = samples(),
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    window_function::Function=DSP.hanning
)
    return plot_ica_component_spectrum(ica_result, dat, collect(comp_indices);
                                     sample_selection=sample_selection,
                                     line_freq=line_freq,
                                     freq_bandwidth=freq_bandwidth,
                                     window_size=window_size,
                                     overlap=overlap,
                                     max_freq=max_freq,
                                     x_scale=x_scale,
                                     y_scale=y_scale,
                                     window_function=window_function)
end

function plot_ica_component_spectrum(
    ica_result::InfoIca,
    dat::ContinuousData,
    comp_indices::Vector{Int};
    sample_selection::Function = samples(),
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    window_size::Int=1024,
    overlap::Real=0.5,
    max_freq::Real=100.0,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    window_function::Function=DSP.hanning,
    display_plot::Bool=true
)
    # If empty vector, use all components
    if isempty(comp_indices)
        comp_indices = 1:size(ica_result.unmixing, 1)
    end
    
    # Apply sample selection
    selected_samples = get_selected_samples(dat, sample_selection)
    
    # Prepare data matrix for selected samples only
    relevant_cols = vcat(ica_result.data_label)
    dat_matrix = permutedims(Matrix(dat.data[selected_samples, relevant_cols]))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate component activations for selected samples
    components = ica_result.unmixing * dat_matrix
    fs = dat.sample_rate
    
    # Create a DataFrame to store all component signals
    component_df = DataFrame(:time => dat.data.time[selected_samples])
    
    # Add each component to the DataFrame with consistent naming
    for comp_idx in comp_indices
        component_name = Symbol("IC_$(comp_idx)")
        component_df[!, component_name] = components[comp_idx, :]
    end
    
    # Create list of component symbols to plot
    components_to_plot = [Symbol("IC_$(comp_idx)") for comp_idx in comp_indices]
    
    # Create figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1])
    _plot_power_spectrum!(fig, ax, component_df, components_to_plot, fs;
                         line_freq = line_freq,
                         freq_bandwidth = freq_bandwidth,
                         window_size = window_size,
                         overlap = overlap,
                         max_freq = max_freq,
                         x_scale = x_scale,
                         y_scale = y_scale,
                         window_function = window_function)
    
        ax.title = "ICA Component Power Spectra"
        
        # Add component variance percentages to the legend
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
        Legend(fig[1, 2], legend_entries, component_labels, "Components", nbanks=ncols)

        if display_plot
            display_figure(fig)
        end

    return fig, ax
end



