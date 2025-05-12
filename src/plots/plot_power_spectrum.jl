
function _plot_power_spectrum_implementation(
    df::DataFrame,
    channels_to_plot::Vector{Symbol},
    fs::Real,
    line_freq::Real,
    freq_bandwidth::Real,
    window_size::Int,
    overlap::Real,
    max_freq::Real,
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
        pgram = DSP.welch_pgram(signal, window_size, noverlap; fs = fs)
        spectra[ch] = (DSP.freq(pgram), DSP.power(pgram))
    end

    # Create figure with appropriate title
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Frequency (Hz)", ylabel = "Power Spectral Density (μV²/Hz)", title = "Channel Power Spectrum")

    # Plot spectra for all channels
    for ch in channels_to_plot
        freqs, psd = spectra[ch]
        if length(channels_to_plot) > 1
            lines!(ax, freqs, psd, label = string(ch))
        else
            lines!(ax, freqs, psd, color = :black)
        end
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
                label = h == 1 && length(channels_to_plot) == 1 ? "Line Frequency" : "",
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

            # Add frequency label
            text!(ax, freq, max_power, text = "$freq Hz", color = :red, align = (:center, :bottom))
        end
    end

    # Set x-axis limits
    xlims!(ax, (0, max_freq))

    # Add legend for multiple channels
    if length(channels_to_plot) > 1
        axislegend(ax, position = (1.0, 1.0))
    end

    return fig, ax
end

"""
    plot_channel_spectrum(df::DataFrame, channel::Union{Symbol,Vector{Symbol}},
                         fs::Real;
                         line_freq::Real=50.0,
                         freq_bandwidth::Real=1.0,
                         window_size::Int=1024,
                         overlap::Real=0.5,
                         max_freq::Real=100.0,
                         display_plot::Bool=false)

Plot the power spectrum of one or more channels from a DataFrame.

# Arguments
- `df::DataFrame`: The data frame containing channel data.
- `channel`: Channel(s) to plot (must be column names in the DataFrame).
- `fs::Real`: Sampling frequency in Hz.

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
- `include_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in the DataFrame marking samples to include.
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in the DataFrame marking samples to exclude.
- `display_plot::Bool`: Whether to display the plot immediately (default: false).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.
"""
function plot_channel_spectrum(
    df::DataFrame,
    channel::Union{Symbol,Vector{Symbol}};
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 100.0,
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
        sample_rate(df),
        line_freq,
        freq_bandwidth,
        window_size,
        overlap,
        max_freq,
    )
end

"""
    plot_channel_spectrum(dat::ContinuousData, channel::Union{Symbol,Vector{Symbol},Nothing}=nothing;
                         line_freq::Real=50.0,
                         freq_bandwidth::Real=1.0,
                         window_size::Int=1024,
                         overlap::Real=0.5,
                         max_freq::Real=100.0)

Plot the power spectrum of one or more EEG channels from ContinuousData.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `channel`: The channel(s) to plot. Can be:
  - `Symbol`: Single channel name
  - `Vector{Symbol}`: Multiple channel names
  - `nothing`: All channels in the layout (default)

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
"""
function plot_channel_spectrum( dat::ContinuousData; kwargs...)
    return plot_channel_spectrum(dat.data, kwargs...)
end

"""
    plot_selected_spectrum(selected_data::DataFrame, channel::Union{Symbol,Vector{Symbol}};
                          line_freq::Real=50.0,
                          freq_bandwidth::Real=1.0,
                          window_size::Int=1024,
                          overlap::Real=0.5,
                          max_freq::Real=100.0)

Plot the power spectrum of a selected time region for specific channel(s).

# Arguments
- `selected_data::DataFrame`: The selected time region data.
- `channel`: Channel(s) to plot (must be column names in selected_data).

# Keyword Arguments
- `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
- `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
- `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
- `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
- `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).

# Returns
- `fig::Figure`: The Makie Figure containing the power spectrum plot.
- `ax::Axis`: The axis containing the plot.
"""
function plot_selected_spectrum(
    selected_data::DataFrame,
    channel::Union{Symbol,Vector{Symbol}};
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    window_size::Int = 1024,
    overlap::Real = 0.5,
    max_freq::Real = 100.0,
)
    # Delegate to the main plotting function with display_plot=true
    fig, ax = plot_channel_spectrum(
        selected_data,
        channel;
        fs = fs,
        line_freq = line_freq,
        freq_bandwidth = freq_bandwidth,
        window_size = window_size,
        overlap = overlap,
        max_freq = max_freq,
    )

    return fig, ax
end






# """
#     plot_component_spectrum(ica_result::InfoIca, dat::ContinuousData, comp_idx::Int;
#                           exclude_samples::Union{Nothing,Vector{Symbol}} = nothing,
#                           line_freq::Real=50.0,
#                           freq_bandwidth::Real=1.0,
#                           window_size::Int=1024,
#                           overlap::Real=0.5,
#                           max_freq::Real=100.0)
# 
# Plot the power spectrum of a specific ICA component.
# 
# # Arguments
# - `ica_result::InfoIca`: The ICA result object.
# - `dat::ContinuousData`: The continuous data.
# - `comp_idx::Int`: Index of the component to plot.
# 
# # Keyword Arguments
# - `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude.
# - `line_freq::Real`: Line frequency in Hz to highlight (default: 50.0).
# - `freq_bandwidth::Real`: Bandwidth around line frequency to highlight (default: 1.0 Hz).
# - `window_size::Int`: Size of the FFT window for spectral estimation (default: 1024).
# - `overlap::Real`: Overlap between windows for Welch's method (default: 0.5).
# - `max_freq::Real`: Maximum frequency to display in Hz (default: 100.0).
# 
# # Returns
# - `fig::Figure`: The Makie Figure containing the power spectrum plot.
# """
# function plot_component_spectrum(
#     ica_result::InfoIca,
#     dat::ContinuousData,
#     comp_idx::Int;
#     exclude_samples::Union{Nothing,Vector{Symbol}} = nothing,
#     line_freq::Real = 50.0,
#     freq_bandwidth::Real = 1.0,
#     window_size::Int = 1024,
#     overlap::Real = 0.5,
#     max_freq::Real = 100.0,
# )
#     # Get samples to use
#     samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
#     if isempty(samples_to_use)
#         @warn "No samples remaining after applying exclude criteria. Cannot plot component spectrum."
#         return Figure()
#     end
# 
#     # Prepare data matrix for valid samples
#     relevant_cols = vcat(ica_result.data_label)
#     data_subset_df = dat.data[samples_to_use, relevant_cols]
#     dat_matrix = permutedims(Matrix(data_subset_df))
#     dat_matrix .-= mean(dat_matrix, dims = 2)
#     dat_matrix ./= ica_result.scale
# 
#     # Calculate component activation
#     components = ica_result.unmixing * dat_matrix
#     fs = dat.sample_rate
# 
#     # Get the component time series and create a temporary DataFrame
#     signal = components[comp_idx, :]
#     component_name = Symbol("Component_$(comp_idx)")
#     component_df = DataFrame(component_name => signal, :time => dat.data.time[samples_to_use])
# 
#     # Create a custom title
#     title = "Power Spectrum of Component $comp_idx"
#     
#     # Call the implementation function directly
#     fig, ax = _plot_power_spectrum_implementation(
#         component_df,
#         [component_name],
#         fs,
#         line_freq,
#         freq_bandwidth,
#         window_size,
#         overlap,
#         max_freq,
#     )
# 
#     return fig
# end

