"""
    create_filter(filter_type::String, filter_method::String, filter_freq::Real, sample_rate::Real; 
                 order::Integer = 3, transition_width::Real = 0.25, plot_filter::Bool = false, print_filter::Bool = false)

Create a digital filter object for the specified parameters.

# Arguments
- `filter_type`: String specifying filter type ("hp"=highpass, "lp"=lowpass)
- `filter_method`: String specifying filter implementation ("iir" or "fir")
- `filter_freq`: Cutoff frequency in Hz
- `sample_rate`: Sampling rate in Hz
- `order`: Filter order for IIR filters (default: 3)
- `transition_width`: Relative width of transition band as fraction of cutoff (default: 0.25)
- `plot_filter`: Boolean to plot frequency response (default: false)
- `print_filter`: Boolean to print filter characteristics (default: false)

# Returns
- Digital filter object that can be applied to data

# Notes
- Validates all input parameters
- Creates appropriate filter prototype based on type and method
- Returns a filter object that can be reused for multiple datasets
- If plot_filter is true, displays the filter frequency response
- If print_filter is true, prints detailed filter characteristics
"""
function create_filter(
    filter_type::String,
    filter_method::String,
    filter_freq::Real,
    sample_rate::Real;
    order::Integer = 3,
    transition_width::Real = 0.25,
    plot_filter::Bool = false,
    print_filter::Bool = false,
)
    # Input validation
    valid_types = ("hp", "lp")
    if !(filter_type in valid_types)
        @minimal_error "filter_type '$filter_type' must be one of: $valid_types"
    end
    
    valid_methods = ("iir", "fir")
    if !(filter_method in valid_methods)
        @minimal_error "filter_method '$filter_method' must be one of: $valid_methods"
    end

    if order <= 0
        @minimal_error "filter order must be positive: $order"
    end

    if sample_rate <= 0
        @minimal_error "sample_rate must be positive: $sample_rate"
    end

    if !(filter_freq > 0 && filter_freq < sample_rate / 2)
        @minimal_error "filter_freq must be between 0 and Nyquist frequency ($(sample_rate/2) Hz): $filter_freq"
    end

    # Create filter prototype based on type
    transition_band = transition_width * filter_freq
    if filter_type == "hp"
        filter_prototype = Highpass(filter_freq - (transition_band/2))  # Subtract for highpass
    elseif filter_type == "lp"
        filter_prototype = Lowpass(filter_freq + (transition_band/2))   # Add for lowpass
    end

    # Create filter with chosen method
    if filter_method == "iir"
        filter = digitalfilter(filter_prototype, Butterworth(order); fs = sample_rate)
    elseif filter_method == "fir"
        # Calculate number of taps ensuring it's a power of 2 for optimal FFT performance
        n_taps = Int(ceil(FIR_TAP_MULTIPLIER * sample_rate / transition_band))
        # Ensure odd number of taps for FIR filter
        if n_taps % 2 == 0  
            n_taps += 1     
        end
        # Ensure n_taps is not too small
        n_taps = max(n_taps, MIN_FIR_TAPS)
        # For FFTW 1.9.0 compatibility, ensure n_taps is a power of 2
        n_taps = nextpow(FFTW_COMPATIBILITY_POWER, n_taps)
        if n_taps % 2 == 0  
            n_taps += 1     
        end
        filter = digitalfilter(filter_prototype, FIRWindow(hamming(n_taps)); fs = sample_rate)
    end

    # Print filter characteristics if requested
    if print_filter
        print_filter_characteristics(filter, sample_rate, filter_freq, transition_band)
    end

    # Plot filter response if requested
    if plot_filter
        plot_filter_response(filter, sample_rate, filter_freq, transition_band)
    end

    return filter
end

"""
    _apply_filter!(dat::DataFrame, channels, filter; filter_func = filtfilt)

Internal helper function to apply a digital filter to specified columns in a DataFrame.
Modifies the data in place.

Arguments:
- `dat`: DataFrame containing the data to filter
- `channels`: Vector of column names to filter
- `filter`: Digital filter object to apply
- `filter_func`: Filtering function to use (default: filtfilt, for two-pass filtering). Use `filt` for one-pass filtering.
"""
function _apply_filter!(dat::DataFrame, channels, filter; filter_func = filtfilt)
    @inbounds for channel in channels
        @views dat[:, channel] .= filter_func(filter, dat[:, channel])
    end
end

"""
    filter_data!(dat::DataFrame, filter, sample_rate::Real; 
                channel_selection::Function = channels(), filter_func=filtfilt)

Apply a pre-generated digital filter to selected channels in a DataFrame. Modifies the data in place.

Arguments:
- `dat`: DataFrame containing the data to filter
- `filter`: Pre-generated digital filter object
- `sample_rate`: Sampling rate in Hz (for filter characteristics display)
- `channel_selection`: Channel selection predicate (default: channels() - all channels)
- `filter_func`: Filtering function to use (default: filtfilt, for two-pass filtering). Use `filt` for one-pass filtering.
"""
function filter_data!(
    dat::DataFrame,
    filter,
    sample_rate::Real;
    channel_selection::Function = channels(),
    filter_func = filtfilt,
)
    selected_channels = get_selected_channels(dat, channel_selection)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for filtering"
        return
    end
    @info "filter_data! applying filter to $(length(selected_channels)) channels"
    
    # Apply filter
    _apply_filter!(dat, selected_channels, filter; filter_func = filter_func)
end

function _update_filter_info!(dat::EegData, filter_type::String, filter_freq::Real)
    if filter_type == "hp"
        dat.analysis_info.hp_filter = filter_freq
    elseif filter_type == "lp"
        dat.analysis_info.lp_filter = filter_freq
    end
end

"""
    filter_data!(dat::EegData, filter_type, filter_method, filter_freq; kwargs...)

Apply a digital filter to EEG data. Modifies the data in place.

Arguments:
- `dat`: EegData object (ContinuousData, ErpData, or EpochData)
- `filter_type`: String specifying filter type ("hp"=highpass, "lp"=lowpass)
- `filter_method`: String specifying filter implementation ("iir" or "fir")
- `filter_freq`: Cutoff frequency in Hz
- `order`: Filter order for IIR filters (default: 3)
- `transition_width`: Relative width of transition band as fraction of cutoff (default: 0.25)
- `channel_selection`: Channel selection predicate (default: channels() - all channels)
- `filter_func`: Filtering function to use (default: filtfilt, for two-pass filtering). Use `filt` for one-pass filtering.
- `plot_filter`: Boolean to plot frequency response (default: false)
- `print_filter`: Boolean to print filter characteristics (default: false)
"""
function filter_data!(
    dat::EegData,
    filter_type,
    filter_method,
    filter_freq;
    order = 3,
    transition_width = 0.25,
    channel_selection::Function = channels(),
    filter_func = filtfilt,
    plot_filter = false,
    print_filter = false,
)
    _update_filter_info!(dat, filter_type, filter_freq)
    
    # Create filter once
    filter = create_filter(filter_type, filter_method, filter_freq, dat.sample_rate; 
                         order, transition_width, plot_filter, print_filter)
    
    # Apply filter (dispatch handles DataFrame vs Vector{DataFrame})
    filter_data!(dat.data, filter, dat.sample_rate; channel_selection, filter_func)
end

# generates all non-mutating versions
@add_nonmutating filter_data!













############################################################
"""
    get_filter_characteristics(filter, sample_rate::Real, transition_width::Real; npoints::Int=1000)

Calculate and return key characteristics of a digital filter.

Arguments:
- `filter`: A digital filter object (FIR coefficients or DSP.jl filter)
- `sample_rate`: Sampling rate in Hz
- `transition_width`: Width of transition band in Hz
- `npoints`: Number of frequency points for analysis (default: 1000)

Returns:
A NamedTuple containing:
- `cutoff_freq_3db`: Frequency at -3dB point(s) in Hz
- `cutoff_freq_6db`: Frequency at -6dB point(s) in Hz
- `passband_ripple`: Maximum ripple in passband in dB
- `stopband_atten`: Mean attenuation in stopband in dB
- `phase_delay`: Mean phase delay in passband (samples)
- `group_delay`: Mean group delay in passband (samples)
- `transition_width`: Width of transition band in Hz
- `filter_type`: Detected filter type ("hp" or "lp")
"""
function get_filter_characteristics(filter, sample_rate::Real, transition_width::Real; npoints::Int=1000)
    # Calculate frequency response
    freqs = exp10.(range(log10(0.01), log10(sample_rate/2), length=npoints))  # log spacing
    # freqs = range(0, sample_rate/2, length=npoints) # linear spacing
    w = 2π * freqs / sample_rate
    
    # Get frequency response based on filter type
    if filter isa Vector  # FIR filter coefficients
        n = 0:(length(filter)-1)
        resp = [sum(filter .* exp.(-im * w_k * n)) for w_k in w]
    else  # Other filter types
        resp = freqresp(filter, w)
    end   

    mag_db = 20 * log10.(abs.(resp))
    phase = angle.(resp)
    
    # Find -3dB and -6dB points using closest points method
    diffs_3db = abs.(mag_db .- (-3))
    diffs_6db = abs.(mag_db .- (-6))
    
    # Find local minima in the differences to detect crossings
    function find_crossings(diffs, threshold=0.5)
        crossings = Int[]
        for i in 2:(length(diffs)-1)
            if diffs[i] < threshold && diffs[i] <= diffs[i-1] && diffs[i] <= diffs[i+1]
                push!(crossings, i)
            end
        end
        return crossings
    end
    
    # Find crossings with more generous threshold
    crossings_3db = find_crossings(diffs_3db, 0.5)  # Within 0.5 dB of -3dB point
    crossings_6db = find_crossings(diffs_6db, 0.5)  # Within 0.5 dB of -6dB point
    
    # If no crossings found, use points closest to target
    if isempty(crossings_3db)
        crossings_3db = [argmin(diffs_3db)]
    end
    if isempty(crossings_6db)
        crossings_6db = [argmin(diffs_6db)]
    end
    
    # Convert to frequencies
    cutoff_freq_3db = freqs[crossings_3db]
    cutoff_freq_6db = freqs[crossings_6db]
    
    # Use -3dB points for filter type detection and other calculations
    cutoff_freq = cutoff_freq_3db
    
    # Determine filter type based on response shape
    start_response = mag_db[2]  # Use second point to avoid DC issues
    end_response = mag_db[end]

    # Determine filter type and masks
    if start_response > -3 && end_response < -20
        filter_type = "lp"  # Lowpass
        passband_mask = freqs .<= cutoff_freq[1]
        stopband_mask = freqs .>= (cutoff_freq[1] + transition_width)
    else  # highpass
        filter_type = "hp"
        passband_mask = freqs .>= cutoff_freq[1]
        stopband_mask = freqs .<= (cutoff_freq[1] - transition_width)
    end
    
    # Calculate passband ripple and stopband attenuation
    passband_ripple = maximum(abs.(mag_db[passband_mask])) # dB
    # Calculate mean stopband attenuation, replacing -Inf with -100 dB
    stopband_db = mag_db[stopband_mask]
    stopband_db = replace(stopband_db, -Inf => -100.0)
    stopband_atten = mean(stopband_db)
    
    # Calculate delays
    group_delay = -diff(unwrap(phase)) ./ diff(w)
    phase_delay = -phase ./ w
    phase_delay[1] = phase_delay[2] # Handle division by zero at DC

    return (
        cutoff_freq_3db = cutoff_freq_3db,
        cutoff_freq_6db = cutoff_freq_6db,
        passband_ripple = passband_ripple,
        stopband_atten = stopband_atten,
        phase_delay = mean(phase_delay[passband_mask]),
        group_delay = mean(group_delay[passband_mask[1:end-1]]),
        transition_width = transition_width,
        filter_type = filter_type
    )
end

"""
    print_filter_characteristics(filter, sample_rate::<:Real, filter_freq::Union{<:Real,Tuple}, transition_width::<:Real; npoints::Int=1000)

Print a formatted summary of filter characteristics.

Arguments:
- `filter`: A digital filter object
- `sample_rate`: Sampling rate in Hz
- `filter_freq`: Cutoff frequency in Hz
- `transition_width`: Width of transition band in Hz
- `npoints`: Number of frequency points for analysis (default: 1000)
"""
function print_filter_characteristics(filter, sample_rate::Real, filter_freq::Union{Real,Tuple}, transition_width::Real; npoints::Int=1000)

    chars = get_filter_characteristics(filter, sample_rate, transition_width; npoints=npoints)
    
    @info "--------------------------------------------------------------------------------" 
    @info "Filter Characteristics: $(uppercase(chars.filter_type)) $(filter_freq) Hz" 
    @info "Cutoff: -3 dB = $(round.(chars.cutoff_freq_3db, digits=2)) Hz, -6 dB = $(round.(chars.cutoff_freq_6db, digits=2)) Hz" 
    @info "Transition width: $(round(chars.transition_width, digits=1)) Hz" 
    @info "Stopband attenuation: $(round(chars.stopband_atten, digits=1)) dB" 
    @info "Single-pass characteristics: Group delay = $(round(chars.group_delay * 1000/sample_rate, digits=1)) ms, Phase delay = $(round(chars.phase_delay * 1000/sample_rate, digits=1)) ms " 

end

"""
    plot_filter_response(filter, sample_rate::Real, filter_freq::Real, transition_band::Real; 
                        ylimit::Tuple=(-100, 5), xlimit::Union{Nothing,Tuple{Real,Real}}=nothing)

Plot the frequency response of a digital filter with ideal response overlay.

Arguments:
- `filter`: A digital filter object (FIR coefficients or DSP.jl filter)
- `sample_rate`: Sampling rate in Hz
- `filter_freq`: Cutoff frequency in Hz
- `transition_band`: Width of transition band in Hz
- `ylimit`: Y-axis limits in dB as (min, max) tuple (default: (-100, 5))
- `xlimit`: Optional X-axis limits in Hz as (min, max) tuple. If nothing, automatically determined.

Returns:
- `fig`: Makie Figure object
- `ax`: Makie Axis object
"""
function plot_filter_response(
    filter,
    sample_rate::Real,
    filter_freq::Real,
    transition_band::Real;
    ylimit::Tuple = (-100, 5),
    xlimit::Union{Nothing,Tuple{Real,Real}} = nothing,
)
    # Determine x-axis limits based on filter type and cutoff
    if isnothing(xlimit)
        if filter_freq < 2
            xlimit = (0, filter_freq * 10)
        else
            xlimit = (0, sample_rate / 2)
        end
    end

    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel="Frequency (Hz)",
        ylabel="Magnitude (dB)",
        title="Filter Frequency Response",
        titlesize=24,
        xlabelsize=22,
        ylabelsize=22,
        xticklabelsize=20,
        yticklabelsize=20,
        xscale=Makie.Symlog10(10.0),
        limits=(xlimit, ylimit)  
    )

    # Simple logarithmic frequency spacing
    n_points = 2000
    freqs = exp10.(range(log10(0.01), log10(sample_rate/2), length=n_points))
    freqs = [0.0; freqs]  # Add DC point
    
    w = 2π * freqs / sample_rate
    
    # Get frequency response
    if filter isa Vector  # FIR filter coefficients
        n = 0:(length(filter)-1)
        resp = [sum(filter .* exp.(-im * w_k * n)) for w_k in w]
    else  # Other filter types
        resp = freqresp(filter, w)
    end
    
    mag_db = 20 * log10.(abs.(resp))
    
    # Get filter type from response
    start_response = mag_db[2]  # Use second point to avoid DC issues
    end_response = mag_db[end]
    
    # Determine filter type
    filter_type = "hp"
    if start_response > -3 && end_response < -20
        filter_type = "lp"
    end

    # Calculate actual stopband attenuation
    if filter_type == "lp"
        stopband_mask = freqs .>= (filter_freq + transition_band)
    elseif filter_type == "hp" 
        stopband_mask = freqs .<= (filter_freq - transition_band)
    end
    
    stopband_db = mag_db[stopband_mask]
    stopband_db = replace(stopband_db, -Inf => -100.0)
    actual_stopband_atten = mean(stopband_db)
    actual_stopband_linear = 10^(actual_stopband_atten/20)
    
    # Calculate ideal response
    ideal_response = zeros(length(freqs))
    for (i, f) in enumerate(freqs)
        if filter_type == "lp" 
            if f <= filter_freq
                ideal_response[i] = 1.0
            elseif f >= filter_freq + transition_band
                ideal_response[i] = actual_stopband_linear
            else
                ideal_response[i] =
                    1.0 * (filter_freq + transition_band - f) / transition_band +
                    actual_stopband_linear * (f - filter_freq) / transition_band
            end
        elseif filter_type == "hp" 
            if f >= filter_freq
                ideal_response[i] = 1.0
            elseif f <= filter_freq - transition_band
                ideal_response[i] = actual_stopband_linear
            else
                ideal_response[i] =
                    1.0 * (f - (filter_freq - transition_band)) / transition_band +
                    actual_stopband_linear * (filter_freq - f) / transition_band
            end
        end
    end

  # Define transition regions first
    if filter_type == "lp"
        f_s = filter_freq + transition_band  # Add for lowpass
        transition_regions = [(filter_freq, f_s)]
    elseif filter_type == "hp"
        f_s = filter_freq - transition_band  # Subtract for highpass
        transition_regions = [(f_s, filter_freq)]
    end

    for (start_f, end_f) in transition_regions
        vspan!(ax, start_f, end_f, color=(:gray, 0.2))
        vlines!(ax, [start_f, end_f], color=:gray, linestyle=:dash)
    end
    
    # Plot responses
    lines!(ax, freqs, mag_db, label="Actual", color=:black, linewidth=4)
    ideal_mag_db = 20 * log10.(ideal_response)
    lines!(ax, freqs, ideal_mag_db, color=:green, linestyle=:dash, label="Ideal",linewidth=3)
    
    # Add reference lines
    hlines!(ax, [-3, -6], color=:gray, linestyle=:dash)
    text!(ax, sample_rate/2, -3, text="-3 dB", align=(:right, :center), fontsize=22)
    text!(ax, sample_rate/2, -6, text="-6 dB", align=(:right, :center), fontsize=22)

    # X-axis ticks based on xlimits
    if xlimit[2] <= 5  # For low frequency (typically highpass) plots
        xticks = [0, 0.1, 0.2, 0.5, 1, 2, filter_freq, f_s, 5]
    else  # For full range plots (typically lowpass)
        xticks = [0, 1, 2, 5, 10, 20, 50, filter_freq, f_s, 100, sample_rate/2]
    end
    
    # Filter out ticks beyond xlimit
    xticks = [x for x in xticks if x >= xlimit[1] && x <= xlimit[2]]
    ax.xticks = (xticks, string.(round.(xticks, digits=1)))    

    legend_position = filter_type == "lp" ? :lb : :rb
    axislegend(;position=legend_position, labelsize=32)
    display(fig)
    
    return fig, ax
end
