###############################################################
# Filter functions 
"""
    _apply_filter!(dat::DataFrame, columns, filter)

Internal helper function to apply a digital filter to specified columns in a DataFrame.
Modifies the data in place.

Arguments:
- `dat`: DataFrame containing the data to filter
- `columns`: Vector of column names to filter
- `filter`: Digital filter object to apply
"""
function _apply_filter!(dat::DataFrame, columns, filter)
    for col in names(dat)
        if col in columns
            dat[:, col] .= filtfilt(filter, dat[:, col])
        end
    end
end

"""
    filter_data!(dat::DataFrame, columns, filter_type, freq, order, filter_method, sample_rate)

Apply a digital filter to specified columns in a DataFrame. Modifies the data in place.

Arguments:
- `dat`: DataFrame containing the data to filter
- `columns`: Vector of column names to filter
- `filter_type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `filter_method`: String specifying filter implementation:
- `sample_rate`: Sampling rate in Hz
"""
function filter_data!(dat::DataFrame, columns, filter_type::String, freq, order::Integer, filter_method::String, sample_rate::Real)
    valid_types = ("hp", "lp", "bp", "bs")
    valid_methods = ("iir", "fir")
    
    # Improved error checking
    if !(filter_type in valid_types)
        throw(ArgumentError("filter_type '$filter_type' must be one of: $valid_types"))
    end
    if !(filter_method in valid_methods)
        throw(ArgumentError("filter_method '$filter_method' must be one of: $valid_methods"))
    end
    if order <= 0
        throw(ArgumentError("filter order must be positive (got $order)"))
    end
    if sample_rate <= 0
        throw(ArgumentError("sample_rate must be positive (got $sample_rate)"))
    end

    # Validate frequency parameters
    if filter_type in ("bp", "bs")
        if !isa(freq, Tuple) || length(freq) != 2
            throw(ArgumentError("freq must be a tuple of (low, high) frequencies for bandpass/bandstop filters"))
        end
        freq_low, freq_high = freq
        if freq_low >= freq_high
            throw(ArgumentError("low frequency ($freq_low) must be less than high frequency ($freq_high)"))
        end
        if freq_low <= 0 || freq_high >= sample_rate/2
            throw(ArgumentError("frequencies must be between 0 and Nyquist frequency ($(sample_rate/2) Hz)"))
        end
    else
        if !(isa(freq, Real) && freq > 0 && freq < sample_rate/2)
            throw(ArgumentError("frequency must be between 0 and Nyquist frequency ($(sample_rate/2) Hz)"))
        end
    end

    # Create filter prototype based on type
    if filter_type == "hp"
        filter_prototype = Highpass(freq)
    elseif filter_type == "lp"
        filter_prototype = Lowpass(freq)
    elseif filter_type == "bp"
        freq_low, freq_high = freq
        filter_prototype = Bandpass(freq_low, freq_high)
    elseif filter_type == "bs"
        freq_low, freq_high = freq
        filter_prototype = Bandstop(freq_low, freq_high)
    end

    # Create filter with chosen method
    if filter_method == "iir"
        filter = digitalfilter(filter_prototype, Butterworth(order); fs=sample_rate)
    else # filter_method == "fir"
        filter = digitalfilter(filter_prototype, FIRWindow(order); fs=sample_rate)
    end
    
    _apply_filter!(dat, columns, filter)
end

"""
    filter_data(dat::DataFrame, columns, type, freq, order, filter_method, sample_rate)

Create a filtered copy of the input DataFrame. Returns a new DataFrame with filtered data.

Arguments:
- `dat`: DataFrame containing the data to filter
- `columns`: Vector of column names to filter
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `sample_rate`: Sampling rate in Hz
- `filter_method`: String specifying filter implementation ("iir" or "fir")

Returns:
- A new DataFrame with filtered data
"""
function filter_data(dat::DataFrame, columns, type, freq, order, filter_method::String, sample_rate::Real)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, columns, type, freq, order, filter_method, sample_rate)
    return dat_out
end

"""
    filter_data!(dat::Union{ContinuousData,ErpData}, type, freq, order)

Apply a filter to ContinuousData or ErpData objects. Modifies the data in place.

Arguments:
- `dat`: ContinuousData or ErpData object
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `filter_method`: String specifying filter implementation ("iir" or "fir")
"""
function filter_data!(dat::Union{ContinuousData,ErpData}, type, freq, order; filter_method::String="iir")
    filter_data!(dat.data, dat.layout.label, type, freq, order, filter_method, dat.sample_rate)
end

"""
    filter_data(dat::Union{ContinuousData,ErpData}, type, freq, order)

Create a filtered copy of ContinuousData or ErpData object.

Arguments:
- `dat`: ContinuousData or ErpData object
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order

Returns:
- A new ContinuousData or ErpData object with filtered data
"""
function filter_data(dat::Union{ContinuousData,ErpData}, type, freq, order; filter_method::String="iir")
    dat_out = deepcopy(dat)
    filter_data!(dat_out.data, dat_out.layout.label, type, freq, order, filter_method, dat_out.sample_rate)
    return dat_out
end



"""
    filter_data(dat::EpochData, type, freq, order, filter_method)

Create a filtered copy of an EpochData object.

Arguments:
- `dat`: EpochData object
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order

Returns:
- A new EpochData object with filtered data
"""
function filter_data(dat::EpochData, type, freq, order; filter_method::String = "iir")
    dat_out = deepcopy(dat)
    for epoch in eachindex(dat_out.data)
        filter_data!(dat_out.data[epoch], dat_out.layout.label, type, freq, order, filter_method, dat_out.sample_rate)
    end
    return dat_out
end


"""
    filter_data(dat::EpochData, type, freq, order, filter_method)

Create a filtered copy of an EpochData object.

Arguments:
- `dat`: EpochData object
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order

Returns:
- A new EpochData object with filtered data
"""
function filter_data!(dat::EpochData, type, freq, order; filter_method::String = "iir")
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], dat.layout.label, type, freq, order, filter_method, dat.sample_rate)
    end
end





"""
    get_filter_characteristics(filter, sample_rate::Real; npoints::Int=1000)

Calculate and return key characteristics of a digital filter.

Arguments:
- `filter`: A digital filter object
- `sample_rate`: Sampling rate in Hz
- `npoints`: Number of frequency points for analysis (default: 1000)

Returns:
A NamedTuple containing:
- `cutoff_freq`: Frequency at -3dB point(s) in Hz
- `passband_ripple`: Maximum ripple in passband in dB
- `stopband_atten`: Minimum attenuation in stopband in dB
- `phase_delay`: Average phase delay in samples
- `group_delay`: Average group delay in samples
- `transition_width`: Width of transition band in Hz
- `filter_type`: Detected filter type ("hp", "lp", "bp", or "bs")
"""
function get_filter_characteristics(filter, sample_rate::Real; npoints::Int=1000)
    # Calculate frequency response
    freqs = range(0, sample_rate/2, length=npoints)
    w = 2π * freqs / sample_rate
    resp = freqresp(filter, w)
    mag_db = 20 * log10.(abs.(resp))
    phase = angle.(resp)
    
    # Find -3dB points (cutoff frequencies)
    cutoff_indices = findall(x -> isapprox(x, -3, atol=0.2), mag_db)
    cutoff_freq = freqs[cutoff_indices]
    
    # Determine filter type based on response shape
    # Compare start, middle, and end of response
    start_response = mag_db[1]
    mid_response = mag_db[floor(Int, npoints/2)]
    end_response = mag_db[end]
    
    # Determine filter type
    if start_response > -3 && end_response < -20
        filter_type = "lp"  # Lowpass: passes low freqs, attenuates high freqs
    elseif start_response < -20 && end_response > -3
        filter_type = "hp"  # Highpass: attenuates low freqs, passes high freqs
    elseif start_response < -20 && end_response < -20 && mid_response > -3
        filter_type = "bp"  # Bandpass: attenuates both low and high freqs
    elseif start_response > -3 && end_response > -3 && mid_response < -20
        filter_type = "bs"  # Bandstop: passes both low and high freqs
    else
        # Default to lowpass if unclear
        filter_type = "lp"
    end
    
    # Calculate passband and stopband characteristics
    passband_mask = if filter_type == "lp"
        freqs .<= (isempty(cutoff_freq) ? sample_rate/4 : cutoff_freq[1])
    elseif filter_type == "hp"
        freqs .>= (isempty(cutoff_freq) ? sample_rate/4 : cutoff_freq[1])
    elseif filter_type == "bp"
        (freqs .>= cutoff_freq[1]) .& (freqs .<= cutoff_freq[end])
    else # bandstop
        (freqs .<= cutoff_freq[1]) .| (freqs .>= cutoff_freq[end])
    end
    
    # Find the stopband (frequencies after -20 dB point)
    stopband_mask = mag_db .<= -20
    
    passband_ripple = maximum(abs.(mag_db[passband_mask])) # dB
    # Calculate mean attenuation in the stopband
    stopband_atten = mean(mag_db[stopband_mask]) # dB
    
    # Calculate delays
    group_delay = -diff(unwrap(phase)) ./ diff(w)
    phase_delay = -phase ./ w
    phase_delay[1] = phase_delay[2] # Handle division by zero at DC
    
    # Calculate transition width
    transition_width = if isempty(cutoff_freq)
        sample_rate/4  # Default to quarter of sampling rate if no cutoff found
    elseif length(cutoff_freq) == 1
        # For lowpass/highpass, find width between -3dB and -20dB points
        db_20_idx = findfirst(x -> x <= -20, mag_db)
        db_20_idx = isnothing(db_20_idx) ? length(mag_db) : db_20_idx
        abs(freqs[db_20_idx] - cutoff_freq[1])
    else
        # For bandpass/bandstop, use the width of both transition regions
        db_20_idx_low = findfirst(x -> x <= -20, mag_db)
        db_20_idx_high = findlast(x -> x <= -20, mag_db)
        db_20_idx_low = isnothing(db_20_idx_low) ? 1 : db_20_idx_low
        db_20_idx_high = isnothing(db_20_idx_high) ? length(mag_db) : db_20_idx_high
        (abs(freqs[db_20_idx_low] - cutoff_freq[1]) + 
         abs(freqs[db_20_idx_high] - cutoff_freq[end])) / 2
    end
    
    return (
        cutoff_freq = cutoff_freq,
        passband_ripple = passband_ripple,
        stopband_atten = stopband_atten,
        phase_delay = mean(phase_delay[passband_mask]),
        group_delay = mean(group_delay[passband_mask[1:end-1]]),
        transition_width = transition_width,
        filter_type = filter_type
    )
end

"""
    print_filter_characteristics(filter, sample_rate::Real; npoints::Int=1000)

Print a formatted summary of filter characteristics.

Arguments:
- `filter`: A digital filter object
- `sample_rate`: Sampling rate in Hz
- `npoints`: Number of frequency points for analysis (default: 1000)
"""
function print_filter_characteristics(filter, sample_rate::Real; npoints::Int=1000)
    chars = get_filter_characteristics(filter, sample_rate; npoints=npoints)
    
    println("Filter Characteristics:")
    println("-------------------------")
    
    # Filter type and cutoff
    println("Type: ", uppercase(chars.filter_type))
    if length(chars.cutoff_freq) == 1
        println("Cutoff (-3 dB): ", round(chars.cutoff_freq[1], digits=1), " Hz")
    elseif !isempty(chars.cutoff_freq)
        println("Cutoff (-3 dB): ", round.(chars.cutoff_freq, digits=1), " Hz")
    end
    
    println("Transition width: ", round(chars.transition_width, digits=1), " Hz")
    println("Stopband attenuation: ", round(chars.stopband_atten, digits=1), " dB")
    
    println("\nTemporal characteristics:")
    # Convert samples to milliseconds using actual sampling rate
    println("Group delay: ", round(chars.group_delay, digits=1), " samples (", 
            round(chars.group_delay * 1000/sample_rate, digits=1), " ms)")
    println("Phase delay: ", round(chars.phase_delay, digits=1), " samples (",
            round(chars.phase_delay * 1000/sample_rate, digits=1), " ms)")
end

"""
    plot_filter_response(filter, fs; npoints::Int=2000)

Create a plot showing the frequency response of a digital filter using Makie,
styled similarly to MNE-Python's filter plotting.

Arguments:
- `filter`: A digital filter object
- `fs`: Sampling frequency in Hz
- `npoints`: Number of frequency points for analysis (default: 2000)

Returns:
- A Makie figure showing magnitude response
"""
function plot_filter_response(filter, fs::Real; npoints::Int=2000)
    # Linear spacing from 0 to Nyquist
    freqs = range(0, fs/2, length=npoints)
    
    # Calculate frequency response
    w = 2π * freqs / fs
    resp = freqresp(filter, w)
    
    # Calculate magnitude response in dB
    mag_db = 20 * log10.(abs.(resp))
    
    # Create figure
    fig = Figure(size=(800, 400))
    
    # Magnitude response plot
    ax = Axis(fig[1, 1],
        xlabel="Frequency (Hz)",
        ylabel="Magnitude (dB)",
        title="Filter frequency response",
        xgridvisible=true,
        ygridvisible=true
    )
    
    # Set y-axis limits to match MNE style (-60 to 10 dB)
    ylims!(ax, -60, 10)
    
    # Set x-axis limits explicitly
    xlims!(ax, 0, fs/2)
    
    # Plot the frequency response
    lines!(ax, freqs, mag_db, linewidth=2, color=:blue)
    
    # Add horizontal lines at 0, -3, -20, and -40 dB
    hlines!(ax, 0, color=:gray, linestyle=:dash, alpha=0.5)
    hlines!(ax, -3, color=:gray, linestyle=:dash, alpha=0.5)
    hlines!(ax, -20, color=:gray, linestyle=:dash, alpha=0.5)
    hlines!(ax, -40, color=:gray, linestyle=:dash, alpha=0.5)
    
    # Add text labels for the lines
    text!(ax, fs/2 * 1.1, 0, text="0 dB", align=(:left, :center))
    text!(ax, fs/2 * 1.1, -3, text="-3 dB", align=(:left, :center))
    text!(ax, fs/2 * 1.1, -20, text="-20 dB", align=(:left, :center))
    text!(ax, fs/2 * 1.1, -40, text="-40 dB", align=(:left, :center))
    
    # Find and mark the -3dB points
    cutoff_indices = findall(x -> isapprox(x, -3, atol=0.2), mag_db)
    if !isempty(cutoff_indices)
        scatter!(ax, freqs[cutoff_indices], mag_db[cutoff_indices],
            color=:red,
            markersize=10
        )
        
        # Add vertical lines at cutoff frequencies
        for idx in cutoff_indices
            vlines!(ax, freqs[idx], color=:red, linestyle=:dash, alpha=0.5)
            # Add frequency labels
            text!(ax, freqs[idx], -55, 
                text="$(round(Int, freqs[idx])) Hz",
                align=(:center, :top)
            )
        end
    end
    
    return fig
end


# Set up filter parameters
sample_rate = 512    # Hz
cutoff_freq = 50      # Hz
filter_order = 12      # Filter order (affects steepness of roll-off)

# Create the filter
responsetype = Lowpass(cutoff_freq)
designmethod = Butterworth(filter_order)
my_filter = digitalfilter(responsetype, designmethod, fs=sample_rate)

# Now you can:
# 1. Plot the filter response
print_filter_characteristics(my_filter, sample_rate)
plot_filter_response(my_filter, sample_rate)




