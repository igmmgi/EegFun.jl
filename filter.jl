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
    filter_data!(dat::DataFrame, columns, filter_type, freq, order, sample_rate; filter_method="iir")

Apply a digital filter to specified columns in a DataFrame. Modifies the data in place.

Arguments:
- `dat`: DataFrame containing the data to filter
- `columns`: Vector of column names to filter
- `filter_type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `sample_rate`: Sampling rate in Hz
- `filter_method`: String specifying filter implementation:
    - "iir": IIR Butterworth filter (default)
    - "fir": FIR filter with Hamming window

Note: FIR filters typically need a higher order than IIR filters for similar performance.
"""
function filter_data!(dat::DataFrame, columns, filter_type::String, freq, order::Integer, 
                     sample_rate::Real; filter_method::String="iir")
    # Input validation
    valid_types = ("hp", "lp", "bp", "bs")
    valid_methods = ("iir", "fir")
    
    if !(filter_type in valid_types)
        throw(ArgumentError("filter_type must be one of: $valid_types"))
    end
    if !(filter_method in valid_methods)
        throw(ArgumentError("filter_method must be one of: $valid_methods"))
    end

    # Create filter prototype based on type
    prototype = if filter_type == "hp"
        Highpass(freq)
    elseif filter_type == "lp"
        Lowpass(freq)
    elseif filter_type == "bp"
        freq_low, freq_high = freq
        Bandpass(freq_low, freq_high)
    elseif filter_type == "bs"
        freq_low, freq_high = freq
        Bandstop(freq_low, freq_high)
    end

    # Create filter with chosen method
    filter = if filter_method == "iir"
        digitalfilter(prototype, Butterworth(order); fs=sample_rate)
    else # fir
        digitalfilter(prototype, FIRWindow(order); fs=sample_rate)
    end
    
    _apply_filter!(dat, columns, filter)
end

"""
    filter_data(dat::DataFrame, columns, type, freq, order, sample_rate; filter_method="iir")

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
function filter_data(dat::DataFrame, columns, type, freq, order, sample_rate; filter_method::String="iir")
    dat_out = deepcopy(dat)
    filter_data!(dat_out, columns, type, freq, order, sample_rate; filter_method=filter_method)
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
"""
function filter_data!(dat::Union{ContinuousData,ErpData}, type, freq, order)
    filter_data!(dat.data, dat.layout.label, type, freq, order, dat.sample_rate)
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
function filter_data(dat::Union{ContinuousData,ErpData}, type, freq, order)
    dat_out = deepcopy(dat)
    filter_data!(dat_out.data, dat_out.layout.label, type, freq, order, dat_out.sample_rate)
    return dat_out
end

"""
    filter_data!(dat::EpochData, args...)

Apply filtering to all epochs in an EpochData object. Modifies the data in place.
Accepts the same arguments as the DataFrame version after the first argument.

Arguments:
- `dat`: EpochData object
- `args...`: Additional arguments passed to the DataFrame filter_data! function
"""
function filter_data!(dat::EpochData, args...)
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], args...)
    end
end

"""
    filter_data(dat::EpochData, type, freq, order)

Create a filtered copy of an EpochData object.

Arguments:
- `dat`: EpochData object
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order

Returns:
- A new EpochData object with filtered data
"""
function filter_data(dat::EpochData, type, freq, order)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, type, freq, order)
    return dat_out
end

"""
    filter_data!(dat::EpochData, columns, type, freq, order, sample_rate)

Apply filtering to specific columns of all epochs in an EpochData object.
Modifies the data in place.

Arguments:
- `dat`: EpochData object
- `columns`: Vector of column names to filter
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `sample_rate`: Sampling rate in Hz
"""
function filter_data!(dat::EpochData, columns, type, freq, order, sample_rate)
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], columns, type, freq, order, sample_rate)
    end
end

"""
    filter_data(dat::EpochData, columns, type, freq, order, sample_rate)

Create a filtered copy of an EpochData object, filtering specific columns.

Arguments:
- `dat`: EpochData object
- `columns`: Vector of column names to filter
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `sample_rate`: Sampling rate in Hz

Returns:
- A new EpochData object with filtered data
"""
function filter_data(dat::EpochData, columns, type, freq, order, sample_rate)
    dat_out = deepcopy(dat)
    for epoch in eachindex(dat.data)
        filter_data!(dat.data[epoch], columns, type, freq, order, sample_rate)
    end
    return dat_out
end

"""
    filter_data(dat::Vector{DataFrame}, columns, type, freq, order, sample_rate)

Create a filtered copy of a vector of DataFrames.

Arguments:
- `dat`: Vector of DataFrames
- `columns`: Vector of column names to filter
- `type`: String specifying filter type ("hp"=highpass, "lp"=lowpass, "bp"=bandpass, "bs"=bandstop)
- `freq`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `order`: Filter order
- `sample_rate`: Sampling rate in Hz

Returns:
- A new Vector of DataFrames with filtered data
"""
function filter_data(dat::Vector{DataFrame}, columns, type, freq, order, sample_rate)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, columns, type, freq, order, sample_rate)
    return dat_out
end

"""
    get_filter_characteristics(filter; npoints::Int=1000)

Calculate and return key characteristics of a digital filter.

Arguments:
- `filter`: A digital filter object
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
function get_filter_characteristics(filter; npoints::Int=1000)
    fs = samplingfreq(filter)
    freqs = range(0, fs/2, length=npoints)
    resp = freqresp(filter, freqs)
    mag_db = 20 * log10.(abs.(resp))
    phase = angle.(resp)
    
    # Find -3dB points (cutoff frequencies)
    cutoff_indices = findall(x -> isapprox(x, -3, atol=0.1), mag_db)
    cutoff_freq = freqs[cutoff_indices]
    
    # Determine filter type based on response shape
    filter_type = if length(cutoff_freq) == 1
        mag_db[1] < -20 ? "hp" : "lp"
    else
        mag_db[floor(Int, npoints/2)] < -20 ? "bs" : "bp"
    end
    
    # Calculate passband and stopband characteristics
    passband_mask = if filter_type == "lp"
        freqs .<= cutoff_freq[1]
    elseif filter_type == "hp"
        freqs .>= cutoff_freq[1]
    elseif filter_type == "bp"
        (freqs .>= cutoff_freq[1]) .& (freqs .<= cutoff_freq[end])
    else # bandstop
        (freqs .<= cutoff_freq[1]) .| (freqs .>= cutoff_freq[end])
    end
    
    passband_ripple = maximum(abs.(mag_db[passband_mask])) # dB
    stopband_atten = minimum(abs.(mag_db[.!passband_mask])) # dB
    
    # Calculate delays
    group_delay = -diff(unwrap(phase)) ./ diff(freqs) .* (fs/(2π))
    phase_delay = -phase ./ (2π .* freqs)
    phase_delay[1] = phase_delay[2] # Handle division by zero at DC
    
    # Calculate transition width
    transition_width = if length(cutoff_freq) == 1
        # For lowpass/highpass, find width between -3dB and -20dB points
        db_20_idx = findfirst(x -> x <= -20, mag_db)
        abs(freqs[db_20_idx] - cutoff_freq[1])
    else
        # For bandpass/bandstop, use the width of both transition regions
        db_20_idx_low = findfirst(x -> x <= -20, mag_db)
        db_20_idx_high = findlast(x -> x <= -20, mag_db)
        (abs(freqs[db_20_idx_low] - cutoff_freq[1]) + 
         abs(freqs[db_20_idx_high] - cutoff_freq[end])) / 2
    end
    
    return (
        cutoff_freq = cutoff_freq,
        passband_ripple = passband_ripple,
        stopband_atten = stopband_atten,
        phase_delay = mean(phase_delay[passband_mask]),
        group_delay = mean(group_delay[passband_mask]),
        transition_width = transition_width,
        filter_type = filter_type
    )
end

"""
    print_filter_characteristics(filter; npoints::Int=1000)

Print a formatted summary of filter characteristics.

Arguments:
- `filter`: A digital filter object
- `npoints`: Number of frequency points for analysis (default: 1000)
"""
function print_filter_characteristics(filter; npoints::Int=1000)
    chars = get_filter_characteristics(filter; npoints=npoints)
    
    println("Filter Characteristics:")
    println("----------------------")
    println("Filter type: ", chars.filter_type)
    if length(chars.cutoff_freq) == 1
        println("Cutoff frequency: ", round(chars.cutoff_freq[1], digits=2), " Hz")
    else
        println("Cutoff frequencies: ", round.(chars.cutoff_freq, digits=2), " Hz")
    end
    println("Passband ripple: ", round(chars.passband_ripple, digits=2), " dB")
    println("Stopband attenuation: ", round(chars.stopband_atten, digits=2), " dB")
    println("Average phase delay: ", round(chars.phase_delay, digits=2), " samples")
    println("Average group delay: ", round(chars.group_delay, digits=2), " samples")
    println("Transition width: ", round(chars.transition_width, digits=2), " Hz")
end
