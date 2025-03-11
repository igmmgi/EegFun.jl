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
function filter_data!(
    dat::DataFrame,
    columns,
    filter_type::String,
    filter_method::String,
    sample_rate::Real;
    order::Integer = 3,
    transition_width::Real = 0.25,
    print_filter_characteristics::Bool = true,
    plot_filter_response::Bool = false,
)
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
        transition_band = transition_width * minimum(freq_low, freq_high)  
    else
        if !(isa(freq, Real) && freq > 0 && freq < sample_rate/2)
            throw(ArgumentError("frequency must be between 0 and Nyquist frequency ($(sample_rate/2) Hz)"))
        end
        transition_band = transition_width * freq  
    end

    # Create filter prototype based on type
    if filter_type == "hp"
        filter_prototype = Highpass(freq+(transition_band/2))
    elseif filter_type == "lp"
        filter_prototype = Lowpass(freq+(transition_band/2))
    elseif filter_type == "bp"
        freq_low, freq_high = freq
        filter_prototype = Bandpass(freq_low, freq_high)
    elseif filter_type == "bs"
        freq_low, freq_high = freq
        filter_prototype = Bandstop(freq_low, freq_high)
    end

    # Create filter with chosen method
    if filter_method == "iir"
        filter = digitalfilter(filter_prototype, Butterworth(order); fs = sample_rate)
    elseif filter_method == "fir"
        n_taps = Int(ceil(3.3 * sample_rate / transition_band))
        if n_taps % 2 == 0  
            n_taps += 1     
        end
        filter = digitalfilter(filter_prototype, FIRWindow(hamming(n_taps)); fs = sample_rate)
    end

    if print_filter_characteristics
        print_filter_characteristics(filter, sample_rate, freq, transition_band)
    end
    if plot_filter_response
        plot_filter_response(filter, sample_rate, freq, transition_band)
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
function filter_data(
    dat::DataFrame,
    columns,
    filter_type::String,
    filter_method::String,
    sample_rate::Real;
    order::Integer = 3,
    transition_width::Real = 0.25,
    print_filter_characteristics::Bool = true,
    plot_filter_response::Bool = false,
)
    dat_out = deepcopy(dat)
    filter_data!(dat_out, columns, filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
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

function filter_data(dat::Union{ContinuousData,ErpData}, filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
    dat_out = deepcopy(dat)
    filter_data!(dat_out.data, dat_out.layout.label, filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
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
function filter_data(dat::EpochData, filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
    dat_out = deepcopy(dat)
    for epoch in eachindex(dat_out.data)
        filter_data!(dat_out.data[epoch], filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
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
function filter_data(dat::EpochData, filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
    dat_out = deepcopy(dat)
    for epoch in eachindex(dat_out.data)
        filter_data!(dat_out.data[epoch], filter_type, filter_method, sample_rate; order, transition_width, print_filter_characteristics, plot_filter_response)
    end
    return dat_out
end



############################################################
"""
    get_filter_characteristics(filter, sample_rate::Real, f_p::Union{Real,Tuple}, transition_width::Real; npoints::Int=1000)

Calculate and return key characteristics of a digital filter.

Arguments:
- `filter`: A digital filter object
- `sample_rate`: Sampling rate in Hz
- `f_p`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `transition_width`: Width of transition band in Hz
- `npoints`: Number of frequency points for analysis (default: 1000)

Returns:
A NamedTuple containing:
- `cutoff_freq_3db`: Frequency at -3dB point(s) in Hz
- `cutoff_freq_6db`: Frequency at -6dB point(s) in Hz
- `passband_ripple`: Maximum ripple in passband in dB
- `stopband_atten`: Minimum attenuation in stopband in dB
- `phase_delay`: Average phase delay in samples
- `group_delay`: Average group delay in samples
- `transition_width`: Width of transition band in Hz
- `filter_type`: Detected filter type ("hp", "lp", "bp", or "bs")
"""
function get_filter_characteristics(filter, sample_rate::Real, f_p::Union{Real,Tuple}, transition_width::Real; npoints::Int=1000)
    # Calculate frequency response
    freqs = range(0, sample_rate/2, length=npoints)
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
    mid_response = mag_db[floor(Int, npoints/2)]
    end_response = mag_db[end]
    
    # Determine filter type and masks
    if f_p isa Tuple  # If given two frequencies, must be bandpass or bandstop
        if length(cutoff_freq_3db) >= 2  # Found two -3dB points
            if mag_db[floor(Int, (crossings_3db[1] + crossings_3db[end])/2)] > -3
                filter_type = "bp"  # Bandpass: passes frequencies between cutoffs
                passband_mask = (freqs .>= cutoff_freq[1]) .& (freqs .<= cutoff_freq[end])
                stopband_mask = (freqs .<= (cutoff_freq[1] - transition_width)) .| 
                               (freqs .>= (cutoff_freq[end] + transition_width))
            else
                filter_type = "bs"  # Bandstop: blocks frequencies between cutoffs
                passband_mask = (freqs .<= cutoff_freq[1]) .| (freqs .>= cutoff_freq[end])
                stopband_mask = (freqs .>= (cutoff_freq[1] + transition_width)) .& 
                               (freqs .<= (cutoff_freq[end] - transition_width))
            end
        else
            filter_type = "bp"  # Default to bandpass if can't determine clearly
            passband_mask = (freqs .>= cutoff_freq[1]) .& (freqs .<= cutoff_freq[end])
            stopband_mask = (freqs .<= (cutoff_freq[1] - transition_width)) .| 
                           (freqs .>= (cutoff_freq[end] + transition_width))
        end
    else  # Single frequency - must be lowpass or highpass
        if start_response > -3 && end_response < -20
            filter_type = "lp"  # Lowpass
            passband_mask = freqs .<= cutoff_freq[1]
            stopband_mask = freqs .>= (cutoff_freq[1] + transition_width)
        else  # highpass
            filter_type = "hp"
            passband_mask = freqs .>= cutoff_freq[1]
            stopband_mask = freqs .<= (cutoff_freq[1] - transition_width)
        end
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
    print_filter_characteristics(filter, sample_rate::Real, f_p::Union{Real,Tuple}, transition_width::Real; npoints::Int=1000)

Print a formatted summary of filter characteristics.

Arguments:
- `filter`: A digital filter object
- `sample_rate`: Sampling rate in Hz
- `f_p`: Cutoff frequency (or tuple of frequencies for bandpass/bandstop)
- `transition_width`: Width of transition band in Hz
- `npoints`: Number of frequency points for analysis (default: 1000)
"""
function filter_characteristics(filter, sample_rate::Real, f_p::Union{Real,Tuple}, transition_width::Real; npoints::Int=1000)
    chars = get_filter_characteristics(filter, sample_rate, f_p, transition_width; npoints=npoints)
    println("Filter Characteristics:")
    println("-------------------------")
    
    println("Type: ", uppercase(chars.filter_type))
    if !isempty(chars.cutoff_freq_3db)
        println("Cutoff (-3 dB): ", round.(chars.cutoff_freq_3db, digits=1), " Hz")
    end
    if !isempty(chars.cutoff_freq_6db)
        println("Cutoff (-6 dB): ", round.(chars.cutoff_freq_6db, digits=1), " Hz")
    end
    
    println("Transition width: ", round(chars.transition_width, digits=1), " Hz")
    println("Stopband attenuation: ", round(chars.stopband_atten, digits=1), " dB")
    
    println("\nSingle-pass characteristics:")
    println("Group delay: ", round(chars.group_delay, digits=1), " samples (", 
            round(chars.group_delay * 1000/sample_rate, digits=1), " ms)")
    println("Phase delay: ", round(chars.phase_delay, digits=1), " samples (",
            round(chars.phase_delay * 1000/sample_rate, digits=1), " ms)")
            
    println("\nNote: Using filtfilt() results in:")
    println("- Zero phase delay")
    println("- Zero group delay")
    println("- Double the stopband attenuation")
end

"""
    plot_filter_response(filter, fs::Real, f_p::Union{Real,Tuple}, transition_band::Real; ylimit::Tuple=(-100, 5), xlimit::Tuple=(1, fs/2))

Plot the frequency response of a digital filter with ideal response overlay.
"""
function plot_filter_response(filter, fs::Real, f_p::Union{Real,Tuple}, transition_band::Real; ylimit::Tuple=(-100, 5), xlimit::Tuple=(1, fs/2))
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
    
    # Calculate frequency response with more points near the transition band
    n_points = 2000
    if f_p isa Tuple  # Bandstop/Bandpass
        f_low, f_high = f_p
        freqs = [0.0; 
                exp10.(range(log10(0.1), log10(f_low/2), length=n_points÷6));
                range(f_low/2, f_low + 2*transition_band, length=n_points÷6);
                range(f_low + 2*transition_band, f_high - 2*transition_band, length=n_points÷3);  # More points in between
                range(f_high - 2*transition_band, f_high + 2*transition_band, length=n_points÷6);
                exp10.(range(log10(f_high + 2*transition_band), log10(fs/2), length=n_points÷6))]
    else  # Lowpass/Highpass
        freqs = [0.0; 
                exp10.(range(log10(0.1), log10(f_p/2), length=n_points÷4));
                range(f_p/2, f_p + 2*transition_band, length=n_points÷2);
                exp10.(range(log10(f_p + 2*transition_band), log10(fs/2), length=n_points÷4))]
    end
    w = 2π * freqs / fs
    
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
    mid_response = mag_db[floor(Int, n_points/2)]
    end_response = mag_db[end]
    
    # Determine filter type
    if start_response > -3 && end_response < -20
        filter_type = "lp"
    elseif start_response < -20 && end_response > -3
        filter_type = "hp"
    elseif start_response < -20 && end_response < -20 && mid_response > -3
        filter_type = "bp"
    elseif start_response > -3 && end_response > -3 && mid_response < -20
        filter_type = "bs"
    else
        filter_type = "lp"
    end
    
    # Calculate actual stopband attenuation
    if f_p isa Tuple  # Bandpass/Bandstop
        f_low, f_high = f_p  # Extract frequencies from tuple
        if filter_type == "bp"
            stopband_mask = (freqs .<= (f_low - transition_band)) .| (freqs .>= (f_high + transition_band))
        else  # bandstop
            stopband_mask = (freqs .>= (f_low + transition_band)) .& (freqs .<= (f_high - transition_band))
        end
    else  # Lowpass/Highpass
        if filter_type == "lp"
            stopband_mask = freqs .>= (f_p + transition_band)
        else  # highpass
            stopband_mask = freqs .<= (f_p - transition_band)
        end
    end
    
    stopband_db = mag_db[stopband_mask]
    stopband_db = replace(stopband_db, -Inf => -100.0)
    actual_stopband_atten = mean(stopband_db)
    actual_stopband_linear = 10^(actual_stopband_atten/20)
    
    # Plot ideal response based on filter type
    ideal_response = zeros(length(freqs))
    for (i, f) in enumerate(freqs)
        if filter_type == "lp" && !(f_p isa Tuple)
            if f <= f_p
                ideal_response[i] = 1.0
            elseif f >= f_p + transition_band
                ideal_response[i] = actual_stopband_linear
            else
                ideal_response[i] = 1.0 * (f_p + transition_band - f) / transition_band + 
                                  actual_stopband_linear * (f - f_p) / transition_band
            end
        elseif filter_type == "hp" && !(f_p isa Tuple)
            if f >= f_p
                ideal_response[i] = 1.0
            elseif f <= f_p - transition_band
                ideal_response[i] = actual_stopband_linear
            else
                ideal_response[i] = 1.0 * (f - (f_p - transition_band)) / transition_band + 
                                  actual_stopband_linear * (f_p - f) / transition_band
            end
        elseif filter_type == "bp" && f_p isa Tuple
            f_low, f_high = f_p  # Extract frequencies here too
            if f >= f_low && f <= f_high
                ideal_response[i] = 1.0
            elseif f <= f_low - transition_band || f >= f_high + transition_band
                ideal_response[i] = actual_stopband_linear
            elseif f > f_low - transition_band && f < f_low  # lower transition
                ideal_response[i] = 1.0 * (f - (f_low - transition_band)) / transition_band + 
                                  actual_stopband_linear * (f_low - f) / transition_band
            elseif f > f_high && f < f_high + transition_band  # upper transition
                ideal_response[i] = 1.0 * (f_high + transition_band - f) / transition_band + 
                                  actual_stopband_linear * (f - f_high) / transition_band
            end
        elseif filter_type == "bs" && f_p isa Tuple
            f_low, f_high = f_p  # Extract frequencies
            if f <= f_low - transition_band || f >= f_high + transition_band
                ideal_response[i] = 1.0  # passband
            elseif f >= f_low && f <= f_high
                ideal_response[i] = actual_stopband_linear  # stopband
            elseif f > f_low - transition_band && f < f_low  # lower transition into stopband
                ideal_response[i] = actual_stopband_linear * (f - (f_low - transition_band)) / transition_band + 
                                  1.0 * (f_low - f) / transition_band
            elseif f > f_high && f < f_high + transition_band  # upper transition out of stopband
                ideal_response[i] = 1.0 * (f - f_high) / transition_band + 
                                  actual_stopband_linear * (f_high + transition_band - f) / transition_band
            end
        end
    end
    
    # Add transition regions
    if f_p isa Tuple
        f_s1 = f_p[1] - transition_band
        f_s2 = f_p[2] + transition_band
        transition_regions = [(f_s1, f_p[1]), (f_p[2], f_s2)]
    else
        f_s = f_p + transition_band
        transition_regions = [(f_p, f_s)]
    end
    
    for (start_f, end_f) in transition_regions
        vspan!(ax, start_f, end_f, color=(:gray, 0.2))
        vlines!(ax, [start_f, end_f], color=:gray, linestyle=:dash)
    end
    
    # Plot responses
    #lines!(ax, freqs, mag_db, label="Actual", color=:black, linewidth=4)
    ideal_mag_db = 20 * log10.(ideal_response)
    lines!(ax, freqs, ideal_mag_db, color=:green, linestyle=:dash, label="Ideal",linewidth=3)
    
    # Add reference lines
    hlines!(ax, [-3, -6], color=:gray, linestyle=:dash)
    text!(ax, fs/2, -3, text="-3 dB", align=(:right, :center), fontsize=22)
    text!(ax, fs/2, -6, text="-6 dB", align=(:right, :center), fontsize=22)
    
    # X-axis ticks
    if f_p isa Tuple
        xticks = [0, 1, 2, 5, 10, 20, 30, f_p[1], f_p[2], f_s1, f_s2, 100, fs/2]
    else
        xticks = [0, 1, 2, 5, 10, 20, 30, f_p, f_s, 100, fs/2]
    end
    ax.xticks = (xticks, string.(round.(xticks, digits=0)))
    
    axislegend(;position=:lb, labelsize=32)
    display(fig)
    
    return fig, ax
end


# Test print_filter_characteristics and plot_filter_response
# Set up filter parameters
sample_rate = 2048 # Hz

# low-pass
cutoff_freq = 50  # Hz  # Where we want the -3dB point
transition_band = 0.25 * cutoff_freq  # = 10 Hz

n_taps = Int(ceil(3.3 * sample_rate / transition_band))
lp_fir_filter = digitalfilter(Lowpass(cutoff_freq+(transition_band/2)), FIRWindow(hamming(n_taps)), fs=sample_rate)
print_filter_characteristics(lp_fir_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(lp_fir_filter, sample_rate, cutoff_freq, transition_band);

lp_irr_filter = digitalfilter(Lowpass(cutoff_freq+(transition_band/2)), Butterworth(6), fs=sample_rate)
print_filter_characteristics(lp_irr_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(lp_irr_filter, sample_rate, cutoff_freq, transition_band);

# high-pass
cutoff_freq = 1  # Hz  # Where we want the -3dB point
transition_band = 0.25 * cutoff_freq  # = 10 Hz

n_taps = Int(ceil(3.3 * sample_rate / transition_band))
if n_taps % 2 == 0  # If even
    n_taps += 1     # Make it odd
end
hp_fir_filter = digitalfilter(Highpass(cutoff_freq+(transition_band/2)), FIRWindow(hamming(n_taps)), fs=sample_rate)
print_filter_characteristics(hp_fir_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(hp_fir_filter, sample_rate, cutoff_freq, transition_band);

hp_irr_filter = digitalfilter(Highpass(cutoff_freq+(transition_band/2)), Butterworth(2), fs=sample_rate)
print_filter_characteristics(hp_irr_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(hp_irr_filter, sample_rate, cutoff_freq, transition_band);

# band-stop
cutoff_freq = (45, 55)  # Hz  # Where we want the -3dB point
transition_band = 0.25 * minimum(cutoff_freq)  # = 10 Hz

n_taps = Int(ceil(3.3 * sample_rate / transition_band))
if n_taps % 2 == 0  # If even
    n_taps += 1     # Make it odd
end
bs_fir_filter = digitalfilter(Bandstop(cutoff_freq[1], cutoff_freq[2]), FIRWindow(hamming(n_taps)), fs=sample_rate)
print_filter_characteristics(bs_fir_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(bs_fir_filter, sample_rate, cutoff_freq, transition_band);

bs_irr_filter = digitalfilter(Bandstop(cutoff_freq[1], cutoff_freq[2]), Butterworth(2), fs=sample_rate)
print_filter_characteristics(bs_irr_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(bs_irr_filter, sample_rate, cutoff_freq, transition_band);


# band-pass
cutoff_freq = (1, 30)  # Hz  # Where we want the -3dB point
transition_band = 0.25 * minimum(cutoff_freq)  # = 10 Hz

n_taps = Int(ceil(3.3 * sample_rate / transition_band))
if n_taps % 2 == 0  # If even
    n_taps += 1     # Make it odd
end 
bp_fir_filter = digitalfilter(Bandpass(cutoff_freq[1], cutoff_freq[2]), FIRWindow(hamming(n_taps)), fs=sample_rate)
print_filter_characteristics(bp_fir_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(bp_fir_filter, sample_rate, cutoff_freq, transition_band);

bp_irr_filter = digitalfilter(Bandpass(cutoff_freq[1], cutoff_freq[2]), Butterworth(6), fs=sample_rate)
print_filter_characteristics(bp_irr_filter, sample_rate, cutoff_freq, transition_band);
plot_filter_response(bp_irr_filter, sample_rate, cutoff_freq, transition_band);






