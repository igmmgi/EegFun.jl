# =============================================================================
# FILTER TYPES AND STRUCTURES
# =============================================================================

"""
    FilterInfo

A struct containing all information about a digital filter for EEG data processing.

# Fields
- `filter_type::String`: Filter type ("hp"=highpass, "lp"=lowpass)
- `filter_object`: The actual DSP.jl filter object (ZeroPoleGain for IIR, Vector{Float64} for FIR)
- `filter_method::String`: Filter implementation ("iir" or "fir")
- `cutoff_freq::Float64`: Cutoff frequency in Hz
- `sample_rate::Float64`: Sampling rate in Hz
- `order::Int`: Filter order (for IIR filters)
- `n_taps::Union{Nothing,Int}`: Number of taps (for FIR filters, nothing for IIR)
- `transition_band::Float64`: Width of transition band in Hz
"""
struct FilterInfo
    filter_type::String
    filter_object::Union{ZeroPoleGain,Vector{Float64}} # TODO: what is going on in DSP?
    filter_method::String
    cutoff_freq::Float64
    sample_rate::Float64
    order::Int
    n_taps::Union{Nothing,Int}
    transition_band::Float64
end

# =============================================================================
# FILTER CREATION
# =============================================================================
"""
    create_lowpass_filter(cutoff_freq::Real, sample_rate::Real; 
    filter_method::String = "iir", 
    order::Integer = 2, 
    transition_width::Real = 0.1)::FilterInfo

Create a digital filter object for the specified parameters.

# Arguments
- `cutoff_freq`: Cutoff frequency in Hz
- `sample_rate`: Sampling rate in Hz
- `filter_method`: String specifying filter implementation ("iir" or "fir")
- `order`: Filter order for IIR filters (default: 2, becomes effective order 4 with filtfilt)
- `transition_width`: Relative width of transition band as fraction of cutoff (default: 0.1 for EEG)

# Returns
- `FilterInfo`: A struct containing all filter information and the actual filter object

# Notes
- NOTE: When using filtfilt (default), effective filter order is approximately doubled

# Examples
```julia
# Create a lowpass filter and plot the response
filter_info = create_lowpass_filter(30.0, 256.0)
plot_filter_response(filter_info)
print_filter_characteristics(filter_info)
```
"""
function create_lowpass_filter(
    cutoff_freq::Real,
    sample_rate::Real;
    filter_method::String = "iir",
    order::Integer = 2,
    transition_width::Real = 0.1,
)
    return _create_filter("lp", cutoff_freq, sample_rate, filter_method = filter_method, order = order, transition_width = transition_width)
end


"""
    create_highpass_filter(cutoff_freq::Real, sample_rate::Real; 
    filter_method::String = "iir", 
    order::Integer = 2, 
    transition_width::Real = 0.1)::FilterInfo

Create a digital filter object for the specified parameters.

# Arguments
- `cutoff_freq`: Cutoff frequency in Hz
- `sample_rate`: Sampling rate in Hz
- `filter_method`: String specifying filter implementation ("iir" or "fir")
- `order`: Filter order for IIR filters (default: 2, becomes effective order 4 with filtfilt)
- `transition_width`: Relative width of transition band as fraction of cutoff (default: 0.1 for EEG)

# Returns
- `FilterInfo`: A struct containing all filter information and the actual filter object

# Notes
- NOTE: When using filtfilt (default), effective filter order is approximately doubled

# Examples
```julia
# Create a highpass filter and plot the response
filter_info = create_highpass_filter(0.1, 256.0)
plot_filter_response(filter_info)
print_filter_characteristics(filter_info)
```
"""
function create_highpass_filter(
    cutoff_freq::Real,
    sample_rate::Real;
    filter_method::String = "iir",
    order::Integer = 2,
    transition_width::Real = 0.1,
)
    return _create_filter("hp", cutoff_freq, sample_rate, filter_method = filter_method, order = order, transition_width = transition_width)
end


"""
    _create_filter(filter_type::String, cutoff_freq::Real, sample_rate::Real; filter_method::String, order::Integer = 2, transition_width::Real = 0.1)

Create a digital filter object for the specified parameters.

# Arguments
- `filter_type`: String specifying filter type ("hp"=highpass, "lp"=lowpass)
- `cutoff_freq`: Cutoff frequency in Hz
- `sample_rate`: Sampling rate in Hz
- `filter_method`: String specifying filter implementation ("iir" or "fir")
- `order`: Filter order for IIR filters (default: 2, becomes effective order 4 with filtfilt)
- `transition_width`: Relative width of transition band as fraction of cutoff (default: 0.1 for EEG)

# Returns
- `FilterInfo`: A struct containing all filter information and the actual filter object

# Notes
- NOTE: When using filtfilt (default), effective filter order is approximately doubled
"""
function _create_filter(
    filter_type::String,
    cutoff_freq::Real,
    sample_rate::Real;
    filter_method::String = "iir",
    order::Integer = 2,
    transition_width::Real = 0.1,
)

    # Input validation
    valid_filter_types = ("hp", "lp")
    if !(filter_type in valid_filter_types)
        @minimal_error "filter_type '$filter_type' must be one of: $valid_filter_types"
    end

    valid_filter_methods = ("iir", "fir")
    if !(filter_method in valid_filter_methods)
        @minimal_error "filter_method '$filter_method' must be one of: $valid_filter_methods"
    end

    if order <= 0
        @minimal_error "filter order must be positive: $order"
    end

    if sample_rate <= 0
        @minimal_error "sample_rate must be positive: $sample_rate"
    end

    if !(cutoff_freq > 0 && cutoff_freq < sample_rate / 2)
        @minimal_error "cutoff_freq must be between 0 and Nyquist frequency ($(sample_rate/2) Hz): $cutoff_freq"
    end

    # Create filter prototype/object based on type
    filter_prototypes = Dict("hp" => Highpass, "lp" => Lowpass)
    filter_prototype = filter_prototypes[filter_type](cutoff_freq)
    transition_band = cutoff_freq * transition_width
    n_taps = nothing

    # Create filter with chosen method
    if filter_method == "iir"
        filter_object = digitalfilter(filter_prototype, Butterworth(order); fs = sample_rate)
    elseif filter_method == "fir"
        n_taps = Int(ceil(4.0 * sample_rate / transition_band))
        n_taps = max(n_taps, 101)
        n_taps = nextpow(2, n_taps)
        n_taps += 1  # Always add 1 to make odd as nextpow2 is even
        filter_object = digitalfilter(filter_prototype, FIRWindow(hamming(n_taps)); fs = sample_rate)
    end

    # Create FilterInfo struct
    filter_info = FilterInfo(filter_type, filter_object, filter_method, cutoff_freq, sample_rate, order, n_taps, transition_band)

    return filter_info
end


# =============================================================================
# FILTER APPLICATION
# =============================================================================

"""
    _apply_filter!(dat::DataFrame, channels::Vector{Symbol}, filter; filter_func::String = "filtfilt")

Internal helper function to apply a digital filter to specified columns in a DataFrame.
Modifies the data in place.

# Arguments
- `dat::DataFrame`: DataFrame containing the data to filter
- `channels::Vector{Symbol}`: Vector of column names to filter
- `filter`: Digital filter object to apply (can be DSP.jl filter or FilterInfo)
- `filter_func::String`: Filtering function to use (default: "filtfilt" for two-pass filtering, or "filt" for one-pass filtering)
"""
function _apply_filter!(dat::DataFrame, channels::Vector{Symbol}, filter::FilterInfo; filter_func::String = "filtfilt")::Nothing
    filter_func = filter_func == "filtfilt" ? filtfilt : filt
    @inbounds for channel in channels
        @views dat[:, channel] .= filter_func(filter.filter_object, dat[:, channel])
    end
    return nothing
end

function _apply_filter!(dat::Vector{DataFrame}, channels::Vector{Symbol}, filter::FilterInfo; filter_func::String = "filtfilt")::Nothing
    _apply_filter!.(dat, Ref(channels), Ref(filter); filter_func = filter_func)
    return nothing
end

function _update_filter_info!(dat::EegData, filter_info::FilterInfo)::Nothing
    if filter_info.filter_type == "hp"
        dat.analysis_info.hp_filter = filter_info.cutoff_freq
    elseif filter_info.filter_type == "lp"
        dat.analysis_info.lp_filter = filter_info.cutoff_freq
    end
    return nothing
end

"""
    lowpass_filter!(dat::EegData, cutoff_freq::Real; kwargs...)

Apply a digital lowpass filter to EEG data. Modifies the data in place.

# Arguments
- `dat::EegData`: EegData object (ContinuousData, ErpData, or EpochData)
- `cutoff_freq::Real`: Cutoff frequency in Hz

# Keyword Arguments
- `order::Integer`: Filter order (default: 3)
- `transition_width::Real`: Relative width of transition band (default: 0.1)
- `filter_method::String`: Filter implementation ("iir" or "fir")
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `filter_func::String`: Filtering function ("filtfilt" or "filt")
"""
function lowpass_filter!(
    dat::EegData,
    cutoff_freq::Real;
    order::Int = 3,
    transition_width::Real = 0.1,
    filter_method::String = "iir",
    channel_selection::Function = channels(),
    filter_func::String = "filtfilt",
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for lowpass filtering"
        return nothing
    end
    @debug "lowpass_filter! applying lp filter at $cutoff_freq Hz to $(length(selected_channels)) channels"

    filter_info = create_lowpass_filter(
        cutoff_freq,
        dat.sample_rate;
        filter_method = filter_method,
        order = order,
        transition_width = transition_width,
    )
    _update_filter_info!(dat, filter_info)
    _apply_filter!(dat.data, selected_channels, filter_info, filter_func = filter_func)
end

"""
    highpass_filter!(dat::EegData, cutoff_freq::Real; kwargs...)

Apply a digital highpass filter to EEG data. Modifies the data in place.

# Arguments
- `dat::EegData`: EegData object (ContinuousData, ErpData, or EpochData)
- `cutoff_freq::Real`: Cutoff frequency in Hz

# Keyword Arguments
- `order::Integer`: Filter order (default: 1)
- `transition_width::Real`: Relative width of transition band (default: 0.25)
- `filter_method::String`: Filter implementation ("iir" or "fir")
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `filter_func::String`: Filtering function ("filtfilt" or "filt")
"""
function highpass_filter!(
    dat::EegData,
    cutoff_freq::Real;
    order::Int = 1,
    transition_width::Real = 0.25,
    filter_method::String = "iir",
    channel_selection::Function = channels(),
    filter_func::String = "filtfilt",
)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for highpass filtering"
        return nothing
    end
    @debug "highpass_filter! applying hp filter at $cutoff_freq Hz to $(length(selected_channels)) channels"

    filter_info = create_highpass_filter(
        cutoff_freq,
        dat.sample_rate;
        filter_method = filter_method,
        order = order,
        transition_width = transition_width,
    )
    _update_filter_info!(dat, filter_info)
    _apply_filter!(dat.data, selected_channels, filter_info, filter_func = filter_func)
end

# Generate non-mutating versions
@add_nonmutating lowpass_filter!
@add_nonmutating highpass_filter!

"""
    lowpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :lowpass)

Apply lowpass filter to EEG data based on configuration.
"""
function lowpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :lowpass)
    sec = getfield(filter_cfg, section)
    if sec.apply
        lowpass_filter!(dat, sec.freq; order = sec.order, filter_method = sec.method, filter_func = sec.func)
    end
end

"""
    highpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :highpass)

Apply highpass filter to EEG data based on configuration.
"""
function highpass_filter!(dat::EegData, filter_cfg::FilterConfig; section::Symbol = :highpass)
    sec = getfield(filter_cfg, section)
    if sec.apply
        highpass_filter!(dat, sec.freq; order = sec.order, filter_method = sec.method, filter_func = sec.func)
    end
end


# =============================================================================
# FILTER ANALYSIS AND CHARACTERIZATION
# =============================================================================

"""
    get_filter_characteristics(filter, sample_rate::Real, transition_width::Real; npoints::Int=1000)

Calculate and return key characteristics of a digital filter.

# Arguments
- `filter`: A digital filter object (FIR coefficients or DSP.jl filter)
- `sample_rate::Real`: Sampling rate in Hz
- `transition_width::Real`: Width of transition band in Hz
- `npoints::Int`: Number of frequency points for analysis (default: 1000)

# Returns
A NamedTuple containing:
- `cutoff_freq_3db`: Frequency at -3dB point(s) in Hz
- `cutoff_freq_6db`: Frequency at -6dB point(s) in Hz
- `passband_ripple`: Maximum ripple in passband in dB
- `stopband_atten`: Mean attenuation in stopband in dB
- `phase_delay`: Mean phase delay in passband (samples)
- `group_delay`: Mean group delay in passband (samples)
- `transition_band`: Width of transition band in Hz
- `filter_type`: Detected filter type ("hp" or "lp")
"""
function get_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000)
    # Calculate frequency response
    freqs = exp10.(range(log10(0.01), log10(filter_info.sample_rate / 2), length = npoints))  # log spacing
    w = 2π * freqs / filter_info.sample_rate

    # Get frequency response based on filter type
    if filter_info.filter_method == "fir"  # FIR filter coefficients
        n = 0:(length(filter_info.filter_object)-1)
        resp = [sum(filter_info.filter_object .* exp.(-im * w_k * n)) for w_k in w]
    else  # Other filter types
        resp = freqresp(filter_info.filter_object, w)
    end

    mag_db = 20 * log10.(abs.(resp))
    phase = angle.(resp)

    # Find -3dB and -6dB points using closest points method
    diffs_3db = abs.(mag_db .- (-3))
    diffs_6db = abs.(mag_db .- (-6))

    # Find local minima in the differences to detect crossings
    function find_crossings(diffs, threshold = 0.5)
        crossings = Int[]
        for i = 2:(length(diffs)-1)
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

    # Determine filter type and masks
    if filter_info.filter_type == "lp"
        passband_mask = freqs .<= cutoff_freq[1]
        stopband_mask = freqs .>= (cutoff_freq[1] + filter_info.transition_band)
    else  # highpass
        passband_mask = freqs .>= cutoff_freq[1]
        stopband_mask = freqs .<= (cutoff_freq[1] - filter_info.transition_band)
    end

    # Calculate passband ripple and stopband attenuation
    passband_ripple = maximum(abs.(mag_db[passband_mask])) # dB
    # Calculate mean stopband attenuation, replacing -Inf with -100 dB
    stopband_db = mag_db[stopband_mask]
    stopband_db = replace(stopband_db, -Inf => -100.0)
    stopband_atten = mean(stopband_db)

    # Calculate delays
    group_delay = -diff(DSP.unwrap(phase)) ./ diff(w)
    phase_delay = -phase ./ w
    phase_delay[1] = phase_delay[2] # Handle division by zero at DC

    return (
        cutoff_freq_3db = cutoff_freq_3db,
        cutoff_freq_6db = cutoff_freq_6db,
        passband_ripple = passband_ripple,
        stopband_atten = stopband_atten,
        phase_delay = mean(phase_delay[passband_mask]),
        group_delay = mean(group_delay[passband_mask[1:(end-1)]]),
        transition_band = filter_info.transition_band,
        filter_type = filter_info.filter_type,
    )
end

"""
    print_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000)

Print a formatted summary of filter characteristics.
"""
function print_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000)

    chars = get_filter_characteristics(filter_info; npoints = npoints)

    @info "--------------------------------------------------------------------------------"
    @info "Filter Characteristics: $(uppercase(chars.filter_type)) $(filter_info.cutoff_freq) Hz"
    @info "Cutoff: -3 dB = $(round.(chars.cutoff_freq_3db, digits=2)) Hz, -6 dB = $(round.(chars.cutoff_freq_6db, digits=2)) Hz"
    @info "Transition band: $(round(chars.transition_band, digits=1)) Hz"
    @info "Stopband attenuation: $(round(chars.stopband_atten, digits=1)) dB"
    @info "Single-pass characteristics: Group delay = $(round(chars.group_delay * 1000/filter_info.sample_rate, digits=1)) ms, Phase delay = $(round(chars.phase_delay * 1000/filter_info.sample_rate, digits=1)) ms "

end


# =============================================================================
# BATCH FILTERING API
# =============================================================================

"""Generate default output directory name for filter operation."""
function _default_filter_output_dir(input_dir::String, pattern::String, filter_type::String, freq::Real)
    joinpath(input_dir, "filtered_$(pattern)_$(filter_type)_$(freq)hz")
end

"""Process a single file through filtering pipeline."""
function _process_filter_file(filepath::String, output_path::String, filter_type::String, cutoff_freq::Real, condition_selection::Function)
    filename = basename(filepath)

    # Read data
    data = read_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Apply filter based on type
    if filter_type == "lp"
        lowpass_filter!.(data, cutoff_freq)
    else
        highpass_filter!.(data, cutoff_freq)
    end

    # Save
    jldsave(output_path; data = data)

    return BatchResult(true, filename, "Filtered successfully")
end

"""
    lowpass_filter(file_pattern::String, cutoff_freq::Real; kwargs...)

Batch lowpass filter EEG/ERP data from JLD2 files.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", etc.)
- `cutoff_freq::Real`: Cutoff frequency in Hz

# Keyword Arguments
- `input_dir::String`: Input directory (default: pwd())
- `participant_selection::Function`: Participant selection (default: participants())
- `condition_selection::Function`: Condition selection (default: conditions())
- `output_dir::Union{String, Nothing}`: Output directory
"""
function lowpass_filter(
    file_pattern::String,
    cutoff_freq::Real;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    return _run_filter_batch(
        file_pattern,
        cutoff_freq,
        "lp";
        input_dir = input_dir,
        participant_selection = participant_selection,
        condition_selection = condition_selection,
        output_dir = output_dir,
    )
end

"""
    highpass_filter(file_pattern::String, cutoff_freq::Real; kwargs...)

Batch highpass filter EEG/ERP data from JLD2 files.
"""
function highpass_filter(
    file_pattern::String,
    cutoff_freq::Real;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    return _run_filter_batch(
        file_pattern,
        cutoff_freq,
        "hp";
        input_dir = input_dir,
        participant_selection = participant_selection,
        condition_selection = condition_selection,
        output_dir = output_dir,
    )
end

# Internal common batch runner
function _run_filter_batch(
    file_pattern::String,
    cutoff_freq::Real,
    filter_type::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    log_file = "filter_$(filter_type).log"
    setup_global_logging(log_file)

    try
        @info "Batch $(filter_type == "lp" ? "lowpass" : "highpass") filtering started at $(now())"
        @info "  cutoff: $cutoff_freq Hz"

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        if cutoff_freq <= 0
            @minimal_error_throw("Cutoff frequency must be positive, got: $cutoff_freq")
        end

        # Setup directories
        output_dir = something(output_dir, _default_filter_output_dir(input_dir, file_pattern, filter_type, cutoff_freq))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        # Execute
        process_fn =
            (input_path, output_path) -> _process_filter_file(input_path, output_path, filter_type, cutoff_freq, condition_selection)
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Filtering")

        return _log_batch_summary(results, output_dir)
    finally
        _cleanup_logging(log_file, output_dir)
    end
end
