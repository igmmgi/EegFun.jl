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

# Notes
- This struct encapsulates all filter parameters and metadata
- Makes it easy to pass filter information to print and plot functions
- Provides a consistent interface regardless of filter type
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
    create_filter(filter_type::String, filter_method::String, filter_freq::Real, sample_rate::Real; 
                 order::Integer = 2, transition_width::Real = 0.1)

Create a digital filter object for the specified parameters.

# Arguments
- `filter_type`: String specifying filter type ("hp"=highpass, "lp"=lowpass)
- `filter_method`: String specifying filter implementation ("iir" or "fir")
- `cutoff_freq`: Cutoff frequency in Hz
- `sample_rate`: Sampling rate in Hz
- `order`: Filter order for IIR filters (default: 2, becomes effective order 4 with filtfilt)
- `transition_width`: Relative width of transition band as fraction of cutoff (default: 0.1 for EEG)

# Returns
- `FilterInfo`: A struct containing all filter information and the actual filter object

# Notes
- Validates all input parameters
- Creates appropriate filter prototype based on type and method
- Returns a FilterInfo object that can be reused for multiple datasets
- Optimized for EEG analysis with appropriate defaults
- NOTE: When using filtfilt (default), effective filter order is approximately doubled

# Examples
```julia
# Create a lowpass filter
filter_info = create_filter("lp", "iir", 30.0, 256.0)

# Create a highpass filter with custom order
filter_info = create_filter("hp", "iir", 1.0, 256.0, order = 4)

# Create a FIR filter
filter_info = create_filter("lp", "fir", 40.0, 256.0)

# Plot the filter response
plot_filter_response(filter_info)

# Print filter characteristics
print_filter_characteristics(filter_info)
```
"""
function create_filter(
    filter_type::String,
    filter_method::String,
    cutoff_freq::Real,
    sample_rate::Real;
    order::Integer = 2,
    transition_width::Real = 0.1,
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
        # For EEG, use higher order Butterworth for steeper rolloff
        filter_object = digitalfilter(filter_prototype, Butterworth(order); fs = sample_rate)
    elseif filter_method == "fir"
        # Calculate number of taps (ensure not too small + next pow2 + odd) for FIR filter
        n_taps = Int(ceil(4.0 * sample_rate / transition_band))
        n_taps = max(n_taps, 101)
        n_taps = nextpow(2, n_taps)
        n_taps += 1  # Always add 1 to make odd as nextpow2 is even
        filter_object = digitalfilter(filter_prototype, FIRWindow(hamming(n_taps)); fs = sample_rate)
    end

    # Create FilterInfo struct
    filter_info = FilterInfo(
        filter_type,
        filter_object,
        filter_method,
        Float64(cutoff_freq),
        Float64(sample_rate),
        order,
        n_taps,
        Float64(transition_band),
    )

    return filter_info
end


# =============================================================================
# FILTER APPLICATION
# =============================================================================

"""
    _apply_filter!(dat::DataFrame, channels::Vector{Symbol}, filter; filter_func::Function = filtfilt)

Internal helper function to apply a digital filter to specified columns in a DataFrame.
Modifies the data in place.

# Arguments
- `dat::DataFrame`: DataFrame containing the data to filter
- `channels::Vector{Symbol}`: Vector of column names to filter
- `filter`: Digital filter object to apply (can be DSP.jl filter or FilterInfo)
- `filter_func::Function`: Filtering function to use (default: filtfilt, for two-pass filtering). Use `filt` for one-pass filtering.
"""
function _apply_filter!(
    dat::DataFrame,
    channels::Vector{Symbol},
    filter::FilterInfo;
    filter_func::String = "filtfilt",
)::Nothing
    filter_func = filter_func == "filtfilt" ? filtfilt : filt
    @inbounds for channel in channels
        @views dat[:, channel] .= filter_func(filter.filter_object, dat[:, channel])
    end
    return nothing
end

function _apply_filter!(
    dat::Vector{DataFrame},
    channels::Vector{Symbol},
    filter::FilterInfo;
    filter_func::String = "filtfilt",
)::Nothing
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
    filter_data!(dat::EegData, filter_type::String, cutoff_freq::Real; kwargs...)

Apply a digital filter to EEG data. Modifies the data in place.

# Arguments
- `dat::EegData`: EegData object (ContinuousData, ErpData, or EpochData)
- `filter_type::String`: Filter type ("hp"=highpass, "lp"=lowpass)
- `cutoff_freq::Real`: Cutoff frequency in Hz

# Keyword Arguments
- `order::Integer`: Filter order for IIR filters (default: 3)
- `transition_width::Real`: Relative width of transition band as fraction of cutoff (default: 0.25)
- `filter_method::String`: Filter implementation ("iir" or "fir")
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)
- `filter_func::String`: Filtering function to use (default: filtfilt, for two-pass filtering). Use `filt` for one-pass filtering.

# Examples
```julia
# Apply highpass filter
filter_data!(dat, "hp", 1.0)

# Apply lowpass filter with custom parameters
filter_data!(dat, "lp", 30.0, order = 4, filter_method = "fir")

# Apply filter to specific channels
filter_data!(dat, "hp", 1.0, channel_selection = channels([:Fp1, :Fp2]))
```
"""
function filter_data!(
    dat::EegData,
    filter_type::String,
    cutoff_freq::Real;
    order::Int = filter_type == "hp" ? 1 : 3,
    transition_width::Real = filter_type == "hp" ? 0.25 : 0.1,
    filter_method::String = "iir",
    channel_selection::Function = channels(),
    filter_func::String = "filtfilt",
)

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for filtering"
        return nothing
    end
    @debug "filter_data! applying $filter_type filter to $(length(selected_channels)) channels"

    # Create and apply filter
    filter_info = create_filter(
        filter_type,
        filter_method,
        cutoff_freq,
        dat.sample_rate;
        order = order,
        transition_width = transition_width,
    )
    _update_filter_info!(dat, filter_info)
    _apply_filter!(dat.data, selected_channels, filter_info, filter_func = filter_func)

end

# generates all non-mutating versions
@add_nonmutating filter_data!

"""
    filter_data!(dat::EegData, filter_cfg::Dict)

Apply filters to EEG data based on configuration dictionary. Modifies the data in place.

# Arguments
- `dat::EegData`: EegData object to filter
- `filter_cfg::Dict`: Configuration dictionary containing filter settings

# Configuration Structure
The `filter_cfg` should contain filter sections with keys like:
- `"highpass"`, `"lowpass"`: Standard filters
- `"ica_highpass"`, `"ica_lowpass"`: ICA-specific filters
- Any other filter section with keys `"apply"`, `"freq"`, `"order"`, `"method"`, `"func"`

# Examples
```julia
# Apply standard filters
filter_data!(dat, cfg["filter"])

# Apply ICA-specific filters
filter_data!(dat, cfg["filter"], filter_sections = ["ica_highpass", "ica_lowpass"])

# Apply only highpass filter
filter_data!(dat, Dict("highpass" => Dict("apply" => true, "freq" => 1.0, "order" => 2, "method" => "iir", "func" => "filtfilt")))
```
"""
# Wrapper method that works directly with FilterSection
function filter_data!(dat::EegData, section_filter::FilterSection)
    if section_filter.apply
        filter_data!(
            dat,
            section_filter.type,
            section_filter.freq;
            order = section_filter.order,
            filter_method = section_filter.method,
            filter_func = section_filter.func,
        )
    end
end

# Main method for FilterConfig
function filter_data!(dat::EegData, filter_cfg::FilterConfig; filter_sections::Vector{Symbol} = [:highpass, :lowpass])
    for section in filter_sections
        section_filter = getfield(filter_cfg, section)
        filter_data!(dat, section_filter)
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
- `transition_width`: Width of transition band in Hz
- `filter_type`: Detected filter type ("hp" or "lp")
"""
function get_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000)
    # Calculate transition_width from transition_band and cutoff_freq
    transition_width = filter_info.transition_band / filter_info.cutoff_freq

    # Calculate frequency response
    freqs = exp10.(range(log10(0.01), log10(filter_info.sample_rate/2), length = npoints))  # log spacing
    # freqs = range(0, sample_rate/2, length=npoints) # linear spacing
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

    # Determine filter type based on response shape
    start_response = mag_db[2]  # Use second point to avoid DC issues
    end_response = mag_db[end]

    # Determine filter type and masks
    if filter_info.filter_type == "lp"
        passband_mask = freqs .<= cutoff_freq[1]
        stopband_mask = freqs .>= (cutoff_freq[1] + transition_width)
    else  # highpass
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
        transition_width = transition_width,
        filter_type = filter_info.filter_type,
    )
end

"""
    print_filter_characteristics(filter, sample_rate::Real, filter_freq::Union{Real,Tuple}, transition_width::Real; npoints::Int=1000)

Print a formatted summary of filter characteristics.

# Arguments
- `filter`: A digital filter object
- `sample_rate::Real`: Sampling rate in Hz
- `filter_freq::Union{Real,Tuple}`: Cutoff frequency in Hz
- `transition_width::Real`: Width of transition band in Hz
- `npoints::Int`: Number of frequency points for analysis (default: 1000)

# Examples
```julia
# Print characteristics of a lowpass filter
filter = create_filter("lp", "fir", 40.0, 1000.0)
print_filter_characteristics(filter, 1000.0, 40.0, 10.0)
```
"""
function print_filter_characteristics(filter_info::FilterInfo; npoints::Int = 1000)

    chars = get_filter_characteristics(filter_info; npoints = npoints)

    @info "--------------------------------------------------------------------------------"
    @info "Filter Characteristics: $(uppercase(chars.filter_type)) $(filter_info.cutoff_freq) Hz"
    @info "Cutoff: -3 dB = $(round.(chars.cutoff_freq_3db, digits=2)) Hz, -6 dB = $(round.(chars.cutoff_freq_6db, digits=2)) Hz"
    @info "Transition width: $(round(chars.transition_width, digits=1)) Hz"
    @info "Stopband attenuation: $(round(chars.stopband_atten, digits=1)) dB"
    @info "Single-pass characteristics: Group delay = $(round(chars.group_delay * 1000/filter_info.sample_rate, digits=1)) ms, Phase delay = $(round(chars.phase_delay * 1000/filter_info.sample_rate, digits=1)) ms "

end


"""
Batch filtering for EEG/ERP data.
"""

#=============================================================================
    FILTER-SPECIFIC VALIDATION
=============================================================================#

"""Validate filter-specific parameters, returning error message or nothing."""
function _validate_filter_params(cutoff_freq::Real, filter_type::String)
    cutoff_freq <= 0 && return "Cutoff frequency must be positive, got: $cutoff_freq"
    filter_type ∉ ["lp", "hp"] && return "Filter type must be one of: lp, hp, got: $filter_type"
    return nothing
end

"""Generate default output directory name for filter operation."""
function _default_filter_output_dir(input_dir::String, pattern::String, filter_type::String, freq::Real)
    joinpath(input_dir, "filtered_$(pattern)_$(filter_type)_$(freq)hz")
end

#=============================================================================
    FILTER-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through filtering pipeline.
Returns BatchResult with success/failure info.
"""
function _process_filter_file(
    filepath::String,
    output_path::String,
    filter_type::String,
    cutoff_freq::Real,
    condition_selection::Function,
)
    filename = basename(filepath)

    # Load data
    data = load_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Apply filter (mutates data in-place)
    filter_data!.(data, filter_type, cutoff_freq)

    # Save (always use "data" as variable name since load_data finds by type)
    jldsave(output_path; data = data)

    return BatchResult(true, filename, "Filtered successfully")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    filter(file_pattern::String, cutoff_freq::Real; 
           input_dir::String = pwd(), 
           filter_type::String = "lp", 
           participant_selection::Function = participants(),
           condition_selection::Function = conditions(),
           output_dir::Union{String, Nothing} = nothing)

Filter EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `cutoff_freq::Real`: Cutoff frequency in Hz
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `filter_type::String`: Type of filter ("lp", "hp")
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on filter settings)

# Example
```julia
# Filter all epochs at 30 Hz
filter("epochs", 30.0)

# Filter specific participant
filter("epochs", 30.0, participants=3)

# Filter specific participants and conditions
filter("epochs", 30.0, participants=[3, 4], conditions=[1, 2])
```
"""
function filter(
    file_pattern::String,
    cutoff_freq::Real;
    input_dir::String = pwd(),
    filter_type::String = "lp",
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "filter.log"
    setup_global_logging(log_file)

    try
        @info ""
        @info "Batch filtering started at $(now())"
        @log_call "filter"

        # Validation 
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        if (error_msg = _validate_filter_params(cutoff_freq, filter_type)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir =
            something(output_dir, _default_filter_output_dir(input_dir, file_pattern, filter_type, cutoff_freq))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Create processing function with captured parameters
        @info "Filter settings: filter_type=\"$filter_type\", cutoff: $cutoff_freq Hz"
        process_fn =
            (input_path, output_path) ->
                _process_filter_file(input_path, output_path, filter_type, cutoff_freq, condition_selection)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Filtering")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
