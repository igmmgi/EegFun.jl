"""
Automatic epoch rejection based on statistical criteria.

This module provides functions to automatically reject epochs based on
statistical measures (variance, max, min, absolute max, range, kurtosis)
using z-score thresholds. This is useful for removing epochs with artifacts
without manual inspection.
"""

using Statistics
using StatsBase

#=============================================================================
    REJECTION INFORMATION STRUCTURE
=============================================================================#

"""
    EpochRejectionInfo

Stores information about which epochs were rejected and why.

# Fields
- `n_original::Int`: Number of epochs before rejection
- `n_remaining::Int`: Number of epochs after rejection
- `rejected_epochs::Vector{Int}`: Indices of rejected epochs
- `rejected_by_variance::Vector{Int}`: Epochs rejected due to high variance
- `rejected_by_max::Vector{Int}`: Epochs rejected due to high maximum values
- `rejected_by_min::Vector{Int}`: Epochs rejected due to low minimum values
- `rejected_by_abs::Vector{Int}`: Epochs rejected due to high absolute values
- `rejected_by_range::Vector{Int}`: Epochs rejected due to large range
- `rejected_by_kurtosis::Vector{Int}`: Epochs rejected due to high kurtosis
- `z_criterion::Float64`: Z-score criterion used for rejection
"""
struct EpochRejectionInfo
    n_original::Int
    n_remaining::Int
    rejected_epochs::Vector{Int}
    rejected_by_variance::Vector{Int}
    rejected_by_max::Vector{Int}
    rejected_by_min::Vector{Int}
    rejected_by_abs::Vector{Int}
    rejected_by_range::Vector{Int}
    rejected_by_kurtosis::Vector{Int}
    z_criterion::Float64
end

#=============================================================================
    CORE REJECTION FUNCTIONS
=============================================================================#

"""
    reject_epochs_automatic!(dat::EpochData, z_criterion::Real;
                             channel_selection::Function = channels())::EpochRejectionInfo

Automatically reject epochs based on statistical criteria using z-score threshold.

This function calculates six statistical measures for each epoch (variance, max, min,
absolute max, range, kurtosis) across selected channels. For each measure, the maximum
across channels is taken for each epoch, then z-scored across epochs. Epochs exceeding
the z-criterion for any measure are rejected.

# Arguments
- `dat::EpochData`: Epoched EEG data to process
- `z_criterion::Real`: Z-score threshold for rejection (typically 2.0 or 3.0)
- `channel_selection::Function`: Channel predicate for selecting channels to analyze (default: all channels)

# Returns
- `EpochRejectionInfo`: Information about which epochs were rejected and why

# Effects
- Modifies the input data in-place by removing rejected epochs

# Examples
```julia
using eegfun, JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Reject epochs with z-score > 2.0
rejection_info = reject_epochs_automatic!(epochs, 2.0)

# Check results
println("Original epochs: \$(rejection_info.n_original)")
println("Remaining epochs: \$(rejection_info.n_remaining)")
println("Rejected epochs: \$(rejection_info.rejected_epochs)")

# Save cleaned data
save("participant_1_epochs_cleaned.jld2", "epochs", epochs)
```

# Notes
- Common z-criteria: 2.0 (more aggressive), 2.5, 3.0 (more conservative)
- Rejection is based on ANY metric exceeding the criterion
- Uses maximum across channels to identify global artifacts
- All six metrics are calculated independently and combined with OR logic
"""
function reject_epochs_automatic!(
    dat::EpochData,
    z_criterion::Real;
    channel_selection::Function = channels(),
)::EpochRejectionInfo

    @info "Starting automatic epoch rejection with z-criterion = $z_criterion"
    
    # Validate inputs
    _validate_rejection_inputs(dat, z_criterion)
    
    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    
    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for epoch rejection")
    end
    
    @info "Using $(length(selected_channels)) channels for rejection analysis"
    
    # Calculate metrics for all epochs
    n_epochs = length(dat.data)
    metrics = _calculate_epoch_metrics(dat, selected_channels)
    
    # Identify rejected epochs for each metric
    rejection_results = _identify_rejected_epochs(metrics, z_criterion)
    
    # Combine all rejected epochs
    rejected_epochs = sort(unique(vcat(
        rejection_results[:variance],
        rejection_results[:max],
        rejection_results[:min],
        rejection_results[:abs],
        rejection_results[:range],
        rejection_results[:kurtosis]
    )))
    
    n_rejected = length(rejected_epochs)
    n_remaining = n_epochs - n_rejected
    
    @info "Rejection analysis complete:"
    @info "  Original epochs: $n_epochs"
    @info "  Rejected epochs: $n_rejected"
    @info "  Remaining epochs: $n_remaining"
    @info "  Rejection breakdown:"
    @info "    Variance: $(length(rejection_results[:variance])) epochs"
    @info "    Maximum: $(length(rejection_results[:max])) epochs"
    @info "    Minimum: $(length(rejection_results[:min])) epochs"
    @info "    Absolute: $(length(rejection_results[:abs])) epochs"
    @info "    Range: $(length(rejection_results[:range])) epochs"
    @info "    Kurtosis: $(length(rejection_results[:kurtosis])) epochs"
    
    # Create rejection info before modifying data
    rejection_info = EpochRejectionInfo(
        n_epochs,
        n_remaining,
        rejected_epochs,
        rejection_results[:variance],
        rejection_results[:max],
        rejection_results[:min],
        rejection_results[:abs],
        rejection_results[:range],
        rejection_results[:kurtosis],
        Float64(z_criterion)
    )
    
    # Remove rejected epochs
    if !isempty(rejected_epochs)
        epochs_to_keep = setdiff(1:n_epochs, rejected_epochs)
        dat.data = dat.data[epochs_to_keep]
        @info "Removed $n_rejected epochs. $(length(dat.data)) epochs remaining."
    else
        @info "No epochs rejected."
    end
    
    return rejection_info
end


"""
    reject_epochs_automatic(dat::EpochData, z_criterion::Real;
                           channel_selection::Function = channels())::Tuple{EpochData, EpochRejectionInfo}

Non-mutating version of reject_epochs_automatic!. Returns new EpochData with rejected epochs removed.

# Returns
- `Tuple{EpochData, EpochRejectionInfo}`: Clean data and rejection information

# Examples
```julia
# Non-mutating version preserves original data
clean_epochs, rejection_info = reject_epochs_automatic(epochs, 2.0)
```
"""
function reject_epochs_automatic(
    dat::EpochData,
    z_criterion::Real;
    channel_selection::Function = channels(),
)::Tuple{EpochData, EpochRejectionInfo}

    # Create a deep copy
    dat_copy = EpochData(
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info)
    )
    
    # Apply rejection
    rejection_info = reject_epochs_automatic!(dat_copy, z_criterion; channel_selection = channel_selection)
    
    return (dat_copy, rejection_info)
end


#=============================================================================
    INTERNAL HELPER FUNCTIONS
=============================================================================#

"""
Validate inputs for epoch rejection.
"""
function _validate_rejection_inputs(dat::EpochData, z_criterion::Real)
    if isempty(dat.data)
        @minimal_error_throw("Cannot reject epochs from empty EpochData")
    end
    
    if z_criterion <= 0
        @minimal_error_throw("Z-criterion must be positive, got $z_criterion")
    end
    
    if length(dat.data) < 3
        @minimal_warning "Only $(length(dat.data)) epochs available. Statistical rejection may not be meaningful with so few epochs."
    end
end


"""
Calculate statistical metrics for all epochs.

Returns a Dict with keys: :variance, :max, :min, :abs, :range, :kurtosis
Each value is a vector of length n_epochs containing the maximum of that metric across channels.
"""
function _calculate_epoch_metrics(dat::EpochData, selected_channels::Vector{Symbol})::Dict{Symbol, Vector{Float64}}
    n_epochs = length(dat.data)
    
    # Pre-allocate arrays for metrics
    variance_max = zeros(Float64, n_epochs)
    max_vals = zeros(Float64, n_epochs)
    min_vals = zeros(Float64, n_epochs)
    abs_max = zeros(Float64, n_epochs)
    range_vals = zeros(Float64, n_epochs)
    kurtosis_max = zeros(Float64, n_epochs)
    
    # Calculate metrics for each epoch
    for (epoch_idx, epoch) in enumerate(dat.data)
        # For each channel, calculate metric, then take max across channels
        channel_variance = Float64[]
        channel_max = Float64[]
        channel_min = Float64[]
        channel_abs = Float64[]
        channel_range = Float64[]
        channel_kurtosis = Float64[]
        
        for ch in selected_channels
            channel_data = epoch[!, ch]
            
            push!(channel_variance, var(channel_data))
            push!(channel_max, maximum(channel_data))
            push!(channel_min, minimum(channel_data))
            push!(channel_abs, maximum(abs.(channel_data)))
            push!(channel_range, maximum(channel_data) - minimum(channel_data))
            push!(channel_kurtosis, kurtosis(channel_data))
        end
        
        # Take maximum across channels for each metric
        variance_max[epoch_idx] = maximum(channel_variance)
        max_vals[epoch_idx] = maximum(channel_max)
        min_vals[epoch_idx] = minimum(channel_min)  # Note: minimum of minimums!
        abs_max[epoch_idx] = maximum(channel_abs)
        range_vals[epoch_idx] = maximum(channel_range)
        kurtosis_max[epoch_idx] = maximum(channel_kurtosis)
    end
    
    return Dict{Symbol, Vector{Float64}}(
        :variance => variance_max,
        :max => max_vals,
        :min => min_vals,
        :abs => abs_max,
        :range => range_vals,
        :kurtosis => kurtosis_max
    )
end


"""
Identify epochs that exceed the z-criterion for each metric.

Returns a Dict with keys: :variance, :max, :min, :abs, :range, :kurtosis
Each value is a vector of epoch indices that were rejected by that criterion.
"""
function _identify_rejected_epochs(metrics::Dict{Symbol, Vector{Float64}}, z_criterion::Real)::Dict{Symbol, Vector{Int}}
    results = Dict{Symbol, Vector{Int}}()
    
    for (metric_name, values) in metrics
        # Calculate z-scores
        z_scores = zscore(values)
        
        # Find epochs exceeding criterion
        rejected = findall(z_scores .> z_criterion)
        
        results[metric_name] = rejected
    end
    
    return results
end


#=============================================================================
    REPORTING FUNCTIONS
=============================================================================#

"""
    Base.show(io::IO, info::EpochRejectionInfo)

Display rejection information in a human-readable format.
"""
function Base.show(io::IO, info::EpochRejectionInfo)
    println(io, "EpochRejectionInfo:")
    println(io, "  Z-criterion: $(info.z_criterion)")
    println(io, "  Original epochs: $(info.n_original)")
    println(io, "  Remaining epochs: $(info.n_remaining)")
    println(io, "  Rejected epochs: $(length(info.rejected_epochs))")
    println(io, "")
    println(io, "  Rejection breakdown:")
    println(io, "    Variance:  $(length(info.rejected_by_variance)) epochs")
    println(io, "    Maximum:   $(length(info.rejected_by_max)) epochs")
    println(io, "    Minimum:   $(length(info.rejected_by_min)) epochs")
    println(io, "    Absolute:  $(length(info.rejected_by_abs)) epochs")
    println(io, "    Range:     $(length(info.rejected_by_range)) epochs")
    println(io, "    Kurtosis:  $(length(info.rejected_by_kurtosis)) epochs")
    
    if length(info.rejected_epochs) <= 20
        println(io, "")
        println(io, "  Rejected epoch indices: $(info.rejected_epochs)")
    end
end


"""
    save_rejection_report(info::EpochRejectionInfo, filename::String)

Save rejection information to a text file.

# Arguments
- `info::EpochRejectionInfo`: Rejection information to save
- `filename::String`: Output filename (typically ends in .txt)

# Examples
```julia
save_rejection_report(rejection_info, "participant_1_rejection_report.txt")
```
"""
function save_rejection_report(info::EpochRejectionInfo, filename::String)
    open(filename, "w") do io
        println(io, repeat("=", 70))
        println(io, "AUTOMATIC EPOCH REJECTION REPORT")
        println(io, repeat("=", 70))
        println(io, "")
        println(io, "Z-criterion: $(info.z_criterion)")
        println(io, "Original epochs: $(info.n_original)")
        println(io, "Remaining epochs: $(info.n_remaining)")
        println(io, "Rejected epochs: $(length(info.rejected_epochs)) ($(round(100 * length(info.rejected_epochs) / info.n_original, digits=1))%)")
        println(io, "")
        println(io, "REJECTION BREAKDOWN:")
        println(io, repeat("-", 70))
        println(io, "Variance:  $(length(info.rejected_by_variance)) epochs - $(info.rejected_by_variance)")
        println(io, "Maximum:   $(length(info.rejected_by_max)) epochs - $(info.rejected_by_max)")
        println(io, "Minimum:   $(length(info.rejected_by_min)) epochs - $(info.rejected_by_min)")
        println(io, "Absolute:  $(length(info.rejected_by_abs)) epochs - $(info.rejected_by_abs)")
        println(io, "Range:     $(length(info.rejected_by_range)) epochs - $(info.rejected_by_range)")
        println(io, "Kurtosis:  $(length(info.rejected_by_kurtosis)) epochs - $(info.rejected_by_kurtosis)")
        println(io, "")
        println(io, "ALL REJECTED EPOCHS:")
        println(io, repeat("-", 70))
        println(io, info.rejected_epochs)
        println(io, "")
        println(io, repeat("=", 70))
    end
    
    @info "Rejection report saved to: $filename"
end


#=============================================================================
    BATCH PROCESSING FUNCTIONS
=============================================================================#

"""Generate default output directory name for rejection operation."""
function _default_rejection_output_dir(input_dir::String, pattern::String, z_criterion::Real)
    z_str = replace(string(z_criterion), "." => "p")  # Replace . with p for filename
    joinpath(input_dir, "rejected_z$(z_str)_$(pattern)")
end


"""
Process a single epoch file through automatic rejection pipeline.
Returns BatchResult with success/failure info.
"""
function _process_rejection_file(
    filepath::String,
    output_path::String,
    z_criterion::Real,
    channel_selection::Function,
)
    filename = basename(filepath)
    
    # Load data
    file_data = load(filepath)
    
    # Try common variable names for epoched data
    epoch_var_names = ["epochs", "epoch_data", "data"]
    epochs_data = nothing
    
    for var_name in epoch_var_names
        if haskey(file_data, var_name)
            epochs_data = file_data[var_name]
            break
        end
    end
    
    if isnothing(epochs_data)
        return BatchResult(false, filename, "No epoched data variable found (tried: $(epoch_var_names))")
    end
    
    if !(epochs_data isa EpochData)
        return BatchResult(false, filename, "Data is not EpochData type")
    end
    
    # Reject epochs
    try
        rejection_info = reject_epochs_automatic!(epochs_data, z_criterion; channel_selection = channel_selection)
        
        # Save cleaned epochs
        save(output_path, "epochs", epochs_data)
        
        # Save rejection report in same directory as output file
        output_dir = dirname(output_path)
        base_filename = splitext(basename(filename))[1]
        report_path = joinpath(output_dir, "$(base_filename)_rejection_report.txt")
        save_rejection_report(rejection_info, report_path)
        
        message = "Rejected $(length(rejection_info.rejected_epochs)) of $(rejection_info.n_original) epochs"
        return BatchResult(true, filename, message)
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end


"""
    reject_epochs_automatic(file_pattern::String, z_criterion::Real;
                           input_dir::String = pwd(),
                           channel_selection::Function = channels(),
                           participants::Union{Int, Vector{Int}, Nothing} = nothing,
                           output_dir::Union{String, Nothing} = nothing)

Batch process epoch files with automatic rejection and save to a new directory.

This function processes multiple participant files at once, automatically rejecting
epochs based on statistical criteria.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "epochs", "epochs_cleaned")
- `z_criterion::Real`: Z-score threshold for rejection (typically 2.0 or 3.0)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)

# Examples
```julia
# Reject epochs with z > 2.0 for all participants
reject_epochs_automatic("epochs", 2.0)

# Specific participants only
reject_epochs_automatic("epochs", 2.5, participants = [1, 2, 3])

# Custom output directory
reject_epochs_automatic("epochs", 3.0, output_dir = "/path/to/output")

# Only use specific channels for rejection analysis
reject_epochs_automatic("epochs", 2.0,
                       channel_selection = channels(x -> startswith.(string.(x), "F")))

# Full example
reject_epochs_automatic("epochs", 2.5,
                       input_dir = "/data/study1",
                       participants = 1:20)
```

# Output
- Creates new directory with cleaned epoch data files
- Each output file contains "epochs" variable with cleaned EpochData
- Rejection reports saved as *_rejection_report.txt files
- Log file saved to output directory

# Notes
- Common z-criteria: 2.0 (aggressive), 2.5 (moderate), 3.0 (conservative)
- Higher z-criterion = fewer epochs rejected = more lenient
- Lower z-criterion = more epochs rejected = more strict
"""
function reject_epochs_automatic(
    file_pattern::String,
    z_criterion::Real;
    input_dir::String = pwd(),
    channel_selection::Function = channels(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)
    
    # Setup logging
    log_file = "reject_epochs_automatic.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch automatic epoch rejection started at $(now())"
        @log_call "reject_epochs_automatic" (file_pattern, z_criterion)
        
        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end
        
        if z_criterion <= 0
            @minimal_error_throw("Z-criterion must be positive, got $z_criterion")
        end
        
        # Setup directories
        output_dir = something(output_dir, _default_rejection_output_dir(input_dir, file_pattern, z_criterion))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Z-criterion: $z_criterion"
        @info "Channel selection: $(channel_selection == channels() ? "all channels" : "custom")"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) ->
            _process_rejection_file(input_path, output_path, z_criterion, channel_selection)
        
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Rejecting epochs")
        
        _log_batch_summary(results, output_dir)
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end

