"""
Batch grand averaging of ERP data across participants.
"""

#=============================================================================
    GRANDAVERAGE-SPECIFIC VALIDATION
=============================================================================#

"""Validate that file pattern is for ERP data."""
function _validate_erps_pattern_grand_average(pattern::String)
    !contains(pattern, "erps") &&
        return "grand_average only works with ERP data. File pattern must contain 'erps', got: '$pattern'"
    return nothing
end

"""Generate default output directory name for grand averaging."""
function _default_grand_average_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "grand_average_$(pattern)")
end

#=============================================================================
    GRANDAVERAGE-SPECIFIC HELPERS
=============================================================================#

"""
Create a grand average by averaging ERP data across participants for a specific condition.
"""
function _create_grand_average(erps::Vector{ErpData}, cond_num::Int)
    if isempty(erps)
        @minimal_error_throw("Cannot create grand average from empty ERP list")
    end

    # Validate that all ERPs have the same structure
    have_same_structure(erps) || @minimal_error_throw("ERPs have inconsistent structure")

    first_erp = erps[1]

    # Get metadata columns and EEG channels
    metadata_cols = meta_labels(first_erp)
    eeg_channels = setdiff(propertynames(first_erp.data), metadata_cols)

    # Create a copy of the first ERP's data as the base
    grand_avg_data = copy(first_erp.data)

    # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
    cols_to_remove = [:condition, :condition_name, :n_epochs]
    for col in cols_to_remove
        if hasproperty(grand_avg_data, col)
            select!(grand_avg_data, Not(col))
        end
    end

    # Average EEG channels across participants
    for ch in eeg_channels
        # Collect data from all participants for this channel
        # Stack as columns: n_timepoints x n_participants
        channel_matrix = hcat([erp.data[!, ch] for erp in erps]...)

        # Average across participants (mean of each time point)
        grand_avg_data[!, ch] = vec(mean(channel_matrix, dims = 2))
    end

    # Calculate total number of epochs across all participants
    total_epochs = sum(erp.n_epochs for erp in erps)
    grand_avg_cond_name = "grand_avg_$(first_erp.condition_name)"

    # Use the layout and analysis info from the first ERP
    return ErpData(
        "grand_avg",
        cond_num,
        grand_avg_cond_name,
        grand_avg_data,
        first_erp.layout,
        first_erp.sample_rate,
        first_erp.analysis_info,
        total_epochs,
    )
end

grand_average(erps::Vector{ErpData}, cond_num::Int) = _create_grand_average(erps, cond_num)

"""
Load and group ERP data by condition from multiple files.
Returns OrderedDict{Int, Vector{ErpData}} mapping condition number to ERPs.
"""
function _load_and_group_erps(files::Vector{String}, input_dir::String, condition_selection::Function)
    # Load all ERPs and group by condition
    all_erps = load_all_data(ErpData, files, input_dir)
    erps_by_condition = group_by_condition(all_erps)

    # Apply condition selection to the sorted condition numbers
    all_cond_nums = collect(keys(erps_by_condition))  # Already sorted
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    # Return only the selected conditions
    return OrderedDict(num => erps_by_condition[num] for num in selected_cond_nums)
end

"""
Create grand averages for all conditions in the grouped data.
Returns Vector{ErpData}.
"""
function _create_all_grand_averages(erps_by_condition::AbstractDict{Int,Vector{ErpData}})
    grand_averages = ErpData[]

    # Sort conditions to ensure consistent ordering
    for cond_num in sort(collect(keys(erps_by_condition)))
        erps = erps_by_condition[cond_num]
        @info "Creating grand average for condition $cond_num (n=$(length(erps)) participants)"

        if length(erps) < 2
            @minimal_warning "Only $(length(erps)) participant(s) for condition $cond_num. Skipping grand average."
            continue
        end

        # Create grand average
        grand_avg = _create_grand_average(erps, cond_num)
        push!(grand_averages, grand_avg)

        @info "  ✓ Created grand average with $(nrow(grand_avg.data)) time points"
    end

    return grand_averages
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    grand_average(file_pattern::String; 
                 input_dir::String = pwd(), 
                 participant_selection::Function = participants(),
                 condition_selection::Function = conditions(),
                 output_dir::Union{String, Nothing} = nothing)

Batch process ERP data files to create grand averages across participants.

This function loads JLD2 files containing ERP data from multiple participants, groups by condition,
and creates grand averages by averaging the EEG channel data across participants for each condition.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "erps_original")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Grand average all ERP files in current directory
grand_average("erps_cleaned")

# Process specific participants and conditions
grand_average("erps_cleaned", 
            input_dir = "/path/to/data", 
            participant_selection = participants([1, 2, 3]), 
            condition_selection = conditions([1, 2]))

# Custom predicate (participants > 5)
grand_average("erps_cleaned", participant_selection = x -> x .> 5)
```
"""
function grand_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "grand_average.log"
    setup_global_logging(log_file)

    try
        @info "Batch grand averaging started at $(now())"
        @log_call "grand_average"

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_erps_pattern_grand_average(file_pattern)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_grand_average_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Load and group ERP data by condition
        erps_by_condition = _load_and_group_erps(files, input_dir, condition_selection)

        if isempty(erps_by_condition)
            @minimal_warning "No valid ERP data found in any files"
            return nothing
        end

        @info "Found conditions: $(sort(collect(keys(erps_by_condition))))"

        # Create grand averages for each condition
        grand_averages = _create_all_grand_averages(erps_by_condition)

        if isempty(grand_averages)
            @minimal_warning "No grand averages created (insufficient participants for any condition)"
            return nothing
        end

        # Save grand averages
        output_file = "grand_average_$(file_pattern).jld2"
        output_path = joinpath(output_dir, output_file)
        jldsave(output_path; data = grand_averages)

        @info "Grand averaging complete! Created $(length(grand_averages)) grand averages"
        @info "Output saved to: $output_path"

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

# ============================================================================
# Grand averaging for TimeFreqData
# ============================================================================

"""
Create a grand average by averaging TimeFreq data across participants for a specific condition.
"""
function _create_tf_grand_average(tf_data::Vector{TimeFreqData}, cond_num::Int)
    if isempty(tf_data)
        @minimal_error_throw("Cannot create grand average from empty TimeFreqData list")
    end

    first_tf = tf_data[1]
    
    # Validate all have same structure (same freqs, times, channels)
    first_times = unique(first_tf.data.time)
    first_freqs = unique(first_tf.data.freq)
    first_channels = channel_labels(first_tf)
    
    for tf in tf_data[2:end]
        if unique(tf.data.time) != first_times
            @minimal_error_throw("TimeFreqData objects have inconsistent time vectors")
        end
        if unique(tf.data.freq) != first_freqs
            @minimal_error_throw("TimeFreqData objects have inconsistent frequency vectors")
        end
        if channel_labels(tf) != first_channels
            @minimal_error_throw("TimeFreqData objects have inconsistent channels")
        end
    end

    # Create a copy of the first TF's data as the base
    grand_avg_data = copy(first_tf.data)

    # Average each channel across participants
    for ch in first_channels
        # Collect data from all participants for this channel
        # Stack as columns: n_rows x n_participants
        channel_matrix = hcat([tf.data[!, ch] for tf in tf_data]...)

        # Average across participants
        grand_avg_data[!, ch] = vec(mean(channel_matrix, dims=2))
    end

    grand_avg_cond_name = "grand_avg_$(first_tf.condition_name)"

    return TimeFreqData(
        "grand_avg",
        cond_num,
        grand_avg_cond_name,
        grand_avg_data,
        first_tf.layout,
        first_tf.sample_rate,
        first_tf.method,
        first_tf.analysis_info
    )
end

grand_average(tf_data::Vector{TimeFreqData}, cond_num::Int) = _create_tf_grand_average(tf_data, cond_num)

"""
Load and group TimeFreqData by condition from multiple files.
Returns OrderedDict{Int, Vector{TimeFreqData}} mapping condition number to TF data.
"""
function _load_and_group_tf_data(files::Vector{String}, input_dir::String, condition_selection::Function)
    # Load all TF data and group by condition
    all_tf = load_all_data(TimeFreqData, files, input_dir)
    tf_by_condition = group_by_condition(all_tf)

    # Apply condition selection to the sorted condition numbers
    all_cond_nums = collect(keys(tf_by_condition))
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    # Return only the selected conditions
    return OrderedDict(num => tf_by_condition[num] for num in selected_cond_nums)
end

"""
Create grand averages for all conditions in the grouped TF data.
Returns Vector{TimeFreqData}.
"""
function _create_all_tf_grand_averages(tf_by_condition::AbstractDict{Int,Vector{TimeFreqData}})
    grand_averages = TimeFreqData[]

    # Sort conditions to ensure consistent ordering
    for cond_num in sort(collect(keys(tf_by_condition)))
        tf_data = tf_by_condition[cond_num]
        @info "Creating TF grand average for condition $cond_num (n=$(length(tf_data)) participants)"

        if length(tf_data) < 2
            @minimal_warning "Only $(length(tf_data)) participant(s) for condition $cond_num. Skipping grand average."
            continue
        end

        # Create grand average
        grand_avg = _create_tf_grand_average(tf_data, cond_num)
        push!(grand_averages, grand_avg)

        n_times = length(unique(grand_avg.data.time))
        n_freqs = length(unique(grand_avg.data.freq))
        @info "  ✓ Created TF grand average: $(n_times) times × $(n_freqs) freqs"
    end

    return grand_averages
end

"""
    tf_grand_average(file_pattern::String; 
                     input_dir::String = pwd(), 
                     participant_selection::Function = participants(),
                     condition_selection::Function = conditions(),
                     output_dir::Union{String, Nothing} = nothing)

Batch process TimeFreqData files to create grand averages across participants.

This function loads JLD2 files containing TimeFreqData from multiple participants, groups by condition,
and creates grand averages by averaging the power data across participants for each condition.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "tf_epochs_wavelet")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Grand average all TF files in current directory
tf_grand_average("tf_epochs_wavelet")

# Process specific participants and conditions
tf_grand_average("tf_epochs_wavelet", 
                 input_dir = "/path/to/data", 
                 participant_selection = participants([1, 2, 3]), 
                 condition_selection = conditions([1, 2]))
```
"""
function tf_grand_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    # Setup logging
    log_file = "tf_grand_average.log"
    setup_global_logging(log_file)

    try
        @info "Batch TF grand averaging started at $(now())"
        @log_call "tf_grand_average"

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, joinpath(input_dir, "tf_grand_average_$(file_pattern)"))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Load and group TF data by condition
        tf_by_condition = _load_and_group_tf_data(files, input_dir, condition_selection)

        if isempty(tf_by_condition)
            @minimal_warning "No valid TimeFreqData found in any files"
            return nothing
        end

        @info "Found conditions: $(sort(collect(keys(tf_by_condition))))"

        # Create grand averages for each condition
        grand_averages = _create_all_tf_grand_averages(tf_by_condition)

        if isempty(grand_averages)
            @minimal_warning "No TF grand averages created (insufficient participants for any condition)"
            return nothing
        end

        # Save grand averages
        output_file = "tf_grand_average_$(file_pattern).jld2"
        output_path = joinpath(output_dir, output_file)
        jldsave(output_path; data=grand_averages)

        @info "TF grand averaging complete! Created $(length(grand_averages)) grand averages"
        @info "Output saved to: $output_path"

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

# Add group_by_condition for TimeFreqData
function group_by_condition(tf_data::Vector{TimeFreqData})
    by_condition = OrderedDict{Int, Vector{TimeFreqData}}()
    
    for tf in tf_data
        cond_num = tf.condition
        if !haskey(by_condition, cond_num)
            by_condition[cond_num] = TimeFreqData[]
        end
        push!(by_condition[cond_num], tf)
    end
    
    # Sort by condition number
    return OrderedDict(k => by_condition[k] for k in sort(collect(keys(by_condition))))
end
