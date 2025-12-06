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
    return ErpData("grand_avg", cond_num, grand_avg_cond_name, grand_avg_data, first_erp.layout, first_erp.sample_rate, first_erp.analysis_info, total_epochs)
end

grand_average(erps::Vector{ErpData}, cond_num::Int) = _create_grand_average(erps, cond_num)

"""
Load and group ERP data by condition from multiple files.
Returns Dict{Int, Vector{ErpData}} mapping condition number to ERPs.
"""
function _load_and_group_erps(files::Vector{String}, input_dir::String, conditions::Vector{Int})

    files = sort(files, by=natural_sort_key)
    all_erps_by_condition = OrderedDict{Int,Vector{ErpData}}()
    for (i, file) in enumerate(files)
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"

        # Load ERP data (using load_data which finds by type)
        erps_data = load_data(input_path)
        if isnothing(erps_data)
            @minimal_warning "No data variables found in $file. Skipping."
            continue
        end
        
        # Validate that data is Vector{ErpData}
        if !(erps_data isa Vector{<:ErpData})
            @minimal_warning "Invalid data type in $file: expected Vector{ErpData}, got $(typeof(erps_data)). Skipping."
            continue
        end

        # Group ERPs by condition number
        for erp in erps_data
            cond_num = erp.condition
            push!(get!(all_erps_by_condition, cond_num, ErpData[]), erp)
        end
    end

    conditions === nothing && return all_erps_by_condition

    # Filter by condition numbers if specified
    return OrderedDict(cond_num => all_erps_by_condition[cond_num] 
                          for cond_num in conditions 
                          if haskey(all_erps_by_condition, cond_num))

end

_load_and_group_erps(files::Vector{String}, input_dir::String, conditions::Int) = _load_and_group_erps(files, input_dir, [conditions])

"""
Create grand averages for all conditions in the grouped data.
Returns Vector{ErpData}.
"""
function _create_all_grand_averages(erps_by_condition::Dict{Int,Vector{ErpData}})
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
                 participants::Union{Int, Vector{Int}, Nothing} = nothing,
                 conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                 output_dir::Union{String, Nothing} = nothing)

Batch process ERP data files to create grand averages across participants.

This function loads JLD2 files containing ERP data from multiple participants, groups by condition,
and creates grand averages by averaging the EEG channel data across participants for each condition.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "erps_original")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Grand average all ERP files in current directory
grand_average("erps_cleaned")

# Process specific participants and conditions
grand_average("erps_cleaned", 
            input_dir = "/path/to/data", 
            participants = [1, 2, 3], 
            conditions = [1, 2])

# Specify custom output directory
grand_average("erps_cleaned", 
            input_dir = "/path/to/data", 
            output_dir = "/path/to/output")
```
"""
function grand_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    conditions::Union{Int,Vector{Int},Nothing} = nothing,
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
        files = _find_batch_files(file_pattern, input_dir, participants)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Load and group ERP data by condition
        erps_by_condition = _load_and_group_erps(files, input_dir, conditions)

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
