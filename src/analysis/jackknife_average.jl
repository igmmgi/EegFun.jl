"""
Jackknife averaging for ERP/LRP data.

The jackknife technique creates leave-one-out averages: for each participant,
compute the average of all other participants (excluding that participant).
This is commonly used with LRP data to reduce variance for statistical testing.
"""

#=============================================================================
    JACKKNIFE-SPECIFIC VALIDATION
=============================================================================#

"""Validate jackknife parameters."""
function _validate_jackknife_params(erps::Vector{ErpData})
    length(erps) < 2 && return "Need at least 2 participants for jackknife averaging"

    # Validate that all ERPs have the same structure
    if !have_same_structure(erps)
        return "ERPs have inconsistent structure (sample rate, number of samples, or channel labels)"
    end

    return nothing
end

"""Generate default output directory name for jackknife operation."""
function _default_jackknife_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "jackknife_$(pattern)")
end

#=============================================================================
    JACKKNIFE-SPECIFIC HELPERS
=============================================================================#

"""
Create jackknife averages: for each participant, average all other participants.

Returns a vector where element i is the average of all participants except participant i.
"""
function _create_jackknife_averages(erps::Vector{ErpData})::Vector{ErpData}
    n_participants = length(erps)

    if n_participants < 2
        @minimal_error_throw("Need at least 2 participants for jackknife averaging, got $(n_participants)")
    end

    # Get metadata columns and EEG channels from first ERP
    first_erp = erps[1]
    metadata_cols = meta_labels(first_erp)
    eeg_channels = setdiff(propertynames(first_erp.data), metadata_cols)

    # Pre-allocate result vector
    jackknife_erps = ErpData[]

    @info "Creating jackknife averages for $n_participants participants"

    for excluded_idx = 1:n_participants
        # Get indices of all participants except the current one
        included_indices = setdiff(1:n_participants, excluded_idx)
        included_erps = erps[included_indices]

        @info "  Participant $excluded_idx: averaging $(length(included_erps)) other participants"

        # Create a copy of the first included ERP's data as the base
        jackknife_data = copy(included_erps[1].data)

        # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
        cols_to_remove = [:condition, :condition_name, :n_epochs]
        for col in cols_to_remove
            if hasproperty(jackknife_data, col)
                select!(jackknife_data, Not(col))
            end
        end

        # Average EEG channels across included participants
        for ch in eeg_channels
            # Collect data from all included participants for this channel
            # Stack as columns: n_timepoints x n_included_participants
            channel_matrix = hcat([erp.data[!, ch] for erp in included_erps]...)

            # Average across participants (mean of each time point)
            jackknife_data[!, ch] = vec(mean(channel_matrix, dims = 2))
        end

        # Calculate total number of epochs across included participants
        total_epochs = sum(erp.n_epochs for erp in included_erps)

        # Get condition info from first ERP and update for jackknife
        cond_name = included_erps[1].condition_name
        condition = included_erps[1].condition
        jackknife_cond_name = "$(cond_name)_jackknife_$(excluded_idx)"

        # Create ErpData object for this jackknife average
        jackknife_erp = ErpData(
            first_erp.file,
            condition,
            jackknife_cond_name,
            jackknife_data,
            first_erp.layout,
            first_erp.sample_rate,
            copy(first_erp.analysis_info),
            total_epochs,
        )

        push!(jackknife_erps, jackknife_erp)
    end

    @info "Created $(length(jackknife_erps)) jackknife averages"
    return jackknife_erps
end

"""
Load ERP/LRP data from multiple files and organize by condition.
Returns Dict{Int, Vector{ErpData}} mapping condition number to ERPs from all participants.
"""
function _load_and_group_for_jackknife(
    files::Vector{String},
    input_dir::String,
    condition_selection::Function,
    data_var::String,
)
    # data_var parameter kept for backwards compatibility but not used - load_data() finds by type
    all_erps_by_condition = Dict{Int,Vector{ErpData}}()
    participant_ids = Int[]

    for (i, file) in enumerate(files)
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"

        # Extract participant ID from filename (assumes format like "1_pattern.jld2")
        m = match(r"^(\d+)_", file)
        participant_id = m !== nothing ? parse(Int, m.captures[1]) : i
        push!(participant_ids, participant_id)

        # Load data (using load_data which finds by type)
        data = load_data(input_path)

        if isnothing(data)
            @minimal_warning "No data variables found in $file. Skipping."
            continue
        end

        # Validate that data is Vector{ErpData} or ErpData
        if !(data isa Union{Vector{<:ErpData},ErpData})
            @minimal_warning "Invalid data type in $file: expected Vector{ErpData} or ErpData, got $(typeof(data)). Skipping."
            continue
        end

        # Handle both single ErpData and Vector{ErpData}
        if data isa ErpData
            data = [data]
        end

        # Select conditions if specified
        data = _condition_select(data, condition_selection)

        # Group by condition
        for erp in data
            # For ErpData, condition is stored in the struct
            cond_num = erp.condition
            if !haskey(all_erps_by_condition, cond_num)
                all_erps_by_condition[cond_num] = ErpData[]
            end
            push!(all_erps_by_condition[cond_num], erp)
        end
    end

    return all_erps_by_condition, participant_ids
end

#=============================================================================
    MAIN API FUNCTIONS
=============================================================================#

"""
    jackknife_average(erps::Vector{ErpData})::Vector{ErpData}

Create jackknife averages from a vector of ERP/LRP data.

For each participant i, creates an average of all other participants (excluding i).
This leave-one-out approach is commonly used with LRP data to reduce variance
for statistical testing.

# Arguments
- `erps::Vector{ErpData}`: Vector of ERP/LRP data, one per participant

# Returns
- `Vector{ErpData}`: Vector of jackknifed averages, where element i is the average
  of all participants except participant i

# Examples
```julia
using JLD2

# Load LRP data from multiple participants
lrp_data = ErpData[]
for participant in 1:20
    data = load("participant_\$(participant)_lrp.jld2", "lrp")
    # If data is a vector, take the first condition, or specify which one
    lrp_result = data isa Vector ? data[1] : data
    push!(lrp_data, lrp_result)
end

# Create jackknife averages
jackknife_results = jackknife_average(lrp_data)

# Now jackknife_results[1] is the average of participants 2-20 (excluding 1)
# jackknife_results[2] is the average of participants 1,3-20 (excluding 2)
# etc.

# Save results
jldsave("jackknife_lrp.jld2"; data = jackknife_results)
```

# Multiple Conditions
```julia
# If you have multiple condition pairs, process each separately
lrp_data_cond1 = [load("participant_\$(i)_lrp.jld2", "lrp")[1] for i in 1:20]
lrp_data_cond2 = [load("participant_\$(i)_lrp.jld2", "lrp")[2] for i in 1:20]

jackknife_cond1 = jackknife_average(lrp_data_cond1)
jackknife_cond2 = jackknife_average(lrp_data_cond2)
```

# Notes
- Requires at least 2 participants
- All ERP/LRP data must have matching structure (same channels, time points, sample rate)
- The resulting data has the same format as the input (ErpData objects)
- Common workflow: Calculate LRP → Jackknife average → Statistical testing
"""
function jackknife_average(erps::Vector{ErpData}; condition_selection::Function = conditions())::Vector{ErpData}
    @info "Starting jackknife averaging"

    # Apply condition_selection first
    erps_filtered = erps[get_selected_conditions(erps, condition_selection)]

    # Validate inputs
    if (error_msg = _validate_jackknife_params(erps_filtered)) !== nothing
        @minimal_error_throw(error_msg)
    end

    # Create jackknife averages
    jackknife_results = _create_jackknife_averages(erps_filtered)

    @info "Jackknife averaging complete"
    return jackknife_results
end


"""
    jackknife_average(file_pattern::String;
                      input_dir::String = pwd(),
                      participant_selection::Function = participants(),
                      condition_selection::Function = conditions(),
                      output_dir::Union{String, Nothing} = nothing,
                      data_var::String = "lrp")

Batch process ERP/LRP files to create jackknife averages across participants.

This function loads data from multiple participant files, groups by condition,
and creates jackknife (leave-one-out) averages for each participant and condition.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "lrp", "erps_cleaned")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `data_var::String`: Deprecated parameter kept for backwards compatibility. Data is now loaded using `load_data()` which finds data by type.

# Examples
```julia
# Jackknife average all LRP files in current directory
jackknife_average("lrp")

# Process specific participants and conditions
jackknife_average("lrp",
                  input_dir = "/path/to/data",
                  participants = 1:20,
                  conditions = [1, 2])

# Process ERP data (data is automatically detected by type)
jackknife_average("erps_cleaned",
                  input_dir = "/path/to/data")

# Specify custom output directory
jackknife_average("lrp",
                  output_dir = "/path/to/output")
```

# Output
The function creates a new directory containing jackknifed data files:
- One file per participant
- Each file contains "data" variable with jackknifed average(s)
- If multiple conditions exist, each file contains a Vector{ErpData} with one element per condition
- Log file saved to output directory

# Typical Workflow
```julia
# 1. Calculate LRP for all participants
pairs = [(i, i+1) for i in 1:2:15]
lrp("erps_cleaned", pairs)

# 2. Create jackknife averages
jackknife_average("lrp")

# 3. Results are in jackknife_lrp/ directory
# Load and analyze:
using JLD2
participant_1_jackknife = load("jackknife_lrp/1_lrp.jld2", "data")
```
"""
function jackknife_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
    data_var::String = "lrp",
)
    # Setup logging
    log_file = "jackknife.log"
    setup_global_logging(log_file)

    try
        @info "Batch jackknife averaging started at $(now())"
        @log_call "jackknife_average"

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_jackknife_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        if length(files) < 2
            @minimal_warning "Need at least 2 participants for jackknife averaging, found $(length(files))"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Load and group data by condition (data_var parameter kept for backwards compatibility but not used)
        erps_by_condition, participant_ids =
            _load_and_group_for_jackknife(files, input_dir, condition_selection, data_var)

        if isempty(erps_by_condition)
            @minimal_warning "No valid data found in any files"
            return nothing
        end

        @info "Found conditions: $(sort(collect(keys(erps_by_condition))))"

        # Create jackknife averages for each condition
        jackknife_by_condition = Dict{Int,Vector{ErpData}}()

        for (cond_num, erps) in erps_by_condition
            @info "Creating jackknife averages for condition $cond_num (n=$(length(erps)) participants)"

            if length(erps) < 2
                @minimal_warning "Only $(length(erps)) participant(s) for condition $cond_num. Skipping jackknife."
                continue
            end

            jackknife_erps = _create_jackknife_averages(erps)
            jackknife_by_condition[cond_num] = jackknife_erps
        end

        if isempty(jackknife_by_condition)
            @minimal_warning "No jackknife averages created"
            return nothing
        end

        # Save jackknife data: one file per participant
        # Each file contains Vector{ErpData} if multiple conditions, or single ErpData if one condition
        n_participants = length(files)

        for (idx, (file, participant_id)) in enumerate(zip(files, participant_ids))
            # Collect jackknife data for this participant across all conditions
            participant_jackknife = ErpData[]

            for cond_num in sort(collect(keys(jackknife_by_condition)))
                jackknife_erps = jackknife_by_condition[cond_num]

                # The idx-th jackknife is the one excluding participant idx
                if idx <= length(jackknife_erps)
                    push!(participant_jackknife, jackknife_erps[idx])
                end
            end

            if !isempty(participant_jackknife)
                # Save with participant ID from filename
                output_file = file  # Keep original filename
                output_path = joinpath(output_dir, output_file)

                # If single condition, save as single ErpData, otherwise as Vector
                data_to_save = length(participant_jackknife) == 1 ? participant_jackknife[1] : participant_jackknife
                jldsave(output_path; data = data_to_save)

                @info "  Saved participant $participant_id: $output_file"
            end
        end

        @info "Jackknife averaging complete!"
        @info "Output saved to: $output_dir"

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
