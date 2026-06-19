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
    if !_have_same_structure(erps)
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
        @minimal_error("Need at least 2 participants for jackknife averaging, got $(n_participants)")
    end

    # Get metadata columns and EEG channels from first ERP
    first_erp = erps[1]
    metadata_cols = meta_labels(first_erp)
    eeg_channels = setdiff(propertynames(first_erp.data), metadata_cols)
    n_timepoints = nrow(first_erp.data)

    # Pre-allocate result vector
    jackknife_erps = ErpData[]

    @info "Creating jackknife averages for $n_participants participants"

    # Step 1: Pre-compute the total sum for each channel to enable O(N) jackknifing
    total_sums = Dict{Symbol,Vector{Float64}}()
    for ch in eeg_channels
        total_sums[ch] = zeros(Float64, n_timepoints)
        for erp in erps
            total_sums[ch] .+= erp.data[!, ch]
        end
    end

    # Pre-calculate total epochs
    total_epochs_all = sum(erp.n_epochs for erp in erps)

    # Step 2: Create jackknife average for each participant
    for excluded_idx = 1:n_participants
        excluded_erp = erps[excluded_idx]

        @info "  Participant $excluded_idx: averaging $(n_participants - 1) other participants"

        # Create a copy of the first ERP's data as the base (just to get the structure and time vector)
        jackknife_data = copy(first_erp.data)

        # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
        cols_to_remove = [:condition, :condition_name, :n_epochs]
        for col in cols_to_remove
            if hasproperty(jackknife_data, col)
                select!(jackknife_data, Not(col))
            end
        end

        # Calculate jackknife average: (Total - Excluded) / (N - 1)
        for ch in eeg_channels
            avg_col = (total_sums[ch] .- excluded_erp.data[!, ch]) ./ (n_participants - 1)
            jackknife_data[!, ch] = avg_col
        end

        # Calculate epochs for this jackknife average
        total_epochs = total_epochs_all - excluded_erp.n_epochs

        # Get condition info
        cond_name = first_erp.condition_name
        condition = first_erp.condition
        jackknife_cond_name = "$(cond_name)_jackknife_$(excluded_idx)"

        # Create ErpData object for this jackknife average
        jackknife_erp = ErpData(
            first_erp.file,
            condition,
            jackknife_cond_name,
            jackknife_data,
            copy(first_erp.layout),
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
function _load_and_group_for_jackknife(files::Vector{String}, input_dir::String, condition_selection::Function)
    # data_var parameter kept for backwards compatibility but not used - read_data() finds by type
    all_erps_by_condition = Dict{Int,Vector{ErpData}}()
    participant_ids = Int[]

    for (i, file) in enumerate(files)
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"

        # Extract participant ID from filename (assumes format like "1_pattern.jld2")
        m = match(r"^(\d+)_", file)
        participant_id = !isnothing(m) ? parse(Int, m.captures[1]) : i
        push!(participant_ids, participant_id)

        # Read data (using read_data which finds by type)
        data = read_data(input_path)

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
    error_msg = _validate_jackknife_params(erps_filtered)
    if !isnothing(error_msg)
        @minimal_error(error_msg)
    end

    # Create jackknife averages
    jackknife_results = _create_jackknife_averages(erps_filtered)

    @info "Jackknife averaging complete"
    return jackknife_results
end

function jackknife_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    # Setup logging
    log_file = "jackknife.log"
    setup_global_logging(log_file)

    try
        @info "Batch jackknife averaging started at $(now())"
        @log_call "jackknife_average"

        # Validation
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
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

        # Load and group data by condition
        erps_by_condition, participant_ids = _load_and_group_for_jackknife(files, input_dir, condition_selection)

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
