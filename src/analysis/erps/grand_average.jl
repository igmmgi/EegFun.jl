"""
Batch grand averaging of ERP data across participants.
"""

#=============================================================================
    GRANDAVERAGE-SPECIFIC HELPERS
=============================================================================#

"""Generate default output directory name for grand averaging."""
function _default_grand_average_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "grand_average_$(pattern)")
end

"""
Create a grand average by averaging ERP data across participants for a specific condition.
"""
function _create_grand_average(erps::Vector{ErpData}, cond_num::Int)
    if isempty(erps)
        @minimal_error("Cannot create grand average from empty ERP list")
    end

    # Validate that all ERPs have the same structure
    _have_same_structure(erps) || @minimal_error("ERPs have inconsistent structure")

    first_erp = erps[1]

    # Get metadata columns and find common EEG channels across all ERPs
    metadata_cols = meta_labels(first_erp)

    # Find intersection of channels across all ERPs (only average channels that exist in all)
    all_channel_sets = [setdiff(propertynames(erp.data), metadata_cols) for erp in erps]
    eeg_channels = collect(intersect(all_channel_sets...))

    # Create a copy of the first ERP's data as the base
    grand_avg_data = copy(first_erp.data)

    # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
    cols_to_remove = [:condition, :condition_name, :n_epochs]
    for col in cols_to_remove
        if hasproperty(grand_avg_data, col)
            select!(grand_avg_data, Not(col))
        end
    end

    n_points = nrow(grand_avg_data)
    n_erps = length(erps)

    # Pre-allocate summation buffer
    avg_buf = Vector{Float64}(undef, n_points)

    # Average EEG channels across participants
    for ch in eeg_channels
        fill!(avg_buf, 0.0)

        # Accumulate without creating matrix or intermediate arrays
        for erp in erps
            col_data = erp.data[!, ch]::Vector{Float64}
            @inbounds @simd for i = 1:n_points
                avg_buf[i] += col_data[i]
            end
        end

        # Divide by N
        @inbounds @simd for i = 1:n_points
            avg_buf[i] /= n_erps
        end

        grand_avg_data[!, ch] = copy(avg_buf)
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
        copy(first_erp.layout),
        first_erp.sample_rate,
        copy(first_erp.analysis_info),
        total_epochs,
    )
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


"""
    grand_average(erps::Vector{ErpData}; condition_selection = conditions()) -> Vector{ErpData}
    grand_average(tfs::Vector{TimeFreqData}; condition_selection = conditions()) -> Vector{TimeFreqData}
    grand_average(rsa_data_list::Vector{RsaData}; compute_noise_ceiling::Bool = true) -> RsaData
    grand_average(decoded_list::Vector{DecodedData}) -> DecodedData
    grand_average(file_pattern::String; input_dir::String = pwd(),
    participant_selection::Function = participants(), condition_selection::Function = conditions(),
    output_dir = nothing)

Create grand averages by averaging across participants.

For `ErpData` / `TimeFreqData`: groups by condition and averages across participants.
For `RsaData`: averages RDMs at each time point with optional noise ceiling.
For `DecodedData`: averages classification accuracy across participants.
The `file_pattern` method batch-processes JLD2 ERP files.

# Examples
```julia
grand_avgs = grand_average([erp1, erp2, erp3])
grand_avgs = grand_average(erps, condition_selection = conditions([1, 3]))
grand_avg_rsa = grand_average([rsa_p1, rsa_p2, rsa_p3])
grand_avg_decoding = grand_average([decoded_p1, decoded_p2])
grand_average("erps_unrejected")
```
"""
function grand_average(erps::Vector{ErpData}; condition_selection::Function = conditions())

    # Group ERPs by condition
    erps_by_condition = group_by_condition(erps)

    # Apply condition selection
    all_cond_nums = collect(keys(erps_by_condition))
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    filtered_groups = OrderedDict(num => erps_by_condition[num] for num in selected_cond_nums)

    if isempty(filtered_groups)
        @minimal_warning "No conditions selected or found for grand averaging"
        return ErpData[]
    end

    # Create grand averages
    return _create_all_grand_averages(filtered_groups)
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

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
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
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

        # Load all data — narrow to concrete type for dispatch
        all_data = _read_all_data_core(EegData, files, input_dir)

        if isempty(all_data)
            @minimal_warning "No valid data found in any files"
            return nothing
        end

        # Narrow abstract Vector{EegData} to concrete type for dispatch
        T = typeof(all_data[1])
        typed_data = T[x for x in all_data]

        # Dispatch to the appropriate Vector method
        grand_averages = grand_average(typed_data; condition_selection = condition_selection)

        if isempty(grand_averages)
            @minimal_warning "No grand averages created"
            return nothing
        end

        # Save
        output_file = "grand_average_$(file_pattern).jld2"
        output_path = joinpath(output_dir, output_file)
        jldsave(output_path; data = grand_averages)

        @info "Grand averaging complete! Created $(length(grand_averages)) grand averages"
        @info "Output saved to: $output_path"

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
