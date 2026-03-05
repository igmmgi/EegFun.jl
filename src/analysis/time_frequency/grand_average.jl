"""
Grand averaging of time-frequency data across participants.
"""

function grand_average(tfs::Vector{TimeFreqData}; condition_selection::Function = conditions())
    isempty(tfs) && @minimal_error("Cannot create grand average from empty TF list")

    # Group by condition
    tfs_by_condition = group_by_condition(tfs)

    # Apply condition selection
    all_cond_nums = collect(keys(tfs_by_condition))
    selected_mask = condition_selection(1:length(all_cond_nums))
    selected_cond_nums = all_cond_nums[selected_mask]

    if isempty(selected_cond_nums)
        @minimal_warning "No conditions selected for TF grand averaging"
        return TimeFreqData[]
    end

    grand_averages = TimeFreqData[]

    for cond_num in sort(selected_cond_nums)
        condition_tfs = tfs_by_condition[cond_num]
        @info "Creating TF grand average for condition $cond_num (n=$(length(condition_tfs)) participants)"

        if length(condition_tfs) < 2
            @minimal_warning "Only $(length(condition_tfs)) participant(s) for condition $cond_num. Skipping."
            continue
        end

        # Validate consistent structure
        _have_same_structure(condition_tfs) || @minimal_error("TF data for condition $cond_num have inconsistent structure")

        # Get shared metadata
        ref_tf = condition_tfs[1]
        electrodes = channel_labels(ref_tf)
        frequencies = Float64.(sort(unique(ref_tf.data_power.freq)))
        time_points = Float64.(sort(unique(ref_tf.data_power.time)))

        # Extract 4D array and create grand average
        data_4d = _extract_tf_array(condition_tfs, electrodes, frequencies, time_points)
        grand_avg = _create_tf_grand_average(ref_tf, data_4d, electrodes, frequencies, time_points, cond_num)
        push!(grand_averages, grand_avg)

        @info "  ✓ Created TF grand average: $(length(frequencies)) freqs × $(length(time_points)) times"
    end

    return grand_averages
end
