"""
Condition averaging for time-frequency data.
"""

#=============================================================================
    TF AVERAGE-SPECIFIC PROCESSING
=============================================================================#

"""
Create an average TF result by averaging power and phase across conditions.
"""
function _create_tf_average_wave(tfs::Vector{TimeFreqData}, conditions::Vector{Int}, avg_cond::Int)
    # Validate that all TFs have the same structure
    for i = 2:length(tfs)
        _have_same_structure(tfs[1], tfs[i]) || @minimal_error("TF data have inconsistent structure")
    end

    # Get EEG channels (exclude metadata columns like :time, :freq)
    metadata_cols = meta_labels(tfs[1])
    eeg_channels = setdiff(propertynames(tfs[1].data_power), metadata_cols)

    # Create copies of the first TF's DataFrames
    avg_power = copy(tfs[1].data_power)
    avg_phase = copy(tfs[1].data_phase)

    # Average EEG channels across conditions
    n = length(tfs)
    for ch in eeg_channels
        all_have_ch = all(hasproperty(tf.data_power, ch) for tf in tfs)
        if all_have_ch
            avg_power[!, ch] = sum(tf.data_power[!, ch] for tf in tfs) ./ n
            avg_phase[!, ch] = sum(tf.data_phase[!, ch] for tf in tfs) ./ n
        else
            @minimal_warning "Channel $ch not found in all conditions, keeping values from first condition"
        end
    end

    avg_condition_name = "average_$(join(conditions, "_"))"

    return TimeFreqData(
        tfs[1].file,
        avg_cond,
        avg_condition_name,
        avg_power,
        avg_phase,
        tfs[1].layout,
        tfs[1].sample_rate,
        tfs[1].method,
        tfs[1].baseline,
        tfs[1].analysis_info,
    )
end


#=============================================================================
    IN-MEMORY API FUNCTION
=============================================================================#

"""
    condition_average(data::Vector{TimeFreqData}, condition_groups::Vector{Vector{Int}})::Vector{TimeFreqData}

Create condition average time-frequency data.

For each group of condition indices, averages the power and phase across those conditions.
Both `data_power` and `data_phase` DataFrames are averaged consistently.

# Arguments
- `data::Vector{TimeFreqData}`: Vector of TimeFreqData, one per condition
- `condition_groups::Vector{Vector{Int}}`: Groups of condition indices to average (e.g., `[[1, 2], [3, 4]]`)

# Returns
- `Vector{TimeFreqData}`: New vector with one averaged TimeFreqData per group

# Examples
```julia
# Average conditions 1 and 2 together, and 3 and 4 together
avg_tf = EegFun.condition_average(tf_data, [[1, 2], [3, 4]])

# Average all four conditions into a single TF result
avg_tf = EegFun.condition_average(tf_data, [[1, 2, 3, 4]])

# Batch: process all participant files
EegFun.condition_average("tf_morlet", [[1, 2], [3, 4]], input_dir = "preprocessed")
```
"""
function condition_average(data::Vector{TimeFreqData}, condition_groups::Vector{Vector{Int}})::Vector{TimeFreqData}

    @info "Creating TF condition average waves..."

    average_waves = TimeFreqData[]

    for (group_idx, conditions) in enumerate(condition_groups)
        # Find TF data for each condition in the group
        group_tfs = TimeFreqData[]
        all_found = true

        for cond in conditions
            tf_found = nothing
            for tf in data
                if tf.condition == cond
                    tf_found = tf
                    break
                end
            end

            if isnothing(tf_found)
                @minimal_warning "Condition $cond not found. Skipping group $conditions."
                all_found = false
                break
            end
            push!(group_tfs, tf_found)
        end

        if !all_found
            continue
        end

        avg_wave = _create_tf_average_wave(group_tfs, conditions, group_idx)
        push!(average_waves, avg_wave)
    end

    if isempty(average_waves)
        @minimal_error("No valid condition groups found")
    end

    @info "Created $(length(average_waves)) TF average(s)"
    return average_waves
end
