"""
Condition difference for time-frequency data.
"""

# === TF DIFFERENCE-SPECIFIC PROCESSING ===

"""
Create a difference TF result by subtracting tf2 from tf1 (power and phase).
"""
function _create_tf_difference_wave(tf1::TimeFreqData, tf2::TimeFreqData, cond1::Int, cond2::Int, diff_cond::Int)
    # Validate that both TFs have the same structure
    _have_same_structure(tf1, tf2) || @minimal_error("TF data have inconsistent structure")

    # Get EEG channels (exclude metadata columns like :time, :freq)
    metadata_cols = meta_labels(tf1)
    eeg_channels = setdiff(propertynames(tf1.data_power), metadata_cols)

    # Create copies of tf1's DataFrames for the difference
    diff_power = copy(tf1.data_power)
    diff_phase = copy(tf1.data_phase)

    # Subtract EEG channels
    for ch in eeg_channels
        if hasproperty(tf2.data_power, ch)
            diff_power[!, ch] = tf1.data_power[!, ch] .- tf2.data_power[!, ch]
            diff_phase[!, ch] = tf1.data_phase[!, ch] .- tf2.data_phase[!, ch]
        else
            @minimal_warning "Channel $ch not found in condition $cond2, keeping original values"
        end
    end

    diff_condition_name = "difference_$(cond1)_$(cond2)"

    return TimeFreqData(
        tf1.file,
        diff_cond,
        diff_condition_name,
        diff_power,
        diff_phase,
        copy(tf1.layout),
        tf1.sample_rate,
        tf1.method,
        tf1.baseline,
        copy(tf1.analysis_info),
    )
end


# === IN-MEMORY API FUNCTION ===

function condition_difference(
    data::Vector{TimeFreqData},
    condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}},
)::Vector{TimeFreqData}

    @info "Creating TF condition difference waves..."

    difference_waves = TimeFreqData[]

    for (pair_idx, pair) in enumerate(condition_pairs)
        cond1, cond2 = pair[1], pair[2]

        # Find TF data for each condition
        tf1 = nothing
        tf2 = nothing
        for tf in data
            cond_num = tf.condition
            if cond_num == cond1
                tf1 = tf
            elseif cond_num == cond2
                tf2 = tf
            end
        end

        if isnothing(tf1)
            @minimal_warning "Condition $cond1 not found. Skipping pair ($cond1, $cond2)."
            continue
        end
        if isnothing(tf2)
            @minimal_warning "Condition $cond2 not found. Skipping pair ($cond1, $cond2)."
            continue
        end

        diff_wave = _create_tf_difference_wave(tf1, tf2, cond1, cond2, pair_idx)
        push!(difference_waves, diff_wave)
    end

    if isempty(difference_waves)
        @minimal_error("No valid condition pairs found")
    end

    @info "Created $(length(difference_waves)) TF difference wave(s)"
    return difference_waves
end
