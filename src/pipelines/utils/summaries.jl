"""
    _epoch_and_file_summary(all_epoch_counts::Vector{DataFrame})

Combine epoch counts from all files and create averaged version.

# Arguments
- `all_epoch_counts::Vector{DataFrame}`: Vector of DataFrames with epoch counts per file

# Returns
- `Tuple{DataFrame, DataFrame}`: (combined_counts, avg_counts)
  - `combined_counts`: All epoch counts concatenated
  - `avg_counts`: Averaged percentage per file (grouped by file)
"""
function _epoch_and_file_summary(all_epoch_counts::Vector{DataFrame})
    epoch_summary = vcat(all_epoch_counts...)
    file_summary = combine(groupby(epoch_summary, [:file]), :percentage => mean => :percentage)
    file_summary.percentage = round.(file_summary.percentage; digits = 1)
    return epoch_summary, file_summary
end

"""
    _merge_summaries(new_epoch_summary::DataFrame, new_file_summary::DataFrame, output_directory::String)

Merge new summaries with existing summaries, replacing data for files that already exist.

# Arguments
- `new_epoch_summary::DataFrame`: New epoch summary data
- `new_file_summary::DataFrame`: New file summary data
- `output_directory::String`: Directory where summary files are stored

# Returns
- `Tuple{DataFrame, DataFrame}`: (merged_epoch_summary, merged_file_summary)
"""
function _merge_summaries(new_epoch_summary::DataFrame, new_file_summary::DataFrame, output_directory::String)
    epoch_path = joinpath(output_directory, "epoch_summary.jld2")
    file_path = joinpath(output_directory, "file_summary.jld2")

    # Get list of files in new data (for replacement)
    new_files = unique(new_epoch_summary.file)

    # Load and merge epoch summary
    if isfile(epoch_path)
        existing_epoch = JLD2.load(epoch_path, "data")
        # Remove rows for files that are being updated
        existing_epoch = existing_epoch[.!in.(existing_epoch.file, Ref(new_files)), :]
        merged_epoch = vcat(existing_epoch, new_epoch_summary)
    else
        merged_epoch = new_epoch_summary
    end

    # Load and merge file summary
    if isfile(file_path)
        existing_file = JLD2.load(file_path, "data")
        # Remove rows for files that are being updated
        existing_file = existing_file[.!in.(existing_file.file, Ref(new_files)), :]
        merged_file = vcat(existing_file, new_file_summary)
    else
        merged_file = new_file_summary
    end

    # Sort by filename using natural sorting (file_1, file_2, ..., file_10)
    merged_epoch = sort(merged_epoch, :file, by = _natural_sort_key)
    merged_file = sort(merged_file, :file, by = _natural_sort_key)

    return merged_epoch, merged_file
end

"""
    compare_rejections(rejection_step1::Vector{EpochRejectionInfo}, rejection_step2::Vector{EpochRejectionInfo})

Compare rejection results before and after repair to show repair effectiveness.

# Arguments
- `rejection_step1::Vector{EpochRejectionInfo}`: Rejection info before repair
- `rejection_step2::Vector{EpochRejectionInfo}`: Rejection info after repair

# Returns
- `DataFrame`: Comparison table with detailed statistics per condition showing:
  - Original number of epochs
  - Rejected before repair (step1)
  - Rejected after repair (step2)
  - Successfully repaired epochs
  - Percentage of step1 rejections that were repaired
  - Final percentage of epochs kept

# Examples
```julia
rejection_step1 = detect_bad_epochs_automatic(epochs, name="step1")
repair_artifacts!(epochs, rejection_step1)
rejection_step2 = detect_bad_epochs_automatic(epochs, name="step2")
comparison = compare_rejections(rejection_step1, rejection_step2)
```
"""
function compare_rejections(rejection1::Vector{EpochRejectionInfo}, rejection2::Vector{EpochRejectionInfo})
    comparison_data = []

    for (info1, info2) in zip(rejection1, rejection2)
        # Get unique epoch indices rejected in each step
        epochs_rejected_step1 = unique_epochs(info1)
        epochs_rejected_step2 = unique_epochs(info2)

        # Calculate statistics
        n_original = info1.info.n
        n_rejected_step1 = length(epochs_rejected_step1)
        n_rejected_step2 = length(epochs_rejected_step2)

        # Epochs that were rejected in step1 but NOT in step2 (successfully repaired)
        epochs_repaired = setdiff(epochs_rejected_step1, epochs_rejected_step2)
        n_repaired = length(epochs_repaired)

        # Epochs rejected in BOTH steps (repair didn't help)
        epochs_still_bad = intersect(epochs_rejected_step1, epochs_rejected_step2)
        n_still_bad = length(epochs_still_bad)

        # New rejections in step2 (epochs that weren't rejected in step1)
        # This can happen if repair introduces artifacts or if thresholds are slightly different
        epochs_new_rejections = setdiff(epochs_rejected_step2, epochs_rejected_step1)
        n_new_rejections = length(epochs_new_rejections)

        # Verify: n_rejected_step2 should equal n_still_bad + n_new_rejections
        # Verify: n_rejected_step1 should equal n_repaired + n_still_bad
        @assert n_rejected_step2 == n_still_bad + n_new_rejections "Logic error: n_rejected_step2 ($n_rejected_step2) != n_still_bad ($n_still_bad) + n_new_rejections ($n_new_rejections)"
        @assert n_rejected_step1 == n_repaired + n_still_bad "Logic error: n_rejected_step1 ($n_rejected_step1) != n_repaired ($n_repaired) + n_still_bad ($n_still_bad)"

        # Calculate key percentages
        pct_repaired = n_rejected_step1 > 0 ? round(100 * n_repaired / n_rejected_step1, digits = 1) : 0.0
        pct_final_kept = round(100 * (n_original - n_rejected_step2) / n_original, digits = 1)

        push!(
            comparison_data,
            (
                condition = info1.info.number,
                condition_name = info1.info.name,
                n_original = n_original,
                rejected_before_repair = n_rejected_step1,
                rejected_after_repair = n_rejected_step2,
                successfully_repaired = n_repaired,
                pct_repaired = pct_repaired,
                pct_final_kept = pct_final_kept,
            ),
        )
    end

    return DataFrame(comparison_data)
end

"""
    summarize_electrode_repairs(continuous_repairs::Vector{ContinuousRepairInfo})

Summarize electrode repairs from a vector of ContinuousRepairInfo objects (continuous level only).
Each ContinuousRepairInfo is treated as coming from a separate participant.

# Arguments
- `continuous_repairs::Vector{ContinuousRepairInfo}`: Vector of ContinuousRepairInfo objects

# Returns
- `DataFrame`: Summary table with columns:
  - `electrode::Symbol`: Electrode name
  - `n_participants::Int`: Number of participants where this electrode was repaired at the continuous level
"""
function summarize_electrode_repairs(continuous_repairs::Vector{ContinuousRepairInfo})
    electrode_participants = Dict{Symbol,Set{Int}}()

    for (participant_idx, repair_info) in enumerate(continuous_repairs)
        # Count continuous repairs (1 per participant)
        for ch in repair_info.repaired
            if !haskey(electrode_participants, ch)
                electrode_participants[ch] = Set{Int}()
            end
            push!(electrode_participants[ch], participant_idx)
        end
    end

    # Create DataFrame sorted by number of participants (descending)
    summary_data = [
        (electrode = ch, n_participants = length(participants)) for
        (ch, participants) in sort(collect(electrode_participants), by = x -> length(x[2]), rev = true)
    ]

    return DataFrame(summary_data)
end
function summarize_electrode_repairs(file_pattern::String; input_dir::String = pwd())
    files = sort(_find_batch_files(file_pattern, input_dir), by = _natural_sort_key)

    empty_return =
        (DataFrame(file = String[], repaired = String[], skipped = String[]), DataFrame(electrode = Symbol[], n_participants = Int[]))
    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return empty_return
    end

    per_file_data = []
    continuous_repairs = ContinuousRepairInfo[]

    for file in files
        file_path = joinpath(input_dir, file)
        isfile(file_path) || continue
        artifact_info = try
            read_data(file_path)
        catch e
            @minimal_warning "Failed to load artifact info from $file_path: $(sprint(showerror, e))"
            continue
        end
        base_filename = replace(file, "_artifact_info.jld2" => "")
        file_repaired = Symbol[]
        file_skipped  = Symbol[]
        if !isnothing(artifact_info)
            for repair_info in artifact_info.continuous_repairs
                append!(file_repaired, repair_info.repaired)
                append!(file_skipped, repair_info.skipped)
                if !isempty(repair_info.repaired)
                    push!(continuous_repairs, repair_info)
                end
            end
        end
        push!(
            per_file_data,
            (
                file     = base_filename,
                repaired = isempty(file_repaired) ? "none" : join(sort(string.(unique(file_repaired))), ", "),
                skipped  = isempty(file_skipped) ? "" : join(sort(string.(unique(file_skipped))), ", "),
            ),
        )
    end

    per_file_df  = DataFrame(per_file_data)
    aggregate_df = summarize_electrode_repairs(continuous_repairs)
    return per_file_df, aggregate_df
end





"""
    summarize_ica_components(ica_components::Vector{ArtifactComponents})

Summarize ICA components removed across all participants.

# Arguments
- `ica_components::Vector{ArtifactComponents}`: Vector of ArtifactComponents objects

# Returns
- `DataFrame`: Summary table with columns:
  - `file::String`: Filename (empty string when called with Vector{ArtifactComponents})
  - `total_components::Int`: Total number of ICA components removed
  - `vEOG::Int`: Number of vEOG components
  - `hEOG::Int`: Number of hEOG components
  - `ECG::Int`: Number of ECG components
  - `line_noise::Int`: Number of line noise components
  - `channel_noise::Int`: Number of channel noise components
- `DataFrame`: Average summary with columns:
  - `component_type::String`: Type of component
  - `avg_per_participant::Float64`: Average number of components per participant
"""
function summarize_ica_components(ica_components::Vector{ArtifactComponents})
    if isempty(ica_components)
        return DataFrame(
            file = String[],
            total_components = Int[],
            vEOG = Int[],
            hEOG = Int[],
            ECG = Int[],
            line_noise = Int[],
            channel_noise = Int[],
        ),
        DataFrame(component_type = String[], avg_per_participant = Float64[])
    end

    # Per-file summary
    per_file_data = []
    for component_artifacts in ica_components
        vEOG_count = length(get(component_artifacts.eog, :vEOG, []))
        hEOG_count = length(get(component_artifacts.eog, :hEOG, []))
        ecg_count = length(component_artifacts.ecg)
        line_noise_count = length(component_artifacts.line_noise)
        channel_noise_count = length(component_artifacts.channel_noise)
        total_count = vEOG_count + hEOG_count + ecg_count + line_noise_count + channel_noise_count

        push!(
            per_file_data,
            (
                file = "",  # No filename when called with Vector{ArtifactComponents}
                total_components = total_count,
                vEOG = vEOG_count,
                hEOG = hEOG_count,
                ECG = ecg_count,
                line_noise = line_noise_count,
                channel_noise = channel_noise_count,
            ),
        )
    end

    per_file_df = DataFrame(per_file_data)

    # Average summary
    avg_data = [
        (component_type = "vEOG", avg_per_participant = round(mean(per_file_df.vEOG), digits = 2)),
        (component_type = "hEOG", avg_per_participant = round(mean(per_file_df.hEOG), digits = 2)),
        (component_type = "ECG", avg_per_participant = round(mean(per_file_df.ECG), digits = 2)),
        (component_type = "line_noise", avg_per_participant = round(mean(per_file_df.line_noise), digits = 2)),
        (component_type = "channel_noise", avg_per_participant = round(mean(per_file_df.channel_noise), digits = 2)),
        (component_type = "total", avg_per_participant = round(mean(per_file_df.total_components), digits = 2)),
    ]

    avg_df = DataFrame(avg_data)

    return per_file_df, avg_df
end
function summarize_ica_components(file_pattern::String; input_dir::String = pwd())
    # Find matching files using _find_batch_files
    files = _find_batch_files(file_pattern, input_dir)

    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return DataFrame(
            file = String[],
            total_components = Int[],
            vEOG = Int[],
            hEOG = Int[],
            ECG = Int[],
            line_noise = Int[],
            channel_noise = Int[],
        ),
        DataFrame(component_type = String[], avg_per_participant = Float64[])
    end

    # Load all ArtifactInfo objects and extract ArtifactComponents with filenames
    ica_components = ArtifactComponents[]
    filenames = String[]

    for file in files
        file_path = joinpath(input_dir, file)
        isfile(file_path) || continue
        artifact_info = try
            read_data(file_path)
        catch e
            @minimal_warning "Failed to load artifact info from $file_path: $(sprint(showerror, e))"
            continue
        end
        if !isnothing(artifact_info) && !isnothing(artifact_info.ica_components)
            push!(ica_components, artifact_info.ica_components)
            push!(filenames, file)
        end
    end

    # Get summary from vector version
    per_file_df, avg_df = summarize_ica_components(ica_components)

    # Update filenames in per_file_df: extract base filename
    for (i, filename) in enumerate(filenames)
        # Extract base filename (remove "_artifact_info.jld2" suffix)
        base_filename = replace(filename, "_artifact_info.jld2" => "")
        per_file_df.file[i] = base_filename
    end

    # Sort by filename using natural sorting (handles numeric parts correctly: Flank_C_3 before Flank_C_12)
    sort!(per_file_df, :file, by = _natural_sort_key)

    return per_file_df, avg_df
end
