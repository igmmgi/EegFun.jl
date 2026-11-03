# Check whether `filename` belongs to `participant_id` using exact token boundaries.
function _matches_participant_file(filename::AbstractString, participant_id)
    # Never move summary, config, or global pipeline log files
    if filename in (
        "epoch_summary.jld2",
        "epoch_summary_subset.jld2",
        "file_summary.jld2",
        "file_summary_subset.jld2",
        "config.toml",
        "neighbors.toml",
    ) || startswith(filename, "preprocess_log_")
        return false
    end

    pid_raw = strip(string(participant_id))
    isempty(pid_raw) && return false
    pid_stem = replace(pid_raw, r"\.(bdf|edf|vhdr|eeg|fif|xdf|set|mat|jld2|csv|txt)$"i => "")
    isempty(pid_stem) && return false

    escaped = replace(pid_stem, r"([.^$*+?()\[\]{}|\\])" => s"\\\1")
    left_boundary = isdigit(first(pid_stem)) ? "(?:^|[^0-9])" : "(?:^|[^A-Za-z0-9])"
    right_boundary = "(?:[^A-Za-z0-9]|\$)"
    pat = Regex(left_boundary * escaped * right_boundary)
    return occursin(pat, filename)
end

"""
    subset_bad_data(data_path::String, threshold::Real;
                    subset_directory::String = "excluded")

Identify and move files from participants with low data retention to a separate directory.

Participants are considered "bad" if they have less than the threshold percentage
of data remaining in ANY condition. All files associated with bad participants
are moved to a subdirectory (default: "excluded").

The function searches for "epoch_summary.jld2" in the specified directory.

# Arguments
- `data_path::String`: Path to the directory containing epoch_summary.jld2 and preprocessed files
- `threshold::Real`: Minimum percentage threshold (e.g., 75.0 means 75% retention required)

# Keyword Arguments
- `subset_directory::String`: Name of subdirectory for excluded participants (default: "excluded")

# Examples
```julia
# Move participants with < 75% data in any condition to "excluded" subdirectory
subset_bad_data("preprocessed_files", 75.0)

# Specify custom subset directory name
subset_bad_data("/path/to/preprocessed", 80.0, subset_directory="excluded")
```
"""
function subset_bad_data(data_path::String, threshold::Real; subset_directory::String = "excluded")

    # Validate inputs/outputs
    (threshold < 0.0 || threshold > 100.0) && @minimal_error("threshold must be between 0 and 100, got $threshold")
    !isdir(data_path) && @minimal_error("data_path must be a directory: $data_path")

    epoch_summary_path = joinpath(data_path, "epoch_summary.jld2")
    !isfile(epoch_summary_path) && @minimal_error("epoch_summary.jld2 not found: $data_path")

    # Load epoch summary (plain DataFrame, not EegFunData — can't use read_data)
    epoch_summary = load(epoch_summary_path, "data")

    # Use data_path as output directory
    output_directory = abspath(data_path)
    subset_dir_path = joinpath(output_directory, subset_directory)
    !isdir(subset_dir_path) && mkpath(subset_dir_path)

    # Find participants with any condition below threshold
    bad_participants = unique(epoch_summary.file[epoch_summary.percentage .< threshold])
    n_remaining = length(unique(epoch_summary.file)) - length(bad_participants)
    @info "Subsetting data: $(length(bad_participants)) removed, $n_remaining remaining"

    # Create subset summary files with only non-excluded participants
    # Filter epoch_summary to exclude bad participants
    epoch_summary_subset = epoch_summary[.!in.(epoch_summary.file, Ref(bad_participants)), :]

    # Load and filter file_summary_subset
    file_summary_path = joinpath(data_path, "file_summary.jld2")
    !isfile(file_summary_path) && @minimal_error("file_summary.jld2 not found: $data_path")

    file_summary = load(file_summary_path, "data")
    file_summary_subset = file_summary[.!in.(file_summary.file, Ref(bad_participants)), :]

    # Save subset summary files
    epoch_summary_subset_path = joinpath(output_directory, "epoch_summary_subset.jld2")
    file_summary_subset_path = joinpath(output_directory, "file_summary_subset.jld2")

    jldsave(epoch_summary_subset_path; data = epoch_summary_subset, format_version = EEGFUN_FORMAT_VERSION)
    jldsave(file_summary_subset_path; data = file_summary_subset, format_version = EEGFUN_FORMAT_VERSION)

    # Find and move all files for bad participants
    all_files = readdir(output_directory)
    for participant_id in bad_participants
        # Find all files matching exact participant_id boundary
        matching_files = filter(all_files) do filename
            _matches_participant_file(filename, participant_id)
        end
        # Move each matching file
        for filename in matching_files
            src_path = joinpath(output_directory, filename)
            dst_path = joinpath(subset_dir_path, filename)
            if isfile(src_path)
                mv(src_path, dst_path, force = true)
            end
        end
    end
end
