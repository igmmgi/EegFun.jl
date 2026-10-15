"""
    subset_bad_data(data_path::String, threshold::Float64; 
                    subset_directory::String = "excluded")

Identify and move files from participants with low data retention to a separate directory.

Participants are considered "bad" if they have less than the threshold percentage
of data remaining in ANY condition. All files associated with bad participants
are moved to a subdirectory (default: "excluded").

The function searches for "epoch_summary.jld2" in the specified directory.

# Arguments
- `data_path::String`: Path to the directory containing epoch_summary.jld2 and preprocessed files
- `threshold::Float64`: Minimum percentage threshold (e.g., 75.0 means 75% retention required)

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
    (threshold < 0.0 || threshold > 100.0) && @minimal_error("threshold must be 0 < threshold < 100, got $threshold")
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
    println("Subsetting data: $(length(bad_participants))")
    println("   N remaining: $(length(unique(epoch_summary.file)) - length(bad_participants))")
    println("   N removed: $(length(bad_participants))")

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

    jldsave(epoch_summary_subset_path; data = epoch_summary_subset)
    jldsave(file_summary_subset_path; data = file_summary_subset)

    # Find and move all files for bad participants
    all_files = readdir(output_directory)
    for participant_id in bad_participants
        # Find all files that contain participant_id
        matching_files = filter(all_files) do filename
            occursin(participant_id, filename)
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
