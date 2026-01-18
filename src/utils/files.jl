"""
    check_files_exist(conditions::Union{Vector{Int}, Int} filetype::String) -> Bool

Check if files exist for all given conditions with specified filetype.

# Arguments
- `conditions::Vector{String}`: List of condition names
- `filetype::String`: Type of file to check

# Returns
- `Bool`: true if all files exist, false otherwise
"""
function check_files_exist(conditions::Vector{Int}, filetype::String)
    all_files_exist = true
    for condition in conditions
        fname = "$(condition)_$(filetype).jld2"
        if !isfile(fname)
            @minimal_warning "File not found: $(fname)"
            all_files_exist = false
        end
    end
    return all_files_exist
end

function check_files_exist(conditions::Int, filetype::String)
    return check_files_exist([conditions], filetype)
end


function check_files_exist(files::Vector{String})
    all_files_exist = true
    for fname in files
        if !isfile(fname)
            @minimal_warning "File not found: $(fname)"
            all_files_exist = false
        end
    end
    return all_files_exist
end

function get_files(directory::String, files::String)
    # replace common wildcard with regex syntax
    matching_files = Base.filter(f -> occursin(Regex(files), f), readdir(directory))
    # Natural order: "file_10" comes after "file_3"
    sorted_files = sort(matching_files, by = x -> (natural_sort_key(x), x))
    return [joinpath(directory, file) for file in sorted_files]
end

function get_files(directory::String, files::Vector{String})
    return [joinpath(directory, file) for file in files]
end

"""
    find_file(filename::String, search_dir::String; 
              recursive::Bool = true, extensions::Vector{String} = []) -> Union{String, Nothing}

Find a file by searching through a directory, optionally recursively.

# Arguments
- `filename::String`: Name of the file to find (e.g., "biosemi64.csv", "config.toml")
- `search_dir::String`: Directory to search in
- `recursive::Bool`: Whether to search subdirectories recursively (default: true)
- `extensions::Vector{String}`: Optional file extensions to match (e.g., [".csv", ".toml"])

# Returns
- `String` or `nothing`: Full path to the file if found, `nothing` otherwise

# Examples
```julia
# Find a layout file in the layouts directory
layout_file = find_file("biosemi64.csv", "data/layouts")

# Find any CSV file with that name
csv_file = find_file("biosemi64", "data/layouts", extensions = [".csv"])

# Non-recursive search (only direct children)
direct_file = find_file("README.txt", "data/layouts", recursive = false)
```
"""
function find_file(filename::String, search_dir::String; recursive::Bool = true, extensions::Vector{String} = String[])

    # Check if directory exists
    if !isdir(search_dir)
        return nothing
    end

    # Add extensions to filename if specified
    if isempty(extensions)
        filenames_to_find = [filename]
    else
        filenames_to_find = [filename * ext for ext in extensions]
    end

    # Try exact match in the directory first
    for target_filename in filenames_to_find
        exact_path = joinpath(search_dir, target_filename)
        if isfile(exact_path)
            return exact_path
        end
    end

    # If recursive search is enabled, search subdirectories
    if recursive
        for (root, dirs, files) in walkdir(search_dir)
            for file in files
                if file in filenames_to_find
                    return joinpath(root, file)
                end
            end
        end
    end

    return nothing
end

"""
    _filter_files(files::Vector{String}; include::Union{Vector{Int}, Int, Nothing} = nothing, 
                  exclude::Union{Vector{Int}, Int, Nothing} = nothing) -> Vector{String}

Filter files by participant number extracted from filename.

# Arguments
- `files::Vector{String}`: List of filenames to filter
- `include::Union{Vector{Int}, Int, Nothing}`: Participant number(s) to include (default: nothing = all)
- `exclude::Union{Vector{Int}, Int, Nothing}`: Participant number(s) to exclude (default: nothing = none)

# Returns
- `Vector{String}`: Filtered list of files matching the criteria

# Note
Assumes filename format like "Flank_C_3_epochs_cleaned.jld2" where participant number is the first numeric part.
"""
function _filter_files(
    files::Vector{String};
    include::Union{Vector{Int},Int,Nothing} = nothing,
    exclude::Union{Vector{Int},Int,Nothing} = nothing,
)

    # Convert single values to vectors
    include_nums = include isa Int ? [include] : include
    exclude_nums = exclude isa Int ? [exclude] : exclude

    return Base.filter(files) do file
        # Extract participant number from filename (assuming format like "Flank_C_3_epochs_cleaned.jld2")
        parts = split(file, "_")
        file_participant = nothing

        for part in parts
            if !isempty(part) && isdigit(part[1])  # First character is a digit
                file_participant = parse(Int, part)
                break
            end
        end

        # If no participant number found, include the file unless explicitly excluded
        if file_participant === nothing
            return exclude_nums === nothing || !(nothing in exclude_nums)
        end

        # Apply include filter
        if include_nums !== nothing && !(file_participant in include_nums)
            return false
        end

        # Apply exclude filter
        if exclude_nums !== nothing && file_participant in exclude_nums
            return false
        end

        return true
    end
end

"""
    check_files_exist(subjects::Union{Vector{Int}, Int}, conditions::Union{Vector{Int}, Int},, filetype::String) -> Bool

Check if files exist for all combinations of subjects and conditions.

# Arguments
- `subjects::Vector{String}`: List of subject identifiers
- `conditions::Vector{String}`: List of condition names
- `filetype::String`: Type of file to check

# Returns
- `Bool`: true if all files exist, false otherwise
"""
function check_files_exist(subjects::Union{Vector{Int},Int}, conditions::Union{Vector{Int},Int}, filetype::String)
    all_files_exist = true
    for subject in subjects
        for condition in conditions
            fname = "$(subject)_$(condition)_$(filetype).jld2"
            if !isfile(fname)
                @minimal_warning "File not found: $(fname)"
                all_files_exist = false
            end
        end
    end
    return all_files_exist
end

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
function subset_bad_data(
    data_path::String,
    threshold::Float64;
    subset_directory::String = "excluded",
)

    # Validate inputs/outputs
    (threshold < 0.0 || threshold > 100.0) && @minimal_error_throw("threshold must be 0 < threshold < 100, got $threshold")
    !isdir(data_path) && @minimal_error_throw("data_path must be a directory: $data_path")

    epoch_summary_path = joinpath(data_path, "epoch_summary.jld2")
    !isfile(epoch_summary_path) && @minimal_error_throw("epoch_summary.jld2 not found: $data_path")
    
    # Load epoch summary and check required columns
    epoch_summary = load_data(epoch_summary_path)
    missing_cols = [col for col in [:file, :percentage] if !(col in propertynames(epoch_summary))]
    !isempty(missing_cols) && @minimal_error_throw("epoch_summary missing required columns: $(missing_cols)")

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
    !isfile(file_summary_path) && @minimal_error_throw("file_summary.jld2 not found: $data_path")
    
    file_summary = load_data(file_summary_path)
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
        matching_files = Base.filter(all_files) do filename
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
