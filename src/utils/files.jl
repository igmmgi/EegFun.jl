"""
    get_files(directory::String, files::String; recursive::Bool = false)

Get files from a directory matching a pattern.
"""
function get_files(directory::String, files::String; recursive::Bool = false)
    if recursive
        # Walk subdirectories (e.g., BIDS sub-XX/eeg/ structure)
        matching_files = String[]
        for (root, _, filenames) in walkdir(directory)
            for f in filenames
                if occursin(Regex(files), f)
                    push!(matching_files, joinpath(root, f))
                end
            end
        end
        return sort(matching_files, by = x -> (_natural_sort_key(basename(x)), x))
    else
        # Original flat directory search
        matching_files = filter(f -> occursin(Regex(files), f), readdir(directory))
        sorted_files = sort(matching_files, by = x -> (_natural_sort_key(x), x))
        return [joinpath(directory, file) for file in sorted_files]
    end
end

"""
    get_files(directory::String, files::Vector{String})

Get files from a directory matching a list of filenames.
"""
function get_files(directory::String, files::Vector{String})
    return [joinpath(directory, file) for file in files]
end

"""
    clean_pattern(pattern::String) -> String

Sanitize a file pattern by taking the substring before the first dot.
Useful for generating clean default output directory and file names.
"""
function clean_pattern(pattern::String)
    return string(first(split(pattern, '.')))
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
    get_file_extension(filename::String) -> String

Get the lowercase file extension including the dot (e.g., ".bdf", ".vhdr").
"""
get_file_extension(filename::String) = lowercase(splitext(filename)[2])
