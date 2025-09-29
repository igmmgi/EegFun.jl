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

function check_files_exist(conditions:: Int, filetype::String)
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
    files = filter(f -> occursin(Regex(files), f), readdir(directory))
    return [joinpath(directory, file) for file in files]
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
function find_file(filename::String, search_dir::String; 
                   recursive::Bool = true, extensions::Vector{String} = String[])
    
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


