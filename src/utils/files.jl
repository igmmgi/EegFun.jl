"""
    get_files(directory::String, files::String; recursive::Bool = false)
    get_files(directory::String, files::Vector{String}; recursive::Bool = false)

Return naturally sorted full paths in `directory` matching regex `files` (or explicit filenames `files`),
optionally searching subdirectories when `recursive = true`.
"""
function get_files(directory::String, files::String; recursive::Bool = false)
    rx = Regex(files)
    if recursive
        matching_files = String[]
        for (root, _, filenames) in walkdir(directory)
            for f in filenames
                occursin(rx, f) && push!(matching_files, joinpath(root, f))
            end
        end
        return sort!(matching_files, by = x -> (_natural_sort_key(basename(x)), x))
    end
    matching_files = filter(f -> occursin(rx, f), readdir(directory))
    sort!(matching_files, by = x -> (_natural_sort_key(x), x))
    return [joinpath(directory, file) for file in matching_files]
end

function get_files(directory::String, files::Vector{String}; recursive::Bool = false)
    recursive && return [something(find_file(file, directory; recursive = true), joinpath(directory, file)) for file in files]
    return [joinpath(directory, file) for file in files]
end

clean_pattern(pattern::String) = string(first(split(pattern, '.')))

"""
    find_file(filename::String, search_dir::String;
              recursive::Bool = true, extensions::Vector{String} = String[]) -> Union{String, Nothing}

Search `search_dir` (and subdirectories when `recursive = true`) for `filename` (or `filename .* extensions`),
returning the first matching full path or `nothing`.
"""
function find_file(filename::String, search_dir::String; recursive::Bool = true, extensions::Vector{String} = String[])
    (!isdir(search_dir) || (isempty(filename) && isempty(extensions))) && return nothing
    targets = isempty(extensions) ? (filename,) : [filename * ext for ext in extensions]

    for target in targets
        exact_path = joinpath(search_dir, target)
        isfile(exact_path) && return exact_path
    end

    if recursive
        for (root, _, files) in walkdir(search_dir)
            for file in files
                file in targets && return joinpath(root, file)
            end
        end
    end
    return nothing
end

function get_file_extension(filename::String)
    endswith(lowercase(filename), ".xdf.gz") && return ".xdf.gz"
    return lowercase(splitext(filename)[2])
end
