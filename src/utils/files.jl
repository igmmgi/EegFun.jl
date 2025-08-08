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


