function _matches_participant_file(filename::AbstractString, participant_id)
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
    pat = Regex(left_boundary * escaped * "(?:[^A-Za-z0-9]|\$)")
    return occursin(pat, filename)
end

"""
    subset_bad_data(data_path::String, threshold::Real; subset_directory::String="excluded")

Move files for participants whose epoch retention falls below `threshold`% in any condition
(from `epoch_summary.jld2` in `data_path`) into `subset_directory`, and write filtered
`epoch_summary_subset.jld2` and `file_summary_subset.jld2` files.
"""
function subset_bad_data(data_path::String, threshold::Real; subset_directory::String = "excluded")
    (threshold < 0.0 || threshold > 100.0) && @minimal_error("threshold must be between 0 and 100, got $threshold")
    !isdir(data_path) && @minimal_error("data_path must be a directory: $data_path")

    epoch_summary_path = joinpath(data_path, "epoch_summary.jld2")
    !isfile(epoch_summary_path) && @minimal_error("epoch_summary.jld2 not found: $data_path")
    file_summary_path = joinpath(data_path, "file_summary.jld2")
    !isfile(file_summary_path) && @minimal_error("file_summary.jld2 not found: $data_path")

    epoch_summary = load(epoch_summary_path, "data")
    file_summary = load(file_summary_path, "data")

    output_directory = abspath(data_path)
    subset_dir_path = mkpath(joinpath(output_directory, subset_directory))

    bad_participants = unique(epoch_summary.file[epoch_summary.percentage .< threshold])
    n_remaining = length(unique(epoch_summary.file)) - length(bad_participants)
    @info "Subsetting data: $(length(bad_participants)) removed, $n_remaining remaining"

    bad_set = Set(bad_participants)
    epoch_summary_subset = epoch_summary[.!in.(epoch_summary.file, Ref(bad_set)), :]
    file_summary_subset = file_summary[.!in.(file_summary.file, Ref(bad_set)), :]

    jldsave(joinpath(output_directory, "epoch_summary_subset.jld2"); data = epoch_summary_subset, format_version = EEGFUN_FORMAT_VERSION)
    jldsave(joinpath(output_directory, "file_summary_subset.jld2"); data = file_summary_subset, format_version = EEGFUN_FORMAT_VERSION)

    all_files = readdir(output_directory)
    for participant_id in bad_participants
        for filename in all_files
            _matches_participant_file(filename, participant_id) || continue
            src_path = joinpath(output_directory, filename)
            isfile(src_path) && mv(src_path, joinpath(subset_dir_path, filename); force = true)
        end
    end
end
